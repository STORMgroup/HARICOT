//! Per read, the raw-signal sample span of every RETAINED event: the coordinate space
//! bins, anchors and paths live in.
//!
//! A retained event is a run of alike-quantising raw events (`emer::retained_quantized`),
//! so retained event `r` covers raw events `[event_pos[r], event_pos[r+1])`. Its span runs
//! from the first of those events' start sample to the last one's end, in indices into the
//! pod5 signal. That is the join key to a basecaller move table (`truth/event_ref.py`),
//! which turns this into an exact event -> read base -> reference map, instead of a
//! uniform `events / rlen` rate.
//!
//! Output, little-endian, one record per read in pod5 order:
//!   `u32 name_len, name bytes, u32 n_samples, u32 n, n x u32 start, n x u32 end`
//!
//! Usage: `haricot-truth-events READS.pod5 OUT.bin [--limit N]`
use std::io::Write;

use anyhow::{Context, Result, ensure};
use haricot::emer;
use haricot::params::Params;
use haricot::signal::{self, RawRead};
use rayon::prelude::*;

struct Rec {
    name: String,
    n_samples: u32,
    start: Vec<u32>,
    end: Vec<u32>,
}

fn spans(r: &RawRead, params: &Params, check: bool) -> Result<Rec> {
    let ev = signal::signal_to_event_spans(&r.signal_pa, params);
    let means: Vec<f32> = ev.iter().map(|e| e.mean).collect();
    if check {
        // The pipeline's own path, so a drift in the span code cannot go unnoticed.
        ensure!(
            means == signal::signal_to_events(&r.signal_pa, params),
            "{}: span means differ from signal_to_events",
            r.id
        );
    }
    let retained = emer::retained_quantized(&means, params);
    let mut start = Vec::with_capacity(retained.len());
    let mut end = Vec::with_capacity(retained.len());
    for (k, q) in retained.iter().enumerate() {
        let first = q.event_pos as usize;
        let last = retained
            .get(k + 1)
            .map_or(ev.len(), |n| n.event_pos as usize)
            - 1;
        start.push(ev[first].start);
        end.push(ev[last].end);
    }
    Ok(Rec {
        name: r.id.clone(),
        n_samples: r.signal_pa.len() as u32,
        start,
        end,
    })
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let usage = "usage: haricot-truth-events READS.pod5 OUT.bin [--limit N]";
    let pod5 = args.get(1).context(usage)?;
    let out = args.get(2).context(usage)?;
    let limit: Option<usize> = args
        .iter()
        .position(|a| a == "--limit")
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok());
    let params = Params::r10_4_1();
    params.validate().map_err(anyhow::Error::msg)?;

    let mut w = std::io::BufWriter::new(
        std::fs::File::create(out).with_context(|| format!("creating {out}"))?,
    );
    let (mut n_reads, mut n_ev, mut chunk) = (0usize, 0usize, Vec::<RawRead>::new());
    let mut flush =
        |chunk: &mut Vec<RawRead>, n_reads: &mut usize, n_ev: &mut usize| -> Result<()> {
            let base = *n_reads;
            let recs: Vec<Rec> = chunk
                .par_iter()
                .enumerate()
                .map(|(i, r)| spans(r, &params, base + i < 200))
                .collect::<Result<_>>()?;
            for r in recs {
                w.write_all(&(r.name.len() as u32).to_le_bytes())?;
                w.write_all(r.name.as_bytes())?;
                w.write_all(&r.n_samples.to_le_bytes())?;
                w.write_all(&(r.start.len() as u32).to_le_bytes())?;
                for v in r.start.iter().chain(&r.end) {
                    w.write_all(&v.to_le_bytes())?;
                }
                *n_ev += r.start.len();
                *n_reads += 1;
            }
            chunk.clear();
            Ok(())
        };
    signal::for_each_raw_read(std::path::Path::new(pod5), limit, |r| {
        chunk.push(r);
        if chunk.len() >= 512 {
            flush(&mut chunk, &mut n_reads, &mut n_ev)?;
        }
        Ok(())
    })?;
    flush(&mut chunk, &mut n_reads, &mut n_ev)?;
    w.flush()?;
    eprintln!("{n_reads} reads, {n_ev} retained events -> {out}");
    Ok(())
}

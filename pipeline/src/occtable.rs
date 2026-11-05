//! The seed-occurrence table: how often each seed hash occurs over every read of the input, and the two
//! dataset statistics the stage-1 features derive from it, written once by the stage-1 features step
//! (`haricot-pair-features --save-occ`) so the contig-end rescue (`haricot-rescue rescue --occ`) and the rescue pool's
//! features (`haricot-pair-features --only .. --occ`) need not segment every read again: they segment only the
//! reads they score (the end reads, the pool's reads). About 12 bytes per distinct hash (~11 MB on a 30x
//! bacterial set), against ~1.3 GB for every read's levels and seeds.
//!
//! Exact: the table is [`SeedCounts::sorted_entries`] of the counts the features step computed over every
//! read's seeds, read back by [`SeedCounts::from_entries`] (the same answer for every hash); the two
//! statistics are written as `f64` bits.
//!
//! Format (little-endian): `RWOCC001`, the median occurrence over seed instances and the mean occurrence
//! over distinct hashes (`f64` each, as haricot-pair-features computes them), the entry count (`u64`), then per
//! entry the hash (`u64`) and its count (`u32`), sorted by hash.
#![allow(dead_code)] // each program that includes this file uses its own part

use std::io::{BufReader, BufWriter, Read, Write};

use anyhow::{Context, Result, bail, ensure};
use haricot::overlap::SeedCounts;

const MAGIC: &[u8; 8] = b"RWOCC001";

/// The table and the two statistics.
pub struct OccTable {
    pub counts: SeedCounts,
    /// Median occurrence over seed instances (`haricot-pair-features`' `occ_med`, after its `max(1)`).
    pub median: f64,
    /// Mean occurrence over distinct hashes (`haricot-pair-features`' `occ_mean`).
    pub mean: f64,
}

/// Write `counts` and the statistics to `path`.
pub fn save(path: &str, counts: &SeedCounts, median: f64, mean: f64) -> Result<()> {
    let f = std::fs::File::create(path).with_context(|| format!("creating {path}"))?;
    let mut w = BufWriter::with_capacity(1 << 20, f);
    let e = counts.sorted_entries();
    w.write_all(MAGIC)?;
    w.write_all(&median.to_le_bytes())?;
    w.write_all(&mean.to_le_bytes())?;
    w.write_all(&(e.len() as u64).to_le_bytes())?;
    for (h, c) in e {
        w.write_all(&h.to_le_bytes())?;
        w.write_all(&c.to_le_bytes())?;
    }
    w.flush()?;
    Ok(())
}

fn take<const N: usize>(r: &mut impl Read) -> Result<[u8; N]> {
    let mut b = [0u8; N];
    r.read_exact(&mut b)?;
    Ok(b)
}

/// The table at `path`.
pub fn load(path: &str) -> Result<OccTable> {
    let f = std::fs::File::open(path).with_context(|| format!("opening {path}"))?;
    let mut r = BufReader::with_capacity(1 << 20, f);
    if &take::<8>(&mut r)? != MAGIC {
        bail!("{path}: not a seed-occurrence table");
    }
    let median = f64::from_le_bytes(take(&mut r)?);
    let mean = f64::from_le_bytes(take(&mut r)?);
    let n = u64::from_le_bytes(take(&mut r)?) as usize;
    let mut entries = Vec::with_capacity(n);
    for _ in 0..n {
        let h = u64::from_le_bytes(take(&mut r)?);
        let c = u32::from_le_bytes(take(&mut r)?);
        entries.push((h, c));
    }
    ensure!(
        r.read(&mut [0u8; 1])? == 0,
        "{path}: trailing bytes after the last entry"
    );
    Ok(OccTable {
        counts: SeedCounts::from_entries(entries),
        median,
        mean,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_table_reads_back_exactly() {
        let lists: Vec<Vec<u64>> = vec![vec![5, 9, 5, 1 << 40], vec![9, 9, 7], vec![]];
        let c = SeedCounts::count(&lists, |l| l.len(), |l, i| l[i]);
        let dir = std::env::temp_dir().join(format!("rwocc_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("occ.bin");
        let p = p.to_str().unwrap();
        save(p, &c, 2.0, 1.75).unwrap();
        let t = load(p).unwrap();
        assert_eq!(
            (t.median.to_bits(), t.mean.to_bits()),
            (2.0f64.to_bits(), 1.75f64.to_bits())
        );
        for h in [1u64, 5, 7, 9, 1 << 40, 12345] {
            assert_eq!(t.counts.get(h), c.get(h), "hash {h}");
        }
        assert_eq!(t.counts.len(), c.len());
        std::fs::remove_dir_all(&dir).unwrap();
    }
}

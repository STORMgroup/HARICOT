//! Raw-signal input: pod5 files (pure-Rust reader via the `pod5lib` crate).

use std::path::Path;

use anyhow::{Result, anyhow, bail};

/// A single read's raw current signal, already calibrated to picoamps (pA).
#[derive(Debug, Clone)]
pub struct RawRead {
    /// Read UUID string.
    pub id: String,
    /// Calibrated current samples in pA, one per acquisition sample.
    pub signal_pa: Vec<f32>,
}

/// Hand each read to `f` as it is decoded, never materialising the whole file.
///
/// Returning a `Vec<RawRead>` would make the entire signal live
/// at once — the single largest term in peak RSS (32.65 GB of D. melanogaster fwd's
/// 49.40 GB). Callers that only need to *transform* each read (sketching, counting) can
/// consume it a read at a time instead and keep the signal's footprint constant.
///
/// The reads in file order, so `read_id` assignment by
/// arrival index is unchanged.
pub fn for_each_raw_read<F>(path: &Path, limit: Option<usize>, mut f: F) -> Result<()>
where
    F: FnMut(RawRead) -> Result<()>,
{
    match path.extension().and_then(|e| e.to_str()) {
        Some("pod5") => {
            use pod5lib::Reader;
            let reader =
                Reader::open(path).map_err(|e| anyhow!("opening pod5 {}: {e}", path.display()))?;
            // Decompressed in parallel batches, yielded in file order: the same samples as
            // `signal_pa()` per read, whose file lock also covers the decompression and so
            // ran it on one thread.
            let threads = rayon::current_num_threads();
            for (n, r) in reader.par_reads_iter(threads).enumerate() {
                if limit.is_some_and(|l| n >= l) {
                    break;
                }
                let (read, adc) =
                    r.map_err(|e| anyhow!("decoding signal in {}: {e}", path.display()))?;
                f(RawRead {
                    id: read.read_id_str().to_string(),
                    signal_pa: read.calibration.to_pa(&adc),
                })?;
            }
            Ok(())
        }
        other => bail!(
            "unrecognized signal file extension {:?} (expected .pod5)",
            other
        ),
    }
}

/// [`for_each_raw_read`] in chunks of about `chunk_bytes / 2` of signal, decoded on a second
/// thread while `f` processes the previous chunk: the same reads in the same order, with the
/// decoding (whose file reads are serial) hidden behind the processing instead of alternating
/// with it. At most two chunks are alive, so the signal's footprint stays about `chunk_bytes`.
pub fn for_each_raw_chunk<F>(
    path: &Path,
    limit: Option<usize>,
    chunk_bytes: usize,
    mut f: F,
) -> Result<()>
where
    F: FnMut(Vec<RawRead>) -> Result<()>,
{
    let half = (chunk_bytes / 2).max(1);
    let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<RawRead>>(0);
    std::thread::scope(|scope| {
        let producer = scope.spawn(move || -> Result<()> {
            let (mut chunk, mut bytes) = (Vec::new(), 0usize);
            for_each_raw_read(path, limit, |r| {
                bytes += r.signal_pa.len() * std::mem::size_of::<f32>();
                chunk.push(r);
                if bytes >= half {
                    bytes = 0;
                    tx.send(std::mem::take(&mut chunk))
                        .map_err(|_| anyhow!("signal consumer stopped"))?;
                }
                Ok(())
            })?;
            if !chunk.is_empty() {
                tx.send(chunk)
                    .map_err(|_| anyhow!("signal consumer stopped"))?;
            }
            Ok(())
        });
        let consumed = rx.iter().try_for_each(&mut f);
        drop(rx);
        let produced = producer.join().expect("the signal decoder panicked");
        // A consumer error stops the producer too; report the consumer's.
        consumed.and(produced)
    })
}

/// A read's raw signal as stored: ADC counts and the calibration that turns them into picoamps.
/// Half the bytes of [`RawRead`]'s `f32` samples, so a chunk of these holds the same reads in half
/// the memory; [`AdcRead::signal_pa`] is the `f32` signal [`RawRead`] holds for the same read (the
/// same `pod5lib` conversion, applied later).
#[derive(Debug, Clone)]
pub struct AdcRead {
    /// Read UUID string.
    pub id: String,
    adc: Vec<i16>,
    calibration: pod5lib::Calibration,
}

impl AdcRead {
    /// The calibrated signal in pA: [`RawRead`]'s `signal_pa` of the same read, sample for sample.
    pub fn signal_pa(&self) -> Vec<f32> {
        self.calibration.to_pa(&self.adc)
    }
}

/// [`for_each_raw_chunk`] for pod5 files with [`AdcRead`] chunks: the same reads in the same order and
/// the same chunk boundaries (a chunk closes when its samples, counted as `f32`, reach
/// `chunk_bytes / 2`, as there), each chunk in half the memory; the consumer makes each read's `f32`
/// signal when it processes it.
pub fn for_each_adc_chunk<F>(
    path: &Path,
    limit: Option<usize>,
    chunk_bytes: usize,
    mut f: F,
) -> Result<()>
where
    F: FnMut(Vec<AdcRead>) -> Result<()>,
{
    use pod5lib::Reader;
    if path.extension().and_then(|e| e.to_str()) != Some("pod5") {
        bail!(
            "{}: ADC chunks are read from pod5 files only",
            path.display()
        );
    }
    let half = (chunk_bytes / 2).max(1);
    let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<AdcRead>>(0);
    std::thread::scope(|scope| {
        let producer = scope.spawn(move || -> Result<()> {
            let reader =
                Reader::open(path).map_err(|e| anyhow!("opening pod5 {}: {e}", path.display()))?;
            let threads = rayon::current_num_threads();
            let (mut chunk, mut bytes) = (Vec::new(), 0usize);
            for (n, r) in reader.par_reads_iter(threads).enumerate() {
                if limit.is_some_and(|l| n >= l) {
                    break;
                }
                let (read, adc) =
                    r.map_err(|e| anyhow!("decoding signal in {}: {e}", path.display()))?;
                bytes += adc.len() * std::mem::size_of::<f32>();
                chunk.push(AdcRead {
                    id: read.read_id_str().to_string(),
                    adc,
                    calibration: read.calibration.clone(),
                });
                if bytes >= half {
                    bytes = 0;
                    tx.send(std::mem::take(&mut chunk))
                        .map_err(|_| anyhow!("signal consumer stopped"))?;
                }
            }
            if !chunk.is_empty() {
                tx.send(chunk)
                    .map_err(|_| anyhow!("signal consumer stopped"))?;
            }
            Ok(())
        });
        let consumed = rx.iter().try_for_each(&mut f);
        drop(rx);
        let produced = producer.join().expect("the signal decoder panicked");
        // A consumer error stops the producer too; report the consumer's.
        consumed.and(produced)
    })
}

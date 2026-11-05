//! The candidate cache, shared by the pipeline programs (`#[path]`-included).
//!
//! Format: bincode of `(counts, names, overlaps)` (every read's length in retained events, its read
//! id, the overlaps), then bincode of every read's signal content hash (`Vec<u64>`,
//! `haricot::readorder::content_hash`), the key of the tie rank. bincode 1's `deserialize` allows
//! trailing bytes, so a reader of the triple alone reads either form; a cache without the hashes
//! decodes with `content: None`, and the steps that break ties (the canonical pair frame,
//! `--canonical-ties`, the rescue's admit) refuse it rather than fall back to the read ids.
#![allow(dead_code)]

use anyhow::{Context, Result, ensure};
use haricot::overlap::Overlap;

/// A decoded cache.
pub struct Cache {
    pub counts: Vec<u32>,
    pub names: Vec<String>,
    pub ovl: Vec<Overlap>,
    /// Every read's signal content hash, by read number; `None` in a cache without them.
    pub content: Option<Vec<u64>>,
}

impl Cache {
    /// Every read's tie rank ([`tie_rank`]).
    pub fn tie_rank(&self) -> Result<Vec<u32>> {
        tie_rank(self.content.as_deref())
    }
}

/// Every read's tie rank from a cache's content hashes (`readorder::rank_by_content`). Errors
/// without them rather than fall back to the read ids.
pub fn tie_rank(content: Option<&[u64]>) -> Result<Vec<u32>> {
    let c = content.context(
        "this candidate cache has no read content hashes: the tie rank needs them; re-run the \
         detection with this build",
    )?;
    haricot::readorder::rank_by_content(c)
}

/// Decode a cache from its bytes.
pub fn decode(bytes: &[u8]) -> Result<Cache> {
    let (counts, names, ovl): (Vec<u32>, Vec<String>, Vec<Overlap>) =
        bincode::deserialize(bytes).context("decoding the candidate cache")?;
    let n = bincode::serialized_size(&(&counts, &names, &ovl))? as usize;
    let content = if bytes.len() > n {
        let c: Vec<u64> =
            bincode::deserialize(&bytes[n..]).context("decoding the cache's content hashes")?;
        ensure!(
            c.len() == names.len(),
            "the cache holds {} content hashes for {} reads",
            c.len(),
            names.len()
        );
        Some(c)
    } else {
        None
    };
    Ok(Cache {
        counts,
        names,
        ovl,
        content,
    })
}

/// Read and decode the cache at `path`.
pub fn read(path: &str) -> Result<Cache> {
    let bytes = std::fs::read(path).with_context(|| format!("reading {path}"))?;
    decode(&bytes).with_context(|| format!("decoding {path}"))
}

/// Encode a cache: the triple, then the content hashes when there are some.
pub fn encode(
    counts: &[u32],
    names: &[String],
    ovl: &[Overlap],
    content: Option<&[u64]>,
) -> Result<Vec<u8>> {
    let mut b = bincode::serialize(&(counts, names, ovl)).context("encoding the cache")?;
    if let Some(c) = content {
        ensure!(
            c.len() == names.len(),
            "{} content hashes for {} reads",
            c.len(),
            names.len()
        );
        b.extend(bincode::serialize(c).context("encoding the content hashes")?);
    }
    Ok(b)
}

/// Encode and write a cache to `path`.
pub fn write(
    path: &str,
    counts: &[u32],
    names: &[String],
    ovl: &[Overlap],
    content: Option<&[u64]>,
) -> Result<()> {
    let b = encode(counts, names, ovl, content)?;
    std::fs::write(path, b).with_context(|| format!("writing {path}"))
}

/// `read_content.tsv`: every read's content hash (16 hex digits) and id, in hash order, so the file
/// is the same whatever the reads' order in the pod5 (and, ids mapped back, whatever their names).
/// Read by the end bridges (`--content-ties`) for their score ties.
pub fn content_table(names: &[String], content: &[u64]) -> String {
    let mut rows: Vec<(u64, &str)> = content
        .iter()
        .zip(names)
        .map(|(&h, n)| (h, n.as_str()))
        .collect();
    rows.sort_unstable();
    let mut s = String::from("content\tname\n");
    for (h, n) in rows {
        s.push_str(&format!("{h:016x}\t{n}\n"));
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use haricot::params::Strand;

    fn ovl(a: u32, b: u32) -> Overlap {
        Overlap {
            a,
            b,
            strand: Strand::Forward,
            n_anchors: 9,
            matched: 400,
            cnt: 9,
            aligned_bins: 5,
            bin_events: 192,
            gap: 1,
            slope: 1.01,
            a_bin_span: (0, 4),
            b_bin_span: (2, 6),
            path: vec![(0, 2), (1, 3), (3, 5), (4, 6)],
            z: 6.5,
        }
    }

    /// A cache with hashes round-trips, a cache without them reads as `None`, and a
    /// reader of the triple alone reads either.
    #[test]
    fn caches_round_trip_with_and_without_content() {
        let counts: Vec<u32> = [900, 1200, 700].to_vec();
        let names: Vec<String> = ["r0", "r1", "r2"].map(String::from).to_vec();
        let o: Vec<Overlap> = [ovl(0, 1), ovl(1, 2)].to_vec();
        let with = encode(&counts, &names, &o, Some(&[7u64, 3, 5][..])).unwrap();
        let c = decode(&with).unwrap();
        assert_eq!(c.counts, counts);
        assert_eq!(c.names, names);
        assert_eq!(c.ovl.len(), 2);
        assert_eq!(c.content.as_deref(), Some(&[7u64, 3, 5][..]));
        assert_eq!(c.tie_rank().unwrap(), [2, 0, 1]);
        let old = bincode::serialize(&(&counts, &names, &o)).unwrap();
        let c = decode(&old).unwrap();
        assert!(c.content.is_none() && c.tie_rank().is_err());
        let t: (Vec<u32>, Vec<String>, Vec<Overlap>) = bincode::deserialize(&with).unwrap();
        assert_eq!(t.1, names);
        assert!(encode(&counts, &names, &o, Some(&[1u64, 2][..])).is_err());
        let dup = encode(&counts, &names, &o, Some(&[1u64, 2, 1][..])).unwrap();
        assert!(decode(&dup).unwrap().tie_rank().is_err());
        assert_eq!(
            content_table(&names, &[7, 3, 5]),
            "content\tname\n0000000000000003\tr1\n0000000000000005\tr2\n0000000000000007\tr0\n"
        );
    }
}

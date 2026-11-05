# datasets/: the manuscript's read sets from public data

The manuscript used 30 training sets and 56 final-test sets. Each set is a list of read ids. These scripts download
the public raw signal (POD5, or fast5 converted to POD5), write each set's POD5 in its saved read order and check
every read. A second script fetches the reference genomes and checks every sequence.

## Read lists

`haricot-read-lists-v1.tar.gz` comes from the Zenodo record [10.5281/zenodo.23963126](https://doi.org/10.5281/zenodo.23963126).
Unpack it in the checkout; it creates `read_lists/`:

```bash
curl -L -o haricot-read-lists-v1.tar.gz \
  'https://zenodo.org/records/23963126/files/haricot-read-lists-v1.tar.gz?download=1'
echo '4c64266e21382a1a47d4e0c75a47b1d27c47b26ec3100a51fc9dd0b8f3b9bd92  haricot-read-lists-v1.tar.gz' | sha256sum -c
tar -xzf haricot-read-lists-v1.tar.gz     # -> read_lists/
```

- `sets.tsv`: one row per set (role, species, reads, public source, parent set of a fragment set, SHA-256 of
  `read_ids.txt`).
- `SHA256SUMS`: every file (`cd read_lists && sha256sum -c SHA256SUMS`).
- `<role>/<set>/read_ids.txt`: the read ids in the order of the set's POD5.
- `<role>/<set>/checks.tsv`: per read, the samples, the MD5 of the raw int16 signal, the calibration offset and
  scale, and the sample rate.
- `<role>/<set>/index.tsv`: per read, its public source and file (as in `files.tsv`), run id and samples.
- `<role>/<set>/fragments.tsv`: 4 kb fragment sets only (cut from the parent set by `truth/fragments.sh`).

## Files

| file | contents |
|---|---|
| `sources.tsv` | one row per public source: type, URL, published MD5, archive, accession, sets, licence as the source states it (and where), citation |
| `files.tsv` | one row per public file: source, path, bytes, offset in the tar (figshare) or S3 ETag |
| `download.py` | downloads the files of some sets (`--set`) or of a source (`--source`), checks sizes and checksums |
| `build_set.py` | writes a set's POD5 from the downloaded files in `read_ids.txt` order; checks every read against `checks.tsv` |
| `subsets.tsv` | the depth subsets (`_c12`, `_c6`): the parent set and the base depth of `truth/subset.py`'s rule |
| `references.tsv` | one row per reference file: public source, URL, licence as stated (and where) |
| `reference_sequences.tsv` | one row per reference sequence: file, name, length, MD5, name in the source |
| `references.py` | downloads the references, checks every sequence's length and MD5, writes the files |

## Usage

Python 3 with the `pod5` package (`envs/haricot.yml`) and `curl`.

```bash
python3 datasets/download.py downloads --set read_lists/final_test/sd_fwd_c6 --list   # files and size, no download
python3 datasets/download.py downloads --set read_lists/final_test/sd_fwd_c6
python3 datasets/build_set.py read_lists/final_test/sd_fwd_c6 downloads -o sd_fwd_c6.pod5
python3 datasets/references.py refs --files sd.fa                             # refs/sd.fa
python3 datasets/references.py refs                                           # all 14 (HG002: 1.8 GB, ~4 GB RAM)
```

`build_set.py` exits with status 1 if a read is missing or differs from `checks.tsv` in samples, signal MD5,
calibration or sample rate. `download.py --source SOURCE` fetches a whole source; with a figshare or ENA archive
this is the only way to check the archive's published MD5 (`--set` fetches figshare members by byte range).

- Fragment sets (`*_L4000*`): cut from the parent set by `truth/fragments.sh`; build the parent set here.
- Depth subsets (`_c12`, `_c6`): `build_set.py` builds them from their own read list. `truth/subset.py` with the base
  depth of `subsets.tsv` (the parent's depth: 30, or what the source held) selects the same reads from the parent
  set, in a different order.

## Sources

| source | data | sets | files | size |
|---|---|---|---:|---:|
| `hall2024_cj`, `_lm`, `_vp`, `_se`, `_sd` | Hall et al. 2024 POD5 tars, University of Melbourne figshare | cj, lm, vp1, vp2, se, sd | one tar each | 42, 60, 85, 198, 11 GB |
| `ena_athal` | ENA ERR17795792, `ATRun3_pod5.tar.gz` | at2, at3, at5 | 1 tar.gz (75 POD5) | 164 GB |
| `ont_dmel` | ONT Open Data `melanogaster_bkim_2023.01`, five runs, fast5 (4 kHz) | dm2L, dm2R, dm3L, dmX | 504 | 228 GB |
| `ont_hg002` | ONT Open Data `giab_2023.05`, HG002 flow cell PAO89685, `pod5_pass` | hg18 | 2,182 | 1,446 GB |
| `ont_nomiss` | ONT Open Data `nomiss_96BC_P2I_SUP_2026`, raw POD5 of all 96 barcodes | sa | 293 | 1,728 GB |
| `basemod` | ONT methylation benchmark (Kulkarni et al.), `ont-basemod-benchmark-data` | ek, hp2, td2, an2 | 4 | 19 GB |

- A set needs most files of its source; `download.py --set ... --list` prints how many. Examples: `sd_fwd_c6` 58 of
  60 files (10.3 GB), `se_fwd_c6` 1,453 of 3,569 (101 GB), `hg18_fwd` all 2,182, `sa_fwd` 282 of 293 (1.65 TB). The
  A. thaliana archive is gzip: `--set` streams all of it (164 GB).
- The ONT Open Data bucket is CC BY-NC 4.0 (non-commercial). Licences of the other sources: `sources.tsv`.
- D. melanogaster: the runs are fast5; `download.py` converts each file with `pod5 convert fast5` and removes the
  fast5. `build_set.py` reads the converted `.pod5` where `index.tsv` names the `.fast5`.
- S. aureus: the sets hold isolate 8 (rapid barcodes 08, 32, 56, 80); the raw POD5 files mix all 96 barcodes.
- ENA: the archive URL is http (ENA serves the same file over http, https and ftp). The archive MD5 (whole download)
  and `build_set.py`'s per-read MD5s check the content.
- `references.py` holds the HG002 sequences in memory (about 4 GB).

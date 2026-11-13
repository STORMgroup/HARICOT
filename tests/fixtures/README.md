# Test fixtures

`debug_reads.pod5` (not distributed): two overlapping R10.4.1 reads (5 kHz) of *E. coli* ATCC 25922, for
`tests/e2e.rs`.

- Read ids: `3c3e1cc6-25e4-4abd-b15e-2eea8c7bcb39`, `97c3ef03-5a6a-469e-a7d0-3f74136c571f`.
- Source: Hall et al. 2024, "ATCC_25922__202309 pod5 data", University of Melbourne figshare,
  [10.26188/25521892](https://doi.org/10.26188/25521892), CC BY 4.0; the archive is
  `https://ndownloader.figshare.com/files/45408628` (`ATCC_25922__202309.tar`).
- Without the file, `cargo test` reports the three e2e tests as ignored.

To make it (the whole archive is downloaded; `pod5` comes with `envs/haricot.yml`):

```bash
curl -L -o ATCC_25922__202309.tar https://ndownloader.figshare.com/files/45408628
tar -xf ATCC_25922__202309.tar
printf '%s\n' 3c3e1cc6-25e4-4abd-b15e-2eea8c7bcb39 97c3ef03-5a6a-469e-a7d0-3f74136c571f > debug_ids.txt
pod5 filter $(tar -tf ATCC_25922__202309.tar | grep '\.pod5$') --ids debug_ids.txt \
    --output tests/fixtures/debug_reads.pod5
```

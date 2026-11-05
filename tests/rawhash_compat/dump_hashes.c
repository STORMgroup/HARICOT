/* dump_hashes.c — golden-fixture generator for the RawHash2 <-> Rust hash
 * cross-check (see tests/hash_compat.rs and README.md).
 *
 * Two modes:
 *
 *   dump_hashes gen  <n> <seed>                       > signal.txt
 *       Emit <n> deterministic pseudo-random normalized event values (one
 *       per line, %.9g so they round-trip to the same f32 in C and Rust).
 *       A tiny self-contained LCG is used so the fixture never depends on
 *       libc rand() or any RNG that could drift between platforms.
 *
 *   dump_hashes dump <signal.txt> <e> <w> <quant_bit> <n_buckets> <k> \
 *                    <diff> <fine_min> <fine_max> <fine_range>        > golden.tsv
 *       Read the signal, run RawHash2's ri_sketch(), and print one
 *       "pos<TAB>hash" line per emitted seed. `hash` is the bare 64-bit
 *       hash value with RawHash2's packing stripped (x >> RI_HASH_SHIFT);
 *       `pos` is the first-event index of the e-mer (low 32 bits of y,
 *       minus the strand bit). With w=0 the seeds come out in event order,
 *       which is exactly the order src/emer/sketch.rs produces them.
 *
 * Build: see regenerate.sh (links extern/RawHash2/src/{rsketch,kalloc}.c).
 * This file is NOT compiled by `cargo test`; it is a one-shot generator whose
 * committed output (signal.txt + *.tsv) is what the Rust test consumes.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rsketch.h" /* ri_sketch, RI_HASH_SHIFT, RI_POS_SHIFT */
#include "rutils.h"  /* mm128_v */

static int gen_signal(int n, unsigned long long seed) {
    /* Numerical Recipes 64-bit LCG; deterministic and platform-independent. */
    unsigned long long s = seed;
    for (int i = 0; i < n; ++i) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        /* top 24 bits -> uniform in [0,1), then map to [-3.5, 3.5] so the
         * value spans the low / fine / high quantization regions (and the
         * clamped extremes beyond +-3). */
        double u = (double)(s >> 40) / (double)(1ULL << 24);
        float v = (float)(u * 7.0 - 3.5);
        printf("%.9g\n", v);
    }
    return 0;
}

static int dump(const char *path, int e, int w, unsigned int quant_bit,
                unsigned int n_buckets, int k, int diff, float fine_min,
                float fine_max, float fine_range) {
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "cannot open %s\n", path);
        return 1;
    }
    size_t cap = 1024, len = 0;
    float *sig = (float *)malloc(cap * sizeof(float));
    while (fscanf(fp, "%f", &sig[len]) == 1) {
        if (++len == cap) {
            cap *= 2;
            sig = (float *)realloc(sig, cap * sizeof(float));
        }
    }
    fclose(fp);

    mm128_v p = {0, 0, 0};
    /* km = NULL -> the kvec macros fall back to plain realloc(). */
    ri_sketch(NULL, sig, /*id=*/0, /*strand=*/0, (uint32_t)len, diff, w, e,
              /*n=*/0, quant_bit, n_buckets, k, fine_min, fine_max, fine_range,
              &p, /*out=*/0);

    for (size_t i = 0; i < p.n; ++i) {
        uint32_t pos = (uint32_t)(p.a[i].y & 0xffffffffULL) >> RI_POS_SHIFT;
        uint64_t hash = p.a[i].x >> RI_HASH_SHIFT;
        printf("%u\t%llu\n", pos, (unsigned long long)hash);
    }

    free(p.a);
    free(sig);
    return 0;
}

int main(int argc, char **argv) {
    if (argc >= 4 && strcmp(argv[1], "gen") == 0) {
        return gen_signal(atoi(argv[2]), strtoull(argv[3], NULL, 10));
    }
    if (argc == 12 && strcmp(argv[1], "dump") == 0) {
        return dump(argv[2], atoi(argv[3]), atoi(argv[4]),
                    (unsigned)atoi(argv[5]), (unsigned)atoi(argv[6]),
                    atoi(argv[7]), atoi(argv[8]), (float)atof(argv[9]),
                    (float)atof(argv[10]), (float)atof(argv[11]));
    }
    fprintf(stderr,
            "usage:\n"
            "  %s gen  <n> <seed>\n"
            "  %s dump <signal.txt> <e> <w> <quant_bit> <n_buckets> <k> "
            "<diff> <fine_min> <fine_max> <fine_range>\n",
            argv[0], argv[0]);
    return 2;
}

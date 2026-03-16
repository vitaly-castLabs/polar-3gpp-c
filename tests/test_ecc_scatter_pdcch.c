#include "polar_3gpp.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "ecc_common.h"

int main(int argc, char **argv) {
    int num_iterations = 1000;
    int num_bits;

    if (argc >= 2) {
        num_iterations = atoi(argv[1]);
        if (num_iterations < 1) {
            num_iterations = 1;
        }
    }
    if (argc >= 3) {
        ecc_seed = (uint32_t)strtoul(argv[2], NULL, 10);
    }

    for (num_bits = 15; num_bits <= 30; num_bits++) {
        int succ = 0;
        int succ_repeat = 0;
        int wrong_dec = 0;
        int i;

        for (i = 0; i < num_iterations; i++) {
            const int A = 13;
            const int E = 8040;
            uint8_t a[A];
            double *f_tilde = (double *)malloc((size_t)E * sizeof(double));
            int *indices = (int *)malloc((size_t)E * sizeof(int));
            uint8_t seen_mod[A];
            polar_u8_vec_t enc = {0};
            polar_u8_vec_t dec = {0};
            polar_status_t st;
            int j;
            int all_mods_seen = 1;

            if (f_tilde == NULL || indices == NULL) {
                fprintf(stderr, "Allocation failed\n");
                free(f_tilde);
                free(indices);
                return 1;
            }

            ecc_fill_random_bits(a, A);

            st = polar_pdcch_encode(a, A, E, NULL, &enc);
            if (st != POLAR_STATUS_OK) {
                fprintf(stderr, "PDCCH encode failed: status=%d\n", st);
                free(f_tilde);
                free(indices);
                return 1;
            }

            ecc_bits_to_llr(enc.data, enc.len, f_tilde);

            for (j = 0; j < E; j++) {
                indices[j] = j;
            }
            ecc_shuffle_i32(indices, E);

            for (j = 0; j < A; j++) {
                seen_mod[j] = 0;
            }
            for (j = 0; j < num_bits; j++) {
                int residue = indices[j] % A;
                seen_mod[residue] = 1;
            }
            for (j = 0; j < A; j++) {
                if (!seen_mod[j]) {
                    all_mods_seen = 0;
                    break;
                }
            }
            if (all_mods_seen) {
                succ_repeat++;
            }

            for (j = num_bits; j < E; j++) {
                f_tilde[indices[j]] = 0.0;
            }

            st = polar_pdcch_decode(f_tilde, E, A, 8, 1, NULL, &dec);
            if (st != POLAR_STATUS_OK) {
                fprintf(stderr, "PDCCH decode failed: status=%d\n", st);
                polar_free(enc.data);
                free(f_tilde);
                free(indices);
                return 1;
            }

            if (dec.len == A && ecc_equal_bits(a, dec.data, A)) {
                succ++;
            } else if (dec.len > 0) {
                wrong_dec++;
            }

            polar_free(enc.data);
            polar_free(dec.data);
            free(f_tilde);
            free(indices);
        }

        printf("Scatter pdcch: %d bits, decoding success: %.1f%%, wrong corrections: %.1f%%, repeat code success: %.1f%%\n",
               num_bits,
               (100.0 * (double)succ) / (double)num_iterations,
               (100.0 * (double)wrong_dec) / (double)num_iterations,
               (100.0 * (double)succ_repeat) / (double)num_iterations);
    }

    return 0;
}

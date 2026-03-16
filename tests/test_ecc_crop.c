#include "polar_3gpp.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "ecc_common.h"

int main(int argc, char **argv) {
    int num_iterations = 1000;
    int crop;

    if (argc >= 2) {
        num_iterations = atoi(argv[1]);
        if (num_iterations < 1) {
            num_iterations = 1;
        }
    }
    if (argc >= 3) {
        ecc_seed = (uint32_t)strtoul(argv[2], NULL, 10);
    }

    for (crop = 15; crop <= 30; crop++) {
        int succ = 0;
        int wrong_dec = 0;
        int i;

        for (i = 0; i < num_iterations; i++) {
            const int A = 13;
            const int E = 8040;
            uint8_t a[A];
            double *f_tilde = (double *)malloc((size_t)E * sizeof(double));
            polar_u8_vec_t enc = {0};
            polar_u8_vec_t dec = {0};
            polar_status_t st;
            int range_start;
            int preserve_start;
            int preserve_end;
            int j;

            if (f_tilde == NULL) {
                fprintf(stderr, "Allocation failed\n");
                return 1;
            }

            ecc_fill_random_bits(a, A);

            st = polar_pucch_encode(a, A, E, &enc);
            if (st != POLAR_STATUS_OK) {
                fprintf(stderr, "Encode failed: status=%d\n", st);
                free(f_tilde);
                return 1;
            }

            ecc_bits_to_llr(enc.data, enc.len, f_tilde);

            range_start = ecc_rand_range(0, E - 1);
            preserve_start = range_start;
            preserve_end = preserve_start + crop;
            if (preserve_end > E) {
                preserve_end = E;
                preserve_start = preserve_end - crop;
            }

            for (j = 0; j < preserve_start; j++) {
                f_tilde[j] = 0.0;
            }
            for (j = preserve_end; j < E; j++) {
                f_tilde[j] = 0.0;
            }

            st = polar_pucch_decode(f_tilde, E, A, 8, 1, &dec);
            if (st != POLAR_STATUS_OK) {
                fprintf(stderr, "Decode failed: status=%d\n", st);
                polar_free(enc.data);
                free(f_tilde);
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
        }

        printf("Crop: %d bits, decoding success: %.1f%%, wrong corrections: %.1f%%\n",
               crop,
               (100.0 * (double)succ) / (double)num_iterations,
               (100.0 * (double)wrong_dec) / (double)num_iterations);
    }

    return 0;
}

#include "polar_3gpp.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "ecc_common.h"

int main(int argc, char **argv) {
    int num_iterations = 3;
    int iter;
    int failures = 0;

    if (argc >= 2) {
        num_iterations = atoi(argv[1]);
        if (num_iterations < 1) {
            num_iterations = 1;
        }
    }
    if (argc >= 3) {
        ecc_seed = (uint32_t)strtoul(argv[2], NULL, 10);
    }

    for (iter = 0; iter < num_iterations; iter++) {
        int A = ecc_rand_range(12, 32);
        int E = A * 3;
        int num_zeros;
        int i;
        uint8_t *a = (uint8_t *)malloc((size_t)A);
        double *f_tilde = (double *)malloc((size_t)E * sizeof(double));
        int *indices = (int *)malloc((size_t)E * sizeof(int));
        polar_u8_vec_t enc = {0};
        polar_u8_vec_t dec = {0};
        polar_status_t st;

        if (a == NULL || f_tilde == NULL || indices == NULL) {
            fprintf(stderr, "Allocation failed\n");
            free(a);
            free(f_tilde);
            free(indices);
            return 1;
        }

        ecc_fill_random_bits(a, A);

        st = polar_pucch_encode(a, A, E, &enc);
        if (st != POLAR_STATUS_OK) {
            fprintf(stderr, "Encode failed at iter %d: status=%d\n", iter + 1, st);
            free(a);
            free(f_tilde);
            free(indices);
            return 1;
        }

        ecc_bits_to_llr(enc.data, enc.len, f_tilde);

        num_zeros = E - (A + 6);
        if (num_zeros < 0) {
            num_zeros = 0;
        }

        for (i = 0; i < E; i++) {
            indices[i] = i;
        }
        ecc_shuffle_i32(indices, E);
        for (i = 0; i < num_zeros; i++) {
            f_tilde[indices[i]] = 0.0;
        }

        st = polar_pucch_decode(f_tilde, E, A, 8, 1, &dec);
        if (st != POLAR_STATUS_OK) {
            fprintf(stderr, "Decode failed at iter %d: status=%d\n", iter + 1, st);
            polar_free(enc.data);
            free(a);
            free(f_tilde);
            free(indices);
            return 1;
        }

        printf("Iter %d: num erasures %d/%d, intact %d -> %s\n",
               iter + 1,
               num_zeros,
               E,
               E - num_zeros,
               (dec.len == A && ecc_equal_bits(a, dec.data, A)) ? "success" : "failure");

        if (!(dec.len == A && ecc_equal_bits(a, dec.data, A))) {
            failures++;
        }

        polar_free(enc.data);
        polar_free(dec.data);
        free(a);
        free(f_tilde);
        free(indices);
    }

    printf("Summary: %d/%d failures\n", failures, num_iterations);
    return 0;
}

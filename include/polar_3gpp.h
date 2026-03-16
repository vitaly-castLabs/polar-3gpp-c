#ifndef POLAR_3GPP_H
#define POLAR_3GPP_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum polar_status {
    POLAR_STATUS_OK = 0,
    POLAR_STATUS_UNSUPPORTED_BLOCK_LENGTH = 1,
    POLAR_STATUS_INVALID_ARGUMENT = 2,
    POLAR_STATUS_DECODING_FAILED = 3,
    POLAR_STATUS_ALLOCATION_FAILED = 4
} polar_status_t;

typedef struct polar_u8_vec {
    uint8_t *data;
    int len;
} polar_u8_vec_t;

typedef struct polar_i32_vec {
    int *data;
    int len;
} polar_i32_vec_t;

void polar_free(void *ptr);

polar_status_t polar_pbch_encode(const uint8_t *a, int A, int E, polar_u8_vec_t *f);
polar_status_t polar_pbch_decode(
    const double *f_tilde,
    int E,
    int A,
    int L,
    int min_sum,
    const int8_t *a_tilde,
    polar_u8_vec_t *a_hat
);

polar_status_t polar_pdcch_encode(
    const uint8_t *a,
    int A,
    int E,
    const uint8_t *rnti,
    polar_u8_vec_t *f
);
polar_status_t polar_pdcch_decode(
    const double *f_tilde,
    int E,
    int A,
    int L,
    int min_sum,
    const uint8_t *rnti,
    polar_u8_vec_t *a_hat
);

polar_status_t polar_pucch_encode(const uint8_t *a, int A, int G, polar_u8_vec_t *f);
polar_status_t polar_pucch_decode(
    const double *f_tilde,
    int G,
    int A,
    int L,
    int min_sum,
    polar_u8_vec_t *a_hat
);

polar_status_t polar_custom1_encode(const uint8_t *a, int A, int E, polar_u8_vec_t *f);
polar_status_t polar_custom1_decode(
    const double *f_tilde,
    int E,
    int A,
    int L,
    int min_sum,
    polar_u8_vec_t *a_hat
);

polar_status_t polar_get_3gpp_n(int K, int E, int n_max, int *N_out);
polar_status_t polar_get_3gpp_rate_matching_pattern(
    int K,
    int N,
    int E,
    polar_i32_vec_t *pattern,
    int *mode_out
);
polar_status_t polar_get_3gpp_channel_interleaver_pattern(int E, polar_i32_vec_t *pattern);

#ifdef __cplusplus
}
#endif

#endif

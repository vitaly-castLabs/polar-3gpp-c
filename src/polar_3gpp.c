#include "polar_3gpp.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {
    RM_MODE_REPETITION = 0,
    RM_MODE_PUNCTURING = 1,
    RM_MODE_SHORTENING = 2
};

typedef struct pm_idx {
    double pm;
    int idx;
} pm_idx_t;

typedef struct u8_matrix {
    uint8_t *data;
    int rows;
    int cols;
} u8_matrix_t;

typedef struct scl_ctx {
    int N;
    int n;
    int n_cols;
    int L_target;
    int L_capacity;
    int L_prime;
    int min_sum;
    int path_stride;

    uint8_t *bits;
    uint8_t *bits_alt;
    uint8_t *bits_updated;

    double *llrs;
    double *llrs_alt;
    uint8_t *llrs_updated;

    double *PM;
    double *PM_alt;
} scl_ctx_t;

#define U8_AT(mat, r, c) ((mat).data[(size_t)(r) * (size_t)(mat).cols + (size_t)(c)])
#define SCL_BITS_IDX(ctx, l, r, c) ((size_t)(l) * (size_t)(ctx)->path_stride + (size_t)(r) * (size_t)(ctx)->n_cols + (size_t)(c))
#define SCL_LLR_IDX(ctx, l, r, c) ((size_t)(l) * (size_t)(ctx)->path_stride + (size_t)(r) * (size_t)(ctx)->n_cols + (size_t)(c))

static const int Q_NMAX[] = {
#include "q_nmax.inc"
};

static const int PI_IL_MAX[] = {
#include "pi_il_max.inc"
};

static const int P32[32] = {
    0, 1, 2, 4, 3, 5, 6, 7,
    8, 16, 9, 17, 10, 18, 11, 19,
    12, 20, 13, 21, 14, 22, 15, 23,
    24, 25, 26, 28, 27, 29, 30, 31
};

static const uint8_t CRC24_POLY[] = {
    1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1
};

static const uint8_t CRC11_POLY[] = {
    1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1
};

static const uint8_t CRC6_POLY[] = {
    1, 1, 0, 0, 0, 0, 1
};

static polar_status_t alloc_u8(size_t n, uint8_t **out) {
    if (n == 0) {
        *out = NULL;
        return POLAR_STATUS_OK;
    }
    *out = (uint8_t *)malloc(n);
    if (*out == NULL) {
        return POLAR_STATUS_ALLOCATION_FAILED;
    }
    return POLAR_STATUS_OK;
}

static polar_status_t alloc_i32(size_t n, int **out) {
    if (n == 0) {
        *out = NULL;
        return POLAR_STATUS_OK;
    }
    *out = (int *)malloc(n * sizeof(int));
    if (*out == NULL) {
        return POLAR_STATUS_ALLOCATION_FAILED;
    }
    return POLAR_STATUS_OK;
}

static polar_status_t alloc_f64(size_t n, double **out) {
    if (n == 0) {
        *out = NULL;
        return POLAR_STATUS_OK;
    }
    *out = (double *)malloc(n * sizeof(double));
    if (*out == NULL) {
        return POLAR_STATUS_ALLOCATION_FAILED;
    }
    return POLAR_STATUS_OK;
}

static polar_status_t calloc_u8(size_t n, uint8_t **out) {
    if (n == 0) {
        *out = NULL;
        return POLAR_STATUS_OK;
    }
    *out = (uint8_t *)calloc(n, 1);
    if (*out == NULL) {
        return POLAR_STATUS_ALLOCATION_FAILED;
    }
    return POLAR_STATUS_OK;
}

static polar_status_t calloc_i32(size_t n, int **out) {
    if (n == 0) {
        *out = NULL;
        return POLAR_STATUS_OK;
    }
    *out = (int *)calloc(n, sizeof(int));
    if (*out == NULL) {
        return POLAR_STATUS_ALLOCATION_FAILED;
    }
    return POLAR_STATUS_OK;
}

static polar_status_t calloc_f64(size_t n, double **out) {
    if (n == 0) {
        *out = NULL;
        return POLAR_STATUS_OK;
    }
    *out = (double *)calloc(n, sizeof(double));
    if (*out == NULL) {
        return POLAR_STATUS_ALLOCATION_FAILED;
    }
    return POLAR_STATUS_OK;
}

static int is_power_of_two_int(int x) {
    return (x > 0) && ((x & (x - 1)) == 0);
}

static int ceil_log2_int(int x) {
    int n = 0;
    int v = 1;
    while (v < x) {
        v <<= 1;
        n++;
    }
    return n;
}

static int log2_int_exact(int x) {
    int n = 0;
    while (x > 1) {
        x >>= 1;
        n++;
    }
    return n;
}

static inline double sign_double(double x) {
    if (x > 0.0) {
        return 1.0;
    }
    if (x < 0.0) {
        return -1.0;
    }
    return 0.0;
}

static inline double minstar_scalar(double a, double b, int approx) {
    if (approx || isinf(fabs(a)) || isinf(fabs(b))) {
        return sign_double(a) * sign_double(b) * fmin(fabs(a), fabs(b));
    }
    return 2.0 * atanh(tanh(a / 2.0) * tanh(b / 2.0));
}

static inline double phi_scalar(double pm_prev, double llr, int u, int approx) {
    if (approx) {
        double hard = 0.5 * (1.0 - sign_double(llr));
        if (hard != (double)u) {
            return pm_prev + fabs(llr);
        }
        return pm_prev;
    }
    return pm_prev + log1p(exp(-(1.0 - 2.0 * (double)u) * llr));
}

static int cmp_pm_idx(const void *a, const void *b) {
    const pm_idx_t *pa = (const pm_idx_t *)a;
    const pm_idx_t *pb = (const pm_idx_t *)b;
    if (pa->pm < pb->pm) {
        return -1;
    }
    if (pa->pm > pb->pm) {
        return 1;
    }
    if (pa->idx < pb->idx) {
        return -1;
    }
    if (pa->idx > pb->idx) {
        return 1;
    }
    return 0;
}

static polar_status_t sorted_path_indices(const double *pm, int count, int **indices_out) {
    pm_idx_t *pairs = NULL;
    int *indices = NULL;
    polar_status_t st;
    int i;

    st = alloc_i32((size_t)count, &indices);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    if (count == 0) {
        *indices_out = indices;
        return POLAR_STATUS_OK;
    }
    pairs = (pm_idx_t *)malloc((size_t)count * sizeof(pm_idx_t));
    if (pairs == NULL) {
        free(indices);
        return POLAR_STATUS_ALLOCATION_FAILED;
    }

    for (i = 0; i < count; i++) {
        pairs[i].pm = pm[i];
        pairs[i].idx = i;
    }
    qsort(pairs, (size_t)count, sizeof(pm_idx_t), cmp_pm_idx);
    for (i = 0; i < count; i++) {
        indices[i] = pairs[i].idx;
    }

    free(pairs);
    *indices_out = indices;
    return POLAR_STATUS_OK;
}

static void copy_path_block_u8(uint8_t *dst, const uint8_t *src, int elem_per_path) {
    memcpy(dst, src, (size_t)elem_per_path);
}

static void copy_path_block_f64(double *dst, const double *src, int elem_per_path) {
    memcpy(dst, src, (size_t)elem_per_path * sizeof(double));
}

static void polar_transform_inplace(uint8_t *d, int N) {
    int stage;
    for (stage = 0; (1 << stage) < N; stage++) {
        int step = 1 << stage;
        int base;
        for (base = 0; base < N; base += 2 * step) {
            int j;
            for (j = 0; j < step; j++) {
                d[base + j] ^= d[base + step + j];
            }
        }
    }
}

static polar_status_t get_crc_generator_matrix(int A, const uint8_t *poly, int poly_len, u8_matrix_t *G) {
    int P = poly_len - 1;
    uint8_t *data = NULL;
    int k;

    if (P < 1) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    if (A < 0) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    if (A == 0) {
        G->data = NULL;
        G->rows = 0;
        G->cols = P;
        return POLAR_STATUS_OK;
    }

    if (alloc_u8((size_t)A * (size_t)P, &data) != POLAR_STATUS_OK) {
        return POLAR_STATUS_ALLOCATION_FAILED;
    }
    memset(data, 0, (size_t)A * (size_t)P);

    for (k = 0; k < P; k++) {
        data[(size_t)(A - 1) * (size_t)P + (size_t)k] = poly[k + 1];
    }

    for (k = A - 2; k >= 0; k--) {
        int j;
        uint8_t leading = data[(size_t)(k + 1) * (size_t)P + 0];
        for (j = 0; j < P; j++) {
            uint8_t shifted = (j + 1 < P) ? data[(size_t)(k + 1) * (size_t)P + (size_t)(j + 1)] : 0;
            uint8_t masked = leading ? poly[j + 1] : 0;
            data[(size_t)k * (size_t)P + (size_t)j] = (uint8_t)(shifted ^ masked);
        }
    }

    G->data = data;
    G->rows = A;
    G->cols = P;
    return POLAR_STATUS_OK;
}

static void free_u8_matrix(u8_matrix_t *M) {
    if (M != NULL) {
        free(M->data);
        M->data = NULL;
        M->rows = 0;
        M->cols = 0;
    }
}

static void vec_mul_matrix_mod2(const uint8_t *vec, int vec_len, const u8_matrix_t *M, uint8_t *out) {
    int c;
    for (c = 0; c < M->cols; c++) {
        uint8_t v = 0;
        int r;
        for (r = 0; r < vec_len; r++) {
            v ^= (uint8_t)(vec[r] & U8_AT((*M), r, c));
        }
        out[c] = (uint8_t)(v & 1u);
    }
}

static polar_status_t get_3gpp_n_impl(int K, int E, int n_max, int *N_out) {
    int n1;
    int n2;
    int n;

    if (K <= 0 || E <= 0 || N_out == NULL) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    if ((double)E <= (9.0 / 8.0) * pow(2.0, (double)ceil_log2_int(E) - 1.0) && ((double)K / (double)E) < (9.0 / 16.0)) {
        n1 = ceil_log2_int(E) - 1;
    } else {
        n1 = ceil_log2_int(E);
    }

    n2 = (int)ceil(log2((double)K / (1.0 / 8.0)));
    if (n_max < 5) {
        n_max = 5;
    }

    n = n1;
    if (n > n2) {
        n = n2;
    }
    if (n > n_max) {
        n = n_max;
    }
    if (n < 5) {
        n = 5;
    }

    *N_out = 1 << n;
    return POLAR_STATUS_OK;
}

static polar_status_t get_3gpp_rate_matching_pattern_impl(int K, int N, int E, int **pattern_out, int *mode_out) {
    int *pattern = NULL;
    int *y = NULL;
    int n32;
    int k;
    polar_status_t st;

    if (pattern_out == NULL || mode_out == NULL) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }
    if (!is_power_of_two_int(N) || N < 32 || E <= 0) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    st = alloc_i32((size_t)E, &pattern);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    st = alloc_i32((size_t)N, &y);
    if (st != POLAR_STATUS_OK) {
        free(pattern);
        return st;
    }

    n32 = N / 32;
    for (k = 0; k < N; k++) {
        int i = (32 * k) / N;
        int j = P32[i] * n32 + (k % n32);
        y[k] = j;
    }

    if (E >= N) {
        for (k = 0; k < E; k++) {
            pattern[k] = y[k % N];
        }
        *mode_out = RM_MODE_REPETITION;
    } else {
        if (((double)K / (double)E) <= (7.0 / 16.0)) {
            for (k = 0; k < E; k++) {
                pattern[k] = y[k + N - E];
            }
            *mode_out = RM_MODE_PUNCTURING;
        } else {
            for (k = 0; k < E; k++) {
                pattern[k] = y[k];
            }
            *mode_out = RM_MODE_SHORTENING;
        }
    }

    free(y);
    *pattern_out = pattern;
    return POLAR_STATUS_OK;
}

static polar_status_t get_bivs_rate_matching_pattern_impl(int N, int E, int **pattern_out, int *mode_out) {
    int *pattern = NULL;
    int n;
    int idx;
    polar_status_t st;

    if (!is_power_of_two_int(N) || E > N || E < 0 || pattern_out == NULL || mode_out == NULL) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }
    n = log2_int_exact(N);

    st = alloc_i32((size_t)E, &pattern);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    for (idx = 0; idx < E; idx++) {
        int x = idx;
        int rev = 0;
        int b;
        for (b = 0; b < n; b++) {
            rev = (rev << 1) | (x & 1);
            x >>= 1;
        }
        pattern[idx] = rev;
    }

    /* Sort ascending (small E in this code path, so O(E^2) is fine). */
    {
        int i;
        for (i = 0; i < E; i++) {
            int j;
            for (j = i + 1; j < E; j++) {
                if (pattern[j] < pattern[i]) {
                    int t = pattern[i];
                    pattern[i] = pattern[j];
                    pattern[j] = t;
                }
            }
        }
    }

    *mode_out = RM_MODE_SHORTENING;
    *pattern_out = pattern;
    return POLAR_STATUS_OK;
}

static polar_status_t get_bivp_rate_matching_pattern_impl(int N, int E, int **pattern_out, int *mode_out) {
    int *pattern = NULL;
    int n;
    int idx;
    polar_status_t st;

    if (!is_power_of_two_int(N) || E > N || E < 0 || pattern_out == NULL || mode_out == NULL) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }
    n = log2_int_exact(N);

    st = alloc_i32((size_t)E, &pattern);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    for (idx = 0; idx < E; idx++) {
        int x = (N - E) + idx;
        int rev = 0;
        int b;
        for (b = 0; b < n; b++) {
            rev = (rev << 1) | (x & 1);
            x >>= 1;
        }
        pattern[idx] = rev;
    }

    {
        int i;
        for (i = 0; i < E; i++) {
            int j;
            for (j = i + 1; j < E; j++) {
                if (pattern[j] < pattern[i]) {
                    int t = pattern[i];
                    pattern[i] = pattern[j];
                    pattern[j] = t;
                }
            }
        }
    }

    *mode_out = RM_MODE_PUNCTURING;
    *pattern_out = pattern;
    return POLAR_STATUS_OK;
}

static polar_status_t get_natp_rate_matching_pattern_impl(int N, int E, int **pattern_out, int *mode_out) {
    int *pattern = NULL;
    int i;
    polar_status_t st;

    if (!is_power_of_two_int(N) || E > N || E < 0 || pattern_out == NULL || mode_out == NULL) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    st = alloc_i32((size_t)E, &pattern);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    for (i = 0; i < E; i++) {
        pattern[i] = (N - E) + i;
    }

    *mode_out = RM_MODE_PUNCTURING;
    *pattern_out = pattern;
    return POLAR_STATUS_OK;
}

static polar_status_t get_nats_rate_matching_pattern_impl(int N, int E, int **pattern_out, int *mode_out) {
    int *pattern = NULL;
    int i;
    polar_status_t st;

    if (!is_power_of_two_int(N) || E > N || E < 0 || pattern_out == NULL || mode_out == NULL) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    st = alloc_i32((size_t)E, &pattern);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    for (i = 0; i < E; i++) {
        pattern[i] = i;
    }

    *mode_out = RM_MODE_SHORTENING;
    *pattern_out = pattern;
    return POLAR_STATUS_OK;
}

static polar_status_t get_3gpp_channel_interleaver_pattern_impl(int E, int **pattern_out) {
    int T = 0;
    int *pattern = NULL;
    int *v = NULL;
    int i;
    int j;
    int k;
    polar_status_t st;

    if (E <= 0 || pattern_out == NULL) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    while ((T * (T + 1)) / 2 < E) {
        T++;
    }

    st = alloc_i32((size_t)E, &pattern);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    st = alloc_i32((size_t)T * (size_t)T, &v);
    if (st != POLAR_STATUS_OK) {
        free(pattern);
        return st;
    }

    for (i = 0; i < T * T; i++) {
        v[i] = -1;
    }

    k = 0;
    for (i = 0; i < T; i++) {
        for (j = 0; j < T - i; j++) {
            if (k < E) {
                v[i * T + j] = k;
            }
            k++;
        }
    }

    k = 0;
    for (j = 0; j < T; j++) {
        for (i = 0; i < T - j; i++) {
            if (v[i * T + j] >= 0) {
                pattern[k++] = v[i * T + j];
            }
        }
    }

    free(v);
    *pattern_out = pattern;
    return POLAR_STATUS_OK;
}

static polar_status_t get_3gpp_crc_interleaver_pattern_impl(int K, int **pattern_out) {
    int *pattern = NULL;
    int k = 0;
    int m;
    int len = (int)(sizeof(PI_IL_MAX) / sizeof(PI_IL_MAX[0]));
    polar_status_t st;

    if (K < 0 || K > len || pattern_out == NULL) {
        return POLAR_STATUS_UNSUPPORTED_BLOCK_LENGTH;
    }

    st = alloc_i32((size_t)K, &pattern);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    for (m = 0; m < len; m++) {
        if (PI_IL_MAX[m] >= len - K) {
            pattern[k++] = PI_IL_MAX[m] - (len - K);
        }
    }

    *pattern_out = pattern;
    return POLAR_STATUS_OK;
}

static polar_status_t get_3gpp_sequence_pattern_impl(int N, int **q_out) {
    int *q = NULL;
    int i;
    int k = 0;
    int len = (int)(sizeof(Q_NMAX) / sizeof(Q_NMAX[0]));
    polar_status_t st;

    if (!is_power_of_two_int(N) || q_out == NULL) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }
    if (N > len) {
        return POLAR_STATUS_UNSUPPORTED_BLOCK_LENGTH;
    }

    st = alloc_i32((size_t)N, &q);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    for (i = 0; i < len; i++) {
        if (Q_NMAX[i] < N) {
            q[k++] = Q_NMAX[i];
        }
    }

    *q_out = q;
    return POLAR_STATUS_OK;
}

static int popcount_u32(uint32_t x) {
    int c = 0;
    while (x != 0u) {
        x &= (x - 1u);
        c++;
    }
    return c;
}

static polar_status_t get_pw_sequence_pattern_impl(int N, int **q_out) {
    double *w = NULL;
    int *idx = NULL;
    int n;
    int i;
    polar_status_t st;

    if (!is_power_of_two_int(N) || q_out == NULL) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    n = log2_int_exact(N);

    st = alloc_f64((size_t)N, &w);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    st = alloc_i32((size_t)N, &idx);
    if (st != POLAR_STATUS_OK) {
        free(w);
        return st;
    }

    for (i = 0; i < N; i++) {
        int b;
        double wi = 0.0;
        for (b = 0; b < n; b++) {
            int bit = (i >> b) & 1;
            if (bit) {
                wi += pow(2.0, (double)b / 4.0);
            }
        }
        w[i] = wi;
        idx[i] = i;
    }

    for (i = 0; i < N; i++) {
        int j;
        for (j = i + 1; j < N; j++) {
            if (w[idx[j]] < w[idx[i]] || (w[idx[j]] == w[idx[i]] && idx[j] < idx[i])) {
                int t = idx[i];
                idx[i] = idx[j];
                idx[j] = t;
            }
        }
    }

    free(w);
    *q_out = idx;
    return POLAR_STATUS_OK;
}

static polar_status_t mode_compat_check(int E, int N, int mode) {
    if (mode == RM_MODE_REPETITION) {
        if (E < N) {
            return POLAR_STATUS_INVALID_ARGUMENT;
        }
    } else if (mode == RM_MODE_PUNCTURING || mode == RM_MODE_SHORTENING) {
        if (E >= N) {
            return POLAR_STATUS_INVALID_ARGUMENT;
        }
    } else {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }
    return POLAR_STATUS_OK;
}

static polar_status_t get_3gpp_info_bit_pattern_impl(
    int I,
    const int *Q_N,
    int N,
    const int *rate_matching_pattern,
    int E,
    int mode,
    uint8_t **info_out
) {
    uint8_t *info = NULL;
    uint8_t *mask_frozen = NULL;
    int *q_itmp = NULL;
    int q_itmp_len = 0;
    int i;
    polar_status_t st;

    if (I < 0 || Q_N == NULL || rate_matching_pattern == NULL || info_out == NULL) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }
    if (!is_power_of_two_int(N)) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }
    if (I > N || I > E) {
        return POLAR_STATUS_UNSUPPORTED_BLOCK_LENGTH;
    }
    for (i = 0; i < E; i++) {
        if (rate_matching_pattern[i] < 0 || rate_matching_pattern[i] >= N) {
            return POLAR_STATUS_INVALID_ARGUMENT;
        }
    }
    st = mode_compat_check(E, N, mode);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    st = calloc_u8((size_t)N, &mask_frozen);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    {
        uint8_t *in_rate = NULL;
        st = calloc_u8((size_t)N, &in_rate);
        if (st != POLAR_STATUS_OK) {
            free(mask_frozen);
            return st;
        }
        for (i = 0; i < E; i++) {
            in_rate[rate_matching_pattern[i]] = 1;
        }
        for (i = 0; i < N; i++) {
            if (!in_rate[i]) {
                mask_frozen[i] = 1;
            }
        }
        free(in_rate);
    }

    if (mode == RM_MODE_PUNCTURING) {
        int extra_end;
        if (E >= (3 * N) / 4) {
            extra_end = (int)ceil(3.0 * (double)N / 4.0 - (double)E / 2.0);
        } else {
            extra_end = (int)ceil(9.0 * (double)N / 16.0 - (double)E / 4.0);
        }
        for (i = 0; i < extra_end && i < N; i++) {
            if (i >= 0) {
                mask_frozen[i] = 1;
            }
        }
    }

    st = alloc_i32((size_t)N, &q_itmp);
    if (st != POLAR_STATUS_OK) {
        free(mask_frozen);
        return st;
    }

    for (i = 0; i < N; i++) {
        int q = Q_N[i];
        if (q < 0 || q >= N) {
            free(mask_frozen);
            free(q_itmp);
            return POLAR_STATUS_INVALID_ARGUMENT;
        }
        if (!mask_frozen[q]) {
            q_itmp[q_itmp_len++] = q;
        }
    }

    if (q_itmp_len < I) {
        free(mask_frozen);
        free(q_itmp);
        return POLAR_STATUS_UNSUPPORTED_BLOCK_LENGTH;
    }

    st = calloc_u8((size_t)N, &info);
    if (st != POLAR_STATUS_OK) {
        free(mask_frozen);
        free(q_itmp);
        return st;
    }

    for (i = q_itmp_len - I; i < q_itmp_len; i++) {
        info[q_itmp[i]] = 1;
    }

    free(mask_frozen);
    free(q_itmp);
    *info_out = info;
    return POLAR_STATUS_OK;
}

static polar_status_t get_info_bit_pattern_impl(
    int I,
    const int *Q_N,
    int N,
    const int *rate_matching_pattern,
    int E,
    uint8_t **info_out
) {
    uint8_t *in_rate = NULL;
    int *q_itmp = NULL;
    int q_itmp_len = 0;
    uint8_t *info = NULL;
    int i;
    polar_status_t st;

    if (I < 0 || Q_N == NULL || rate_matching_pattern == NULL || info_out == NULL) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }
    if (!is_power_of_two_int(N)) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }
    if (I > N || I > E) {
        return POLAR_STATUS_UNSUPPORTED_BLOCK_LENGTH;
    }

    st = calloc_u8((size_t)N, &in_rate);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    for (i = 0; i < E; i++) {
        int r = rate_matching_pattern[i];
        if (r < 0 || r >= N) {
            free(in_rate);
            return POLAR_STATUS_INVALID_ARGUMENT;
        }
        in_rate[r] = 1;
    }

    st = alloc_i32((size_t)N, &q_itmp);
    if (st != POLAR_STATUS_OK) {
        free(in_rate);
        return st;
    }
    for (i = 0; i < N; i++) {
        int q = Q_N[i];
        if (q < 0 || q >= N) {
            free(in_rate);
            free(q_itmp);
            return POLAR_STATUS_INVALID_ARGUMENT;
        }
        if (in_rate[q]) {
            q_itmp[q_itmp_len++] = q;
        }
    }

    if (q_itmp_len < I) {
        free(in_rate);
        free(q_itmp);
        return POLAR_STATUS_UNSUPPORTED_BLOCK_LENGTH;
    }

    st = calloc_u8((size_t)N, &info);
    if (st != POLAR_STATUS_OK) {
        free(in_rate);
        free(q_itmp);
        return st;
    }

    for (i = q_itmp_len - I; i < q_itmp_len; i++) {
        info[q_itmp[i]] = 1;
    }

    free(in_rate);
    free(q_itmp);
    *info_out = info;
    return POLAR_STATUS_OK;
}

static polar_status_t get_pc_bit_pattern_impl(
    const uint8_t *info_bit_pattern,
    const int *Q_N,
    int N,
    int n_PC,
    int n_PC_wm,
    uint8_t **pc_out
) {
    int I = 0;
    int *q_n_i = NULL;
    int q_n_i_len = 0;
    int *q_tilde_flip = NULL;
    int q_tilde_len = 0;
    int *pc_positions = NULL;
    uint8_t *pc = NULL;
    int i;
    polar_status_t st;

    if (info_bit_pattern == NULL || Q_N == NULL || pc_out == NULL) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }
    if (!is_power_of_two_int(N)) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    for (i = 0; i < N; i++) {
        if (info_bit_pattern[i]) {
            I++;
        }
    }

    if (n_PC > I || n_PC_wm > n_PC) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    st = alloc_i32((size_t)N, &q_n_i);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    for (i = 0; i < N; i++) {
        int q = Q_N[i];
        if (q < 0 || q >= N) {
            free(q_n_i);
            return POLAR_STATUS_INVALID_ARGUMENT;
        }
        if (info_bit_pattern[q]) {
            q_n_i[q_n_i_len++] = q;
        }
    }

    q_tilde_len = q_n_i_len - n_PC;
    if (q_tilde_len < 0) {
        q_tilde_len = 0;
    }

    st = alloc_i32((size_t)q_tilde_len, &q_tilde_flip);
    if (st != POLAR_STATUS_OK) {
        free(q_n_i);
        return st;
    }

    for (i = 0; i < q_tilde_len; i++) {
        q_tilde_flip[i] = q_n_i[q_n_i_len - 1 - i];
    }

    if (n_PC_wm > 0) {
        int *weights = NULL;
        int *indices = NULL;

        st = alloc_i32((size_t)q_tilde_len, &weights);
        if (st != POLAR_STATUS_OK) {
            free(q_n_i);
            free(q_tilde_flip);
            return st;
        }
        st = alloc_i32((size_t)q_tilde_len, &indices);
        if (st != POLAR_STATUS_OK) {
            free(q_n_i);
            free(q_tilde_flip);
            free(weights);
            return st;
        }

        for (i = 0; i < q_tilde_len; i++) {
            int row = q_tilde_flip[i];
            weights[i] = 1 << popcount_u32((uint32_t)row);
            indices[i] = i;
        }

        for (i = 0; i < q_tilde_len; i++) {
            int j;
            for (j = i + 1; j < q_tilde_len; j++) {
                int ai = indices[i];
                int aj = indices[j];
                if (weights[aj] < weights[ai] || (weights[aj] == weights[ai] && aj < ai)) {
                    int t = indices[i];
                    indices[i] = indices[j];
                    indices[j] = t;
                }
            }
        }

        st = alloc_i32((size_t)n_PC, &pc_positions);
        if (st != POLAR_STATUS_OK) {
            free(q_n_i);
            free(q_tilde_flip);
            free(weights);
            free(indices);
            return st;
        }

        for (i = 0; i < (n_PC - n_PC_wm); i++) {
            pc_positions[i] = q_n_i[i];
        }
        for (i = 0; i < n_PC_wm; i++) {
            pc_positions[n_PC - n_PC_wm + i] = q_tilde_flip[indices[i]];
        }

        free(weights);
        free(indices);
    } else {
        st = alloc_i32((size_t)n_PC, &pc_positions);
        if (st != POLAR_STATUS_OK) {
            free(q_n_i);
            free(q_tilde_flip);
            return st;
        }
        for (i = 0; i < n_PC; i++) {
            pc_positions[i] = q_n_i[i];
        }
    }

    st = calloc_u8((size_t)N, &pc);
    if (st != POLAR_STATUS_OK) {
        free(q_n_i);
        free(q_tilde_flip);
        free(pc_positions);
        return st;
    }

    for (i = 0; i < n_PC; i++) {
        pc[pc_positions[i]] = 1;
    }

    free(q_n_i);
    free(q_tilde_flip);
    free(pc_positions);
    *pc_out = pc;
    return POLAR_STATUS_OK;
}

static polar_status_t rate_dematch(
    const double *e_tilde,
    int E,
    int N,
    const int *rate_matching_pattern,
    int mode,
    double **d_tilde_out
) {
    double *d_tilde = NULL;
    int i;
    polar_status_t st;

    st = alloc_f64((size_t)N, &d_tilde);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    if (mode == RM_MODE_REPETITION) {
        for (i = 0; i < N; i++) {
            d_tilde[i] = 0.0;
        }
        for (i = 0; i < E; i++) {
            d_tilde[rate_matching_pattern[i]] += e_tilde[i];
        }
    } else if (mode == RM_MODE_PUNCTURING) {
        for (i = 0; i < N; i++) {
            d_tilde[i] = 0.0;
        }
        for (i = 0; i < E; i++) {
            d_tilde[rate_matching_pattern[i]] = e_tilde[i];
        }
    } else if (mode == RM_MODE_SHORTENING) {
        for (i = 0; i < N; i++) {
            d_tilde[i] = INFINITY;
        }
        for (i = 0; i < E; i++) {
            d_tilde[rate_matching_pattern[i]] = e_tilde[i];
        }
    } else {
        free(d_tilde);
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    *d_tilde_out = d_tilde;
    return POLAR_STATUS_OK;
}

static polar_status_t extract_info_bits(
    const uint8_t *u_hat,
    int N,
    const uint8_t *info_pattern,
    uint8_t **out_bits,
    int *out_len
) {
    uint8_t *bits;
    int count = 0;
    int i;
    int k = 0;
    polar_status_t st;

    for (i = 0; i < N; i++) {
        if (info_pattern[i]) {
            count++;
        }
    }

    st = alloc_u8((size_t)count, &bits);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    for (i = 0; i < N; i++) {
        if (info_pattern[i]) {
            bits[k++] = u_hat[i];
        }
    }

    *out_bits = bits;
    *out_len = count;
    return POLAR_STATUS_OK;
}

static polar_status_t extract_info_bits_from_scl(
    const scl_ctx_t *ctx,
    int path,
    const uint8_t *info_pattern,
    uint8_t **out_bits,
    int *out_len
) {
    uint8_t *bits;
    int count = 0;
    int i;
    int k = 0;
    polar_status_t st;

    for (i = 0; i < ctx->N; i++) {
        if (info_pattern[i]) {
            count++;
        }
    }

    st = alloc_u8((size_t)count, &bits);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    for (i = 0; i < ctx->N; i++) {
        if (info_pattern[i]) {
            bits[k++] = ctx->bits[SCL_BITS_IDX(ctx, path, i, 0)];
        }
    }

    *out_bits = bits;
    *out_len = count;
    return POLAR_STATUS_OK;
}

static polar_status_t polar_kernel_rate_match(
    const uint8_t *u,
    int N,
    const int *rate_matching_pattern,
    int E,
    uint8_t **e_out
) {
    uint8_t *d = NULL;
    uint8_t *e = NULL;
    polar_status_t st;
    int i;

    if (u == NULL || rate_matching_pattern == NULL || e_out == NULL || !is_power_of_two_int(N) || E < 0) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    st = alloc_u8((size_t)N, &d);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    memcpy(d, u, (size_t)N);
    polar_transform_inplace(d, N);

    st = alloc_u8((size_t)E, &e);
    if (st != POLAR_STATUS_OK) {
        free(d);
        return st;
    }

    for (i = 0; i < E; i++) {
        int idx = rate_matching_pattern[i];
        if (idx < 0 || idx >= N) {
            free(d);
            free(e);
            return POLAR_STATUS_INVALID_ARGUMENT;
        }
        e[i] = d[idx];
    }

    free(d);
    *e_out = e;
    return POLAR_STATUS_OK;
}

static polar_status_t polar_encoder_impl(
    const uint8_t *a,
    int A,
    const uint8_t *info_bit_pattern,
    int N,
    const int *rate_matching_pattern,
    int E,
    uint8_t **e_out
) {
    uint8_t *u = NULL;
    int i;
    int k = 0;
    polar_status_t st;

    if (!is_power_of_two_int(N)) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }
    if (a == NULL || info_bit_pattern == NULL || rate_matching_pattern == NULL || e_out == NULL) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    st = calloc_u8((size_t)N, &u);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    for (i = 0; i < N; i++) {
        if (info_bit_pattern[i]) {
            if (k >= A) {
                free(u);
                return POLAR_STATUS_INVALID_ARGUMENT;
            }
            u[i] = a[k++];
        }
    }
    if (k != A) {
        free(u);
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    st = polar_kernel_rate_match(u, N, rate_matching_pattern, E, e_out);
    free(u);
    return st;
}

static polar_status_t ca_polar_encoder_impl(
    const uint8_t *a,
    int A,
    const uint8_t *crc_poly,
    int crc_poly_len,
    const uint8_t *info_bit_pattern,
    int N,
    const int *rate_matching_pattern,
    int E,
    uint8_t **e_out
) {
    int P = crc_poly_len - 1;
    int K = A + P;
    u8_matrix_t Gp = {0};
    uint8_t *crc = NULL;
    uint8_t *b = NULL;
    polar_status_t st;

    st = get_crc_generator_matrix(A, crc_poly, crc_poly_len, &Gp);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    st = alloc_u8((size_t)P, &crc);
    if (st != POLAR_STATUS_OK) {
        free_u8_matrix(&Gp);
        return st;
    }
    vec_mul_matrix_mod2(a, A, &Gp, crc);

    st = alloc_u8((size_t)K, &b);
    if (st != POLAR_STATUS_OK) {
        free_u8_matrix(&Gp);
        free(crc);
        return st;
    }
    memcpy(b, a, (size_t)A);
    memcpy(b + A, crc, (size_t)P);

    st = polar_encoder_impl(b, K, info_bit_pattern, N, rate_matching_pattern, E, e_out);

    free_u8_matrix(&Gp);
    free(crc);
    free(b);
    return st;
}

static polar_status_t dca_polar_encoder_impl(
    const uint8_t *a,
    int A,
    const uint8_t *crc_poly,
    int crc_poly_len,
    const int *crc_interleaver,
    int K,
    const uint8_t *info_bit_pattern,
    int N,
    const int *rate_matching_pattern,
    int E,
    uint8_t **e_out
) {
    int P = crc_poly_len - 1;
    u8_matrix_t Gp = {0};
    uint8_t *crc = NULL;
    uint8_t *b = NULL;
    uint8_t *c = NULL;
    polar_status_t st;
    int i;

    if (A + P != K) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    st = get_crc_generator_matrix(A, crc_poly, crc_poly_len, &Gp);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    st = alloc_u8((size_t)P, &crc);
    if (st != POLAR_STATUS_OK) {
        free_u8_matrix(&Gp);
        return st;
    }
    vec_mul_matrix_mod2(a, A, &Gp, crc);

    st = alloc_u8((size_t)K, &b);
    if (st != POLAR_STATUS_OK) {
        free_u8_matrix(&Gp);
        free(crc);
        return st;
    }
    memcpy(b, a, (size_t)A);
    memcpy(b + A, crc, (size_t)P);

    st = alloc_u8((size_t)K, &c);
    if (st != POLAR_STATUS_OK) {
        free_u8_matrix(&Gp);
        free(crc);
        free(b);
        return st;
    }

    for (i = 0; i < K; i++) {
        int idx = crc_interleaver[i];
        if (idx < 0 || idx >= K) {
            free_u8_matrix(&Gp);
            free(crc);
            free(b);
            free(c);
            return POLAR_STATUS_INVALID_ARGUMENT;
        }
        c[i] = b[idx];
    }

    st = polar_encoder_impl(c, K, info_bit_pattern, N, rate_matching_pattern, E, e_out);

    free_u8_matrix(&Gp);
    free(crc);
    free(b);
    free(c);
    return st;
}

static polar_status_t ds1ca_polar_encoder_impl(
    const uint8_t *a,
    int A,
    const uint8_t *crc_poly,
    int crc_poly_len,
    const uint8_t *crc_scrambling_pattern,
    int scrambling_len,
    const int *crc_interleaver,
    int K,
    const uint8_t *info_bit_pattern,
    int N,
    const int *rate_matching_pattern,
    int E,
    uint8_t **e_out
) {
    int P = crc_poly_len - 1;
    u8_matrix_t Gp = {0};
    uint8_t *in = NULL;
    uint8_t *crc = NULL;
    uint8_t *b = NULL;
    uint8_t *c = NULL;
    polar_status_t st;
    int i;

    if (A + P != K) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }
    if (P < scrambling_len) {
        return POLAR_STATUS_UNSUPPORTED_BLOCK_LENGTH;
    }

    st = get_crc_generator_matrix(A + P, crc_poly, crc_poly_len, &Gp);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    st = alloc_u8((size_t)(A + P), &in);
    if (st != POLAR_STATUS_OK) {
        free_u8_matrix(&Gp);
        return st;
    }
    for (i = 0; i < P; i++) {
        in[i] = 1;
    }
    memcpy(in + P, a, (size_t)A);

    st = alloc_u8((size_t)P, &crc);
    if (st != POLAR_STATUS_OK) {
        free_u8_matrix(&Gp);
        free(in);
        return st;
    }
    vec_mul_matrix_mod2(in, A + P, &Gp, crc);

    for (i = 0; i < scrambling_len; i++) {
        crc[P - scrambling_len + i] ^= crc_scrambling_pattern[i];
    }

    st = alloc_u8((size_t)K, &b);
    if (st != POLAR_STATUS_OK) {
        free_u8_matrix(&Gp);
        free(in);
        free(crc);
        return st;
    }
    memcpy(b, a, (size_t)A);
    memcpy(b + A, crc, (size_t)P);

    st = alloc_u8((size_t)K, &c);
    if (st != POLAR_STATUS_OK) {
        free_u8_matrix(&Gp);
        free(in);
        free(crc);
        free(b);
        return st;
    }
    for (i = 0; i < K; i++) {
        int idx = crc_interleaver[i];
        if (idx < 0 || idx >= K) {
            free_u8_matrix(&Gp);
            free(in);
            free(crc);
            free(b);
            free(c);
            return POLAR_STATUS_INVALID_ARGUMENT;
        }
        c[i] = b[idx];
    }

    st = polar_encoder_impl(c, K, info_bit_pattern, N, rate_matching_pattern, E, e_out);

    free_u8_matrix(&Gp);
    free(in);
    free(crc);
    free(b);
    free(c);
    return st;
}

static polar_status_t pc_polar_encoder_impl(
    const uint8_t *a,
    int A,
    const uint8_t *info_bit_pattern,
    const uint8_t *pc_bit_pattern,
    int N,
    int pc_buf_len,
    const int *rate_matching_pattern,
    int E,
    uint8_t **e_out
) {
    uint8_t *u = NULL;
    uint8_t *y = NULL;
    int k = 0;
    int n;
    polar_status_t st;

    st = calloc_u8((size_t)N, &u);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    st = calloc_u8((size_t)pc_buf_len, &y);
    if (st != POLAR_STATUS_OK) {
        free(u);
        return st;
    }

    for (n = 0; n < N; n++) {
        uint8_t first = y[0];
        memmove(y, y + 1, (size_t)(pc_buf_len - 1));
        y[pc_buf_len - 1] = first;

        if (info_bit_pattern[n]) {
            if (pc_bit_pattern[n]) {
                u[n] = y[0];
            } else {
                if (k >= A) {
                    free(u);
                    free(y);
                    return POLAR_STATUS_INVALID_ARGUMENT;
                }
                u[n] = a[k++];
                y[0] ^= u[n];
            }
        }
    }

    st = polar_kernel_rate_match(u, N, rate_matching_pattern, E, e_out);

    free(u);
    free(y);
    return st;
}

static polar_status_t pcca_polar_encoder_impl(
    const uint8_t *a,
    int A,
    const uint8_t *crc_poly,
    int crc_poly_len,
    const uint8_t *info_bit_pattern,
    const uint8_t *pc_bit_pattern,
    int N,
    int pc_buf_len,
    const int *rate_matching_pattern,
    int E,
    uint8_t **e_out
) {
    int P = crc_poly_len - 1;
    int K = A + P;
    u8_matrix_t Gp = {0};
    uint8_t *crc = NULL;
    uint8_t *b = NULL;
    uint8_t *u = NULL;
    uint8_t *y = NULL;
    int k = 0;
    int n;
    polar_status_t st;

    st = get_crc_generator_matrix(A, crc_poly, crc_poly_len, &Gp);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    st = alloc_u8((size_t)P, &crc);
    if (st != POLAR_STATUS_OK) {
        free_u8_matrix(&Gp);
        return st;
    }
    vec_mul_matrix_mod2(a, A, &Gp, crc);

    st = alloc_u8((size_t)K, &b);
    if (st != POLAR_STATUS_OK) {
        free_u8_matrix(&Gp);
        free(crc);
        return st;
    }
    memcpy(b, a, (size_t)A);
    memcpy(b + A, crc, (size_t)P);

    st = calloc_u8((size_t)N, &u);
    if (st != POLAR_STATUS_OK) {
        free_u8_matrix(&Gp);
        free(crc);
        free(b);
        return st;
    }
    st = calloc_u8((size_t)pc_buf_len, &y);
    if (st != POLAR_STATUS_OK) {
        free_u8_matrix(&Gp);
        free(crc);
        free(b);
        free(u);
        return st;
    }

    for (n = 0; n < N; n++) {
        uint8_t first = y[0];
        memmove(y, y + 1, (size_t)(pc_buf_len - 1));
        y[pc_buf_len - 1] = first;

        if (info_bit_pattern[n]) {
            if (pc_bit_pattern[n]) {
                u[n] = y[0];
            } else {
                if (k >= K) {
                    free_u8_matrix(&Gp);
                    free(crc);
                    free(b);
                    free(u);
                    free(y);
                    return POLAR_STATUS_INVALID_ARGUMENT;
                }
                u[n] = b[k++];
                y[0] ^= u[n];
            }
        }
    }

    st = polar_kernel_rate_match(u, N, rate_matching_pattern, E, e_out);

    free_u8_matrix(&Gp);
    free(crc);
    free(b);
    free(u);
    free(y);
    return st;
}

static polar_status_t scl_ctx_init(
    scl_ctx_t *ctx,
    int N,
    int L,
    int min_sum,
    const uint8_t *info_bit_pattern,
    const double *d_tilde
) {
    int r;
    size_t bits_count;
    size_t llr_count;

    memset(ctx, 0, sizeof(*ctx));

    if (!is_power_of_two_int(N) || L <= 0 || info_bit_pattern == NULL || d_tilde == NULL) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    ctx->N = N;
    ctx->n = log2_int_exact(N);
    ctx->n_cols = ctx->n + 1;
    ctx->L_target = L;
    ctx->L_capacity = (2 * L > 2) ? (2 * L) : 2;
    ctx->L_prime = 1;
    ctx->min_sum = min_sum ? 1 : 0;
    ctx->path_stride = N * ctx->n_cols;

    bits_count = (size_t)ctx->L_capacity * (size_t)ctx->path_stride;
    llr_count = bits_count;

    if (calloc_u8(bits_count, &ctx->bits) != POLAR_STATUS_OK) {
        return POLAR_STATUS_ALLOCATION_FAILED;
    }
    if (calloc_u8(bits_count, &ctx->bits_alt) != POLAR_STATUS_OK) {
        return POLAR_STATUS_ALLOCATION_FAILED;
    }
    if (calloc_u8((size_t)ctx->N * (size_t)ctx->n_cols, &ctx->bits_updated) != POLAR_STATUS_OK) {
        return POLAR_STATUS_ALLOCATION_FAILED;
    }

    if (calloc_f64(llr_count, &ctx->llrs) != POLAR_STATUS_OK) {
        return POLAR_STATUS_ALLOCATION_FAILED;
    }
    if (calloc_f64(llr_count, &ctx->llrs_alt) != POLAR_STATUS_OK) {
        return POLAR_STATUS_ALLOCATION_FAILED;
    }
    if (calloc_u8((size_t)ctx->N * (size_t)ctx->n_cols, &ctx->llrs_updated) != POLAR_STATUS_OK) {
        return POLAR_STATUS_ALLOCATION_FAILED;
    }

    if (calloc_f64((size_t)ctx->L_capacity, &ctx->PM) != POLAR_STATUS_OK) {
        return POLAR_STATUS_ALLOCATION_FAILED;
    }
    if (calloc_f64((size_t)ctx->L_capacity, &ctx->PM_alt) != POLAR_STATUS_OK) {
        return POLAR_STATUS_ALLOCATION_FAILED;
    }

    for (r = 0; r < ctx->N; r++) {
        ctx->bits_updated[(size_t)r * (size_t)ctx->n_cols + 0] = (uint8_t)(info_bit_pattern[r] ? 0 : 1);
        ctx->llrs[SCL_LLR_IDX(ctx, 0, r, ctx->n)] = d_tilde[r];
        ctx->llrs_updated[(size_t)r * (size_t)ctx->n_cols + (size_t)ctx->n] = 1;
    }
    ctx->PM[0] = 0.0;

    return POLAR_STATUS_OK;
}

static void scl_ctx_free(scl_ctx_t *ctx) {
    if (ctx == NULL) {
        return;
    }
    free(ctx->bits);
    free(ctx->bits_alt);
    free(ctx->bits_updated);
    free(ctx->llrs);
    free(ctx->llrs_alt);
    free(ctx->llrs_updated);
    free(ctx->PM);
    free(ctx->PM_alt);
    memset(ctx, 0, sizeof(*ctx));
}

static void update_bit_rec(scl_ctx_t *ctx, int row, int col);
static void update_llr_rec(scl_ctx_t *ctx, int row, int col);

static void update_bit_rec(scl_ctx_t *ctx, int row, int col) {
    int offset;
    int l;

    if (col <= 0) {
        return;
    }

    offset = ctx->N >> (ctx->n_cols - col);

    for (l = 0; l < ctx->L_prime; l++) {
        if ((row % (2 * offset)) >= offset) {
            if (!ctx->bits_updated[(size_t)row * (size_t)ctx->n_cols + (size_t)(col - 1)]) {
                update_bit_rec(ctx, row, col - 1);
            }
            ctx->bits[SCL_BITS_IDX(ctx, l, row, col)] = ctx->bits[SCL_BITS_IDX(ctx, l, row, col - 1)];
        } else {
            if (!ctx->bits_updated[(size_t)row * (size_t)ctx->n_cols + (size_t)(col - 1)]) {
                update_bit_rec(ctx, row, col - 1);
            }
            if (!ctx->bits_updated[(size_t)(row + offset) * (size_t)ctx->n_cols + (size_t)(col - 1)]) {
                update_bit_rec(ctx, row + offset, col - 1);
            }
            ctx->bits[SCL_BITS_IDX(ctx, l, row, col)] =
                (uint8_t)((ctx->bits[SCL_BITS_IDX(ctx, l, row, col - 1)] + ctx->bits[SCL_BITS_IDX(ctx, l, row + offset, col - 1)]) & 1u);
        }
    }

    ctx->bits_updated[(size_t)row * (size_t)ctx->n_cols + (size_t)col] = 1;
}

static void update_llr_rec(scl_ctx_t *ctx, int row, int col) {
    int offset;
    int l;

    if (col >= ctx->n) {
        return;
    }

    offset = ctx->N >> (ctx->n_cols - col - 1);

    for (l = 0; l < ctx->L_prime; l++) {
        if ((row % (2 * offset)) >= offset) {
            if (!ctx->bits_updated[(size_t)(row - offset) * (size_t)ctx->n_cols + (size_t)col]) {
                update_bit_rec(ctx, row - offset, col);
            }
            if (!ctx->llrs_updated[(size_t)(row - offset) * (size_t)ctx->n_cols + (size_t)(col + 1)]) {
                update_llr_rec(ctx, row - offset, col + 1);
            }
            if (!ctx->llrs_updated[(size_t)row * (size_t)ctx->n_cols + (size_t)(col + 1)]) {
                update_llr_rec(ctx, row, col + 1);
            }
            ctx->llrs[SCL_LLR_IDX(ctx, l, row, col)] =
                (ctx->bits[SCL_BITS_IDX(ctx, l, row - offset, col)] ? -1.0 : 1.0) *
                    ctx->llrs[SCL_LLR_IDX(ctx, l, row - offset, col + 1)] +
                ctx->llrs[SCL_LLR_IDX(ctx, l, row, col + 1)];
        } else {
            if (!ctx->llrs_updated[(size_t)row * (size_t)ctx->n_cols + (size_t)(col + 1)]) {
                update_llr_rec(ctx, row, col + 1);
            }
            if (!ctx->llrs_updated[(size_t)(row + offset) * (size_t)ctx->n_cols + (size_t)(col + 1)]) {
                update_llr_rec(ctx, row + offset, col + 1);
            }
            ctx->llrs[SCL_LLR_IDX(ctx, l, row, col)] = minstar_scalar(
                ctx->llrs[SCL_LLR_IDX(ctx, l, row, col + 1)],
                ctx->llrs[SCL_LLR_IDX(ctx, l, row + offset, col + 1)],
                ctx->min_sum);
        }
    }

    ctx->llrs_updated[(size_t)row * (size_t)ctx->n_cols + (size_t)col] = 1;
}

static void scl_duplicate_paths_main(scl_ctx_t *ctx, int old_l) {
    int l;
    for (l = 0; l < old_l; l++) {
        copy_path_block_u8(
            ctx->bits + (size_t)(old_l + l) * (size_t)ctx->path_stride,
            ctx->bits + (size_t)l * (size_t)ctx->path_stride,
            ctx->path_stride);
        copy_path_block_f64(
            ctx->llrs + (size_t)(old_l + l) * (size_t)ctx->path_stride,
            ctx->llrs + (size_t)l * (size_t)ctx->path_stride,
            ctx->path_stride);
    }
}

static polar_status_t scl_split_info_bit(scl_ctx_t *ctx, int row, int *old_l_out) {
    int old_l = ctx->L_prime;
    int l;

    if (2 * old_l > ctx->L_capacity) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    scl_duplicate_paths_main(ctx, old_l);

    for (l = 0; l < old_l; l++) {
        double pm_old = ctx->PM[l];
        double llr = ctx->llrs[SCL_LLR_IDX(ctx, l, row, 0)];
        ctx->PM[l] = phi_scalar(pm_old, llr, 0, ctx->min_sum);
        ctx->PM[old_l + l] = phi_scalar(pm_old, llr, 1, ctx->min_sum);

        ctx->bits[SCL_BITS_IDX(ctx, l, row, 0)] = 0;
        ctx->bits[SCL_BITS_IDX(ctx, old_l + l, row, 0)] = 1;
    }

    ctx->bits_updated[(size_t)row * (size_t)ctx->n_cols + 0] = 1;
    ctx->L_prime = 2 * old_l;

    if (old_l_out != NULL) {
        *old_l_out = old_l;
    }
    return POLAR_STATUS_OK;
}

static void scl_force_bit(scl_ctx_t *ctx, int row, int bit, int mark_updated) {
    int l;
    for (l = 0; l < ctx->L_prime; l++) {
        double llr = ctx->llrs[SCL_LLR_IDX(ctx, l, row, 0)];
        ctx->PM[l] = phi_scalar(ctx->PM[l], llr, bit, ctx->min_sum);
        if (bit != 0) {
            ctx->bits[SCL_BITS_IDX(ctx, l, row, 0)] = 1;
        }
    }
    if (mark_updated) {
        ctx->bits_updated[(size_t)row * (size_t)ctx->n_cols + 0] = 1;
    }
}

static void prune_u8_paths(
    uint8_t *arr,
    uint8_t *tmp,
    int elem_per_path,
    const int *keep_indices,
    int keep_count
) {
    int j;
    for (j = 0; j < keep_count; j++) {
        int src = keep_indices[j];
        copy_path_block_u8(tmp + (size_t)j * (size_t)elem_per_path, arr + (size_t)src * (size_t)elem_per_path, elem_per_path);
    }
    memcpy(arr, tmp, (size_t)keep_count * (size_t)elem_per_path);
}

static void prune_f64_paths(
    double *arr,
    double *tmp,
    int elem_per_path,
    const int *keep_indices,
    int keep_count
) {
    int j;
    for (j = 0; j < keep_count; j++) {
        int src = keep_indices[j];
        copy_path_block_f64(tmp + (size_t)j * (size_t)elem_per_path, arr + (size_t)src * (size_t)elem_per_path, elem_per_path);
    }
    memcpy(arr, tmp, (size_t)keep_count * (size_t)elem_per_path * sizeof(double));
}

static void prune_u8_scalars(uint8_t *arr, uint8_t *tmp, const int *keep_indices, int keep_count) {
    int j;
    for (j = 0; j < keep_count; j++) {
        tmp[j] = arr[keep_indices[j]];
    }
    memcpy(arr, tmp, (size_t)keep_count);
}

static polar_status_t scl_prune_main(scl_ctx_t *ctx, int keep_count, int **indices_out) {
    int *indices = NULL;
    polar_status_t st;

    st = sorted_path_indices(ctx->PM, ctx->L_prime, &indices);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    prune_u8_paths(ctx->bits, ctx->bits_alt, ctx->path_stride, indices, keep_count);
    prune_f64_paths(ctx->llrs, ctx->llrs_alt, ctx->path_stride, indices, keep_count);
    prune_f64_paths(ctx->PM, ctx->PM_alt, 1, indices, keep_count);
    ctx->L_prime = keep_count;

    if (indices_out != NULL) {
        *indices_out = indices;
    } else {
        free(indices);
    }

    return POLAR_STATUS_OK;
}

static polar_status_t decode_plain_polar(
    const double *e_tilde,
    int E,
    const uint8_t *info_bit_pattern,
    int N,
    const int *rate_matching_pattern,
    int mode,
    int L,
    int min_sum,
    uint8_t **a_hat,
    int *a_len
) {
    double *d_tilde = NULL;
    scl_ctx_t ctx;
    int i;
    int *sorted = NULL;
    polar_status_t st;

    *a_hat = NULL;
    *a_len = 0;

    st = rate_dematch(e_tilde, E, N, rate_matching_pattern, mode, &d_tilde);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    st = scl_ctx_init(&ctx, N, L, min_sum, info_bit_pattern, d_tilde);
    free(d_tilde);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    for (i = 0; i < N; i++) {
        update_llr_rec(&ctx, i, 0);
        if (!info_bit_pattern[i]) {
            scl_force_bit(&ctx, i, 0, 0);
        } else {
            st = scl_split_info_bit(&ctx, i, NULL);
            if (st != POLAR_STATUS_OK) {
                scl_ctx_free(&ctx);
                return st;
            }
            if (ctx.L_prime > L) {
                st = scl_prune_main(&ctx, L, NULL);
                if (st != POLAR_STATUS_OK) {
                    scl_ctx_free(&ctx);
                    return st;
                }
            }
        }
    }

    st = sorted_path_indices(ctx.PM, ctx.L_prime, &sorted);
    if (st != POLAR_STATUS_OK) {
        scl_ctx_free(&ctx);
        return st;
    }

    st = extract_info_bits_from_scl(
        &ctx,
        sorted[0],
        info_bit_pattern,
        a_hat,
        a_len);

    free(sorted);
    scl_ctx_free(&ctx);
    return st;
}

static int crc_pass_with_matrix(const uint8_t *b_hat, int K, const u8_matrix_t *G_crc) {
    int p;
    for (p = 0; p < G_crc->cols; p++) {
        int k;
        uint8_t v = 0;
        for (k = 0; k < K; k++) {
            v ^= (uint8_t)(b_hat[k] & U8_AT((*G_crc), k, p));
        }
        if ((v & 1u) != 0u) {
            return 0;
        }
    }
    return 1;
}

static polar_status_t extract_info_not_pc(
    const scl_ctx_t *ctx,
    int path,
    const uint8_t *info_bit_pattern,
    const uint8_t *pc_bit_pattern,
    uint8_t **b_hat,
    int *k_out
) {
    int i;
    int K = 0;
    uint8_t *b;
    int pos = 0;
    polar_status_t st;

    for (i = 0; i < ctx->N; i++) {
        if (info_bit_pattern[i] && !pc_bit_pattern[i]) {
            K++;
        }
    }

    st = alloc_u8((size_t)K, &b);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    for (i = 0; i < ctx->N; i++) {
        if (info_bit_pattern[i] && !pc_bit_pattern[i]) {
            b[pos++] = ctx->bits[SCL_BITS_IDX(ctx, path, i, 0)];
        }
    }

    *b_hat = b;
    *k_out = K;
    return POLAR_STATUS_OK;
}

static void deinterleave_crc_bits(const uint8_t *c_hat, int K, const int *crc_interleaver, uint8_t *b_hat) {
    int i;
    for (i = 0; i < K; i++) {
        b_hat[crc_interleaver[i]] = c_hat[i];
    }
}

/* Decoder implementations are below. They follow the MATLAB structure closely. */

static polar_status_t ca_polar_decoder_impl(
    const double *e_tilde,
    int E,
    const uint8_t *crc_poly,
    int crc_poly_len,
    const uint8_t *info_bit_pattern,
    int N,
    const int *rate_matching_pattern,
    int mode,
    int L,
    int min_sum,
    int P2,
    uint8_t **a_hat,
    int *a_len
) {
    int K = 0;
    int P = crc_poly_len - 1;
    double *d_tilde = NULL;
    scl_ctx_t ctx;
    polar_status_t st;
    int i;
    int *sorted = NULL;
    u8_matrix_t G_crc = {0};
    int max_candidates;

    *a_hat = NULL;
    *a_len = 0;

    for (i = 0; i < N; i++) {
        if (info_bit_pattern[i]) {
            K++;
        }
    }

    st = rate_dematch(e_tilde, E, N, rate_matching_pattern, mode, &d_tilde);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    st = scl_ctx_init(&ctx, N, L, min_sum, info_bit_pattern, d_tilde);
    free(d_tilde);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    for (i = 0; i < N; i++) {
        update_llr_rec(&ctx, i, 0);
        if (!info_bit_pattern[i]) {
            scl_force_bit(&ctx, i, 0, 0);
        } else {
            st = scl_split_info_bit(&ctx, i, NULL);
            if (st != POLAR_STATUS_OK) {
                scl_ctx_free(&ctx);
                return st;
            }
            if (ctx.L_prime > L) {
                st = scl_prune_main(&ctx, L, NULL);
                if (st != POLAR_STATUS_OK) {
                    scl_ctx_free(&ctx);
                    return st;
                }
            }
        }
    }

    st = get_crc_generator_matrix(K, crc_poly, crc_poly_len, &G_crc);
    if (st != POLAR_STATUS_OK) {
        scl_ctx_free(&ctx);
        return st;
    }

    st = sorted_path_indices(ctx.PM, ctx.L_prime, &sorted);
    if (st != POLAR_STATUS_OK) {
        free_u8_matrix(&G_crc);
        scl_ctx_free(&ctx);
        return st;
    }

    max_candidates = L;
    if ((1 << P2) < max_candidates) {
        max_candidates = (1 << P2);
    }
    if (max_candidates > ctx.L_prime) {
        max_candidates = ctx.L_prime;
    }

    for (i = 0; i < max_candidates; i++) {
        uint8_t *b_hat = NULL;
        int b_len = 0;
        uint8_t *out = NULL;

        st = extract_info_bits_from_scl(
            &ctx,
            sorted[i],
            info_bit_pattern,
            &b_hat,
            &b_len);
        if (st != POLAR_STATUS_OK) {
            continue;
        }

        if (crc_pass_with_matrix(b_hat, b_len, &G_crc)) {
            st = alloc_u8((size_t)(b_len - P), &out);
            if (st == POLAR_STATUS_OK) {
                memcpy(out, b_hat, (size_t)(b_len - P));
                *a_hat = out;
                *a_len = b_len - P;
                free(b_hat);
                break;
            }
        }
        free(b_hat);
    }

    free(sorted);
    free_u8_matrix(&G_crc);
    scl_ctx_free(&ctx);
    return POLAR_STATUS_OK;
}

static polar_status_t dca_family_prepare(
    int A,
    int P,
    const uint8_t *crc_poly,
    int crc_poly_len,
    const int *crc_interleaver,
    int K,
    uint8_t **G_row_major,
    int **last_one_index,
    uint8_t ds1_mode,
    const uint8_t *crc_scrambling_pattern,
    int scrambling_len,
    uint8_t **extended_scrambling
) {
    u8_matrix_t Gp = {0};
    uint8_t *Gp2 = NULL;
    uint8_t *Gp3 = NULL;
    int *last = NULL;
    int p;
    int i;
    polar_status_t st;

    st = get_crc_generator_matrix(ds1_mode ? (A + P) : A, crc_poly, crc_poly_len, &Gp);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    st = alloc_u8((size_t)K * (size_t)P, &Gp2);
    if (st != POLAR_STATUS_OK) {
        free_u8_matrix(&Gp);
        return st;
    }

    if (!ds1_mode) {
        for (i = 0; i < A; i++) {
            memcpy(Gp2 + (size_t)i * (size_t)P, Gp.data + (size_t)i * (size_t)P, (size_t)P);
        }
        for (i = 0; i < P; i++) {
            int j;
            for (j = 0; j < P; j++) {
                Gp2[(size_t)(A + i) * (size_t)P + (size_t)j] = (uint8_t)(i == j);
            }
        }
    } else {
        for (i = 0; i < A; i++) {
            memcpy(Gp2 + (size_t)i * (size_t)P, Gp.data + (size_t)(P + i) * (size_t)P, (size_t)P);
        }
        for (i = 0; i < P; i++) {
            int j;
            for (j = 0; j < P; j++) {
                Gp2[(size_t)(A + i) * (size_t)P + (size_t)j] = (uint8_t)(i == j);
            }
        }
    }

    st = alloc_u8((size_t)K * (size_t)P, &Gp3);
    if (st != POLAR_STATUS_OK) {
        free_u8_matrix(&Gp);
        free(Gp2);
        return st;
    }

    for (i = 0; i < K; i++) {
        int idx = crc_interleaver[i];
        if (idx < 0 || idx >= K) {
            free_u8_matrix(&Gp);
            free(Gp2);
            free(Gp3);
            return POLAR_STATUS_INVALID_ARGUMENT;
        }
        memcpy(Gp3 + (size_t)i * (size_t)P, Gp2 + (size_t)idx * (size_t)P, (size_t)P);
    }

    st = alloc_i32((size_t)P, &last);
    if (st != POLAR_STATUS_OK) {
        free_u8_matrix(&Gp);
        free(Gp2);
        free(Gp3);
        return st;
    }

    for (p = 0; p < P; p++) {
        int found = 0;
        for (i = K - 1; i >= 0; i--) {
            if (Gp3[(size_t)i * (size_t)P + (size_t)p]) {
                last[p] = i;
                found = 1;
                break;
            }
        }
        if (!found) {
            last[p] = -1;
        }
    }

    *G_row_major = Gp3;
    *last_one_index = last;

    if (ds1_mode) {
        uint8_t *ext = NULL;
        st = calloc_u8((size_t)P, &ext);
        if (st != POLAR_STATUS_OK) {
            free_u8_matrix(&Gp);
            free(Gp2);
            free(Gp3);
            free(last);
            return st;
        }
        for (i = 0; i < scrambling_len; i++) {
            ext[P - scrambling_len + i] = crc_scrambling_pattern[i];
        }
        *extended_scrambling = ext;
    }

    free_u8_matrix(&Gp);
    free(Gp2);
    return POLAR_STATUS_OK;
}

static polar_status_t dca_like_decoder_impl(
    const double *e_tilde,
    int E,
    const uint8_t *crc_poly,
    int crc_poly_len,
    const int *crc_interleaver,
    int K,
    const uint8_t *info_bit_pattern,
    int N,
    const int *rate_matching_pattern,
    int mode,
    int L,
    int min_sum,
    int P2,
    int allow_known_bits,
    const int8_t *a_tilde,
    int ds1_mode,
    const uint8_t *crc_scrambling_pattern,
    int scrambling_len,
    uint8_t **a_hat,
    int *a_len
) {
    int P = crc_poly_len - 1;
    int A = K - P;
    double *d_tilde = NULL;
    scl_ctx_t ctx;
    polar_status_t st;
    uint8_t *G_rows = NULL;
    int *last_one_index = NULL;
    uint8_t *extended_scrambling = NULL;
    uint8_t *crc_checksums = NULL;
    uint8_t *crc_checksums_tmp = NULL;
    uint8_t *crc_okay = NULL;
    uint8_t *crc_okay_tmp = NULL;
    int i;
    int i2 = 0;
    int *sorted = NULL;
    int max_candidates;

    *a_hat = NULL;
    *a_len = 0;

    st = rate_dematch(e_tilde, E, N, rate_matching_pattern, mode, &d_tilde);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    st = scl_ctx_init(&ctx, N, L, min_sum, info_bit_pattern, d_tilde);
    free(d_tilde);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    st = dca_family_prepare(
        A,
        P,
        crc_poly,
        crc_poly_len,
        crc_interleaver,
        K,
        &G_rows,
        &last_one_index,
        (uint8_t)ds1_mode,
        crc_scrambling_pattern,
        scrambling_len,
        &extended_scrambling);
    if (st != POLAR_STATUS_OK) {
        scl_ctx_free(&ctx);
        return st;
    }

    st = calloc_u8((size_t)ctx.L_capacity * (size_t)P, &crc_checksums);
    if (st != POLAR_STATUS_OK) {
        scl_ctx_free(&ctx);
        free(G_rows);
        free(last_one_index);
        free(extended_scrambling);
        return st;
    }
    st = calloc_u8((size_t)ctx.L_capacity * (size_t)P, &crc_checksums_tmp);
    if (st != POLAR_STATUS_OK) {
        scl_ctx_free(&ctx);
        free(G_rows);
        free(last_one_index);
        free(extended_scrambling);
        free(crc_checksums);
        return st;
    }
    st = alloc_u8((size_t)ctx.L_capacity, &crc_okay);
    if (st != POLAR_STATUS_OK) {
        scl_ctx_free(&ctx);
        free(G_rows);
        free(last_one_index);
        free(extended_scrambling);
        free(crc_checksums);
        free(crc_checksums_tmp);
        return st;
    }
    st = alloc_u8((size_t)ctx.L_capacity, &crc_okay_tmp);
    if (st != POLAR_STATUS_OK) {
        scl_ctx_free(&ctx);
        free(G_rows);
        free(last_one_index);
        free(extended_scrambling);
        free(crc_checksums);
        free(crc_checksums_tmp);
        free(crc_okay);
        return st;
    }

    for (i = 0; i < ctx.L_capacity; i++) {
        crc_okay[i] = 1;
    }

    if (ds1_mode) {
        u8_matrix_t Ginit = {0};
        st = get_crc_generator_matrix(A + P, crc_poly, crc_poly_len, &Ginit);
        if (st != POLAR_STATUS_OK) {
            scl_ctx_free(&ctx);
            free(G_rows);
            free(last_one_index);
            free(extended_scrambling);
            free(crc_checksums);
            free(crc_checksums_tmp);
            free(crc_okay);
            free(crc_okay_tmp);
            return st;
        }

        for (i = 0; i < P; i++) {
            int r;
            uint8_t sum = 0;
            for (r = 0; r < P; r++) {
                sum ^= U8_AT(Ginit, r, i);
            }
            crc_checksums[i] = sum;
        }

        free_u8_matrix(&Ginit);
    }

    for (i = 0; i < N; i++) {
        update_llr_rec(&ctx, i, 0);

        if (!info_bit_pattern[i]) {
            scl_force_bit(&ctx, i, 0, 0);
            continue;
        }

        {
            int forced = 0;
            int known = -1;
            int bit_idx = crc_interleaver[i2];

            if (allow_known_bits && bit_idx < A && a_tilde != NULL) {
                if (a_tilde[bit_idx] == 0) {
                    known = 0;
                    forced = 1;
                } else if (a_tilde[bit_idx] == 1) {
                    known = 1;
                    forced = 1;
                }
            }

            if (forced) {
                int l;
                scl_force_bit(&ctx, i, known, 1);
                if (known == 1) {
                    for (l = 0; l < ctx.L_prime; l++) {
                        int p;
                        for (p = 0; p < P; p++) {
                            crc_checksums[(size_t)l * (size_t)P + (size_t)p] ^=
                                G_rows[(size_t)i2 * (size_t)P + (size_t)p];
                        }
                    }
                }
            } else {
                int old_l;
                st = scl_split_info_bit(&ctx, i, &old_l);
                if (st != POLAR_STATUS_OK) {
                    scl_ctx_free(&ctx);
                    free(G_rows);
                    free(last_one_index);
                    free(extended_scrambling);
                    free(crc_checksums);
                    free(crc_checksums_tmp);
                    free(crc_okay);
                    free(crc_okay_tmp);
                    return st;
                }

                for (int l = 0; l < old_l; l++) {
                    memcpy(
                        crc_checksums + (size_t)(old_l + l) * (size_t)P,
                        crc_checksums + (size_t)l * (size_t)P,
                        (size_t)P);
                    for (int p = 0; p < P; p++) {
                        crc_checksums[(size_t)(old_l + l) * (size_t)P + (size_t)p] ^=
                            G_rows[(size_t)i2 * (size_t)P + (size_t)p];
                    }
                    crc_okay[old_l + l] = crc_okay[l];
                }

                if (ctx.L_prime > L) {
                    int *indices = NULL;
                    st = scl_prune_main(&ctx, L, &indices);
                    if (st != POLAR_STATUS_OK) {
                        scl_ctx_free(&ctx);
                        free(G_rows);
                        free(last_one_index);
                        free(extended_scrambling);
                        free(crc_checksums);
                        free(crc_checksums_tmp);
                        free(crc_okay);
                        free(crc_okay_tmp);
                        return st;
                    }

                    prune_u8_paths(crc_checksums, crc_checksums_tmp, P, indices, L);
                    prune_u8_scalars(crc_okay, crc_okay_tmp, indices, L);
                    free(indices);
                }
            }
        }

        for (int p = 0; p < P; p++) {
            if (last_one_index[p] == i2) {
                int l;
                for (l = 0; l < ctx.L_prime; l++) {
                    uint8_t expected = ds1_mode ? extended_scrambling[p] : 0;
                    if (crc_checksums[(size_t)l * (size_t)P + (size_t)p] != expected) {
                        crc_okay[l] = 0;
                    }
                }
            }
        }

        {
            int any_ok = 0;
            for (int l = 0; l < ctx.L_prime; l++) {
                if (crc_okay[l]) {
                    any_ok = 1;
                    break;
                }
            }
            if (!any_ok) {
                free(G_rows);
                free(last_one_index);
                free(extended_scrambling);
                free(crc_checksums);
                free(crc_checksums_tmp);
                free(crc_okay);
                free(crc_okay_tmp);
                scl_ctx_free(&ctx);
                return POLAR_STATUS_OK;
            }
        }

        i2++;
    }

    st = sorted_path_indices(ctx.PM, ctx.L_prime, &sorted);
    if (st != POLAR_STATUS_OK) {
        free(G_rows);
        free(last_one_index);
        free(extended_scrambling);
        free(crc_checksums);
        free(crc_checksums_tmp);
        free(crc_okay);
        free(crc_okay_tmp);
        scl_ctx_free(&ctx);
        return st;
    }

    max_candidates = L;
    if ((1 << P2) < max_candidates) {
        max_candidates = (1 << P2);
    }
    if (max_candidates > ctx.L_prime) {
        max_candidates = ctx.L_prime;
    }

    for (i = 0; i < max_candidates; i++) {
        int path = sorted[i];
        uint8_t *c_hat = NULL;
        int c_len = 0;
        uint8_t *b_hat = NULL;
        uint8_t *out = NULL;

        if (!crc_okay[path]) {
            continue;
        }

        st = extract_info_bits_from_scl(
            &ctx,
            path,
            info_bit_pattern,
            &c_hat,
            &c_len);
        if (st != POLAR_STATUS_OK) {
            continue;
        }

        st = alloc_u8((size_t)K, &b_hat);
        if (st != POLAR_STATUS_OK) {
            free(c_hat);
            continue;
        }
        memset(b_hat, 0, (size_t)K);
        deinterleave_crc_bits(c_hat, K, crc_interleaver, b_hat);

        st = alloc_u8((size_t)A, &out);
        if (st == POLAR_STATUS_OK) {
            memcpy(out, b_hat, (size_t)A);
            *a_hat = out;
            *a_len = A;
            free(c_hat);
            free(b_hat);
            break;
        }

        free(c_hat);
        free(b_hat);
    }

    free(sorted);
    free(G_rows);
    free(last_one_index);
    free(extended_scrambling);
    free(crc_checksums);
    free(crc_checksums_tmp);
    free(crc_okay);
    free(crc_okay_tmp);
    scl_ctx_free(&ctx);
    return POLAR_STATUS_OK;
}

static polar_status_t pc_family_decoder_impl(
    const double *e_tilde,
    int E,
    const uint8_t *crc_poly,
    int crc_poly_len,
    int with_crc,
    const uint8_t *info_bit_pattern,
    const uint8_t *pc_bit_pattern,
    int N,
    int pc_buf_len,
    const int *rate_matching_pattern,
    int mode,
    int L,
    int min_sum,
    int P2,
    uint8_t **a_hat,
    int *a_len
) {
    double *d_tilde = NULL;
    scl_ctx_t ctx;
    uint8_t *y = NULL;
    uint8_t *y_tmp = NULL;
    polar_status_t st;
    int i;
    int *sorted = NULL;
    u8_matrix_t G_crc = {0};
    int max_candidates;

    *a_hat = NULL;
    *a_len = 0;

    st = rate_dematch(e_tilde, E, N, rate_matching_pattern, mode, &d_tilde);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    st = scl_ctx_init(&ctx, N, L, min_sum, info_bit_pattern, d_tilde);
    free(d_tilde);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    st = calloc_u8((size_t)ctx.L_capacity * (size_t)pc_buf_len, &y);
    if (st != POLAR_STATUS_OK) {
        scl_ctx_free(&ctx);
        return st;
    }
    st = calloc_u8((size_t)ctx.L_capacity * (size_t)pc_buf_len, &y_tmp);
    if (st != POLAR_STATUS_OK) {
        scl_ctx_free(&ctx);
        free(y);
        return st;
    }

    for (i = 0; i < N; i++) {
        int l;
        for (l = 0; l < ctx.L_prime; l++) {
            uint8_t first = y[(size_t)l * (size_t)pc_buf_len + 0];
            memmove(
                y + (size_t)l * (size_t)pc_buf_len,
                y + (size_t)l * (size_t)pc_buf_len + 1,
                (size_t)(pc_buf_len - 1));
            y[(size_t)l * (size_t)pc_buf_len + (size_t)(pc_buf_len - 1)] = first;
        }

        update_llr_rec(&ctx, i, 0);

        if (!info_bit_pattern[i]) {
            scl_force_bit(&ctx, i, 0, 0);
        } else if (pc_bit_pattern[i]) {
            for (l = 0; l < ctx.L_prime; l++) {
                int bit = y[(size_t)l * (size_t)pc_buf_len + 0] ? 1 : 0;
                double llr = ctx.llrs[SCL_LLR_IDX(&ctx, l, i, 0)];
                ctx.PM[l] = phi_scalar(ctx.PM[l], llr, bit, ctx.min_sum);
                ctx.bits[SCL_BITS_IDX(&ctx, l, i, 0)] = (uint8_t)bit;
            }
            ctx.bits_updated[(size_t)i * (size_t)ctx.n_cols + 0] = 1;
        } else {
            int old_l;
            st = scl_split_info_bit(&ctx, i, &old_l);
            if (st != POLAR_STATUS_OK) {
                scl_ctx_free(&ctx);
                free(y);
                free(y_tmp);
                return st;
            }

            for (l = 0; l < old_l; l++) {
                memcpy(
                    y + (size_t)(old_l + l) * (size_t)pc_buf_len,
                    y + (size_t)l * (size_t)pc_buf_len,
                    (size_t)pc_buf_len);
            }

            for (l = 0; l < ctx.L_prime; l++) {
                y[(size_t)l * (size_t)pc_buf_len + 0] ^=
                    ctx.bits[SCL_BITS_IDX(&ctx, l, i, 0)];
            }

            if (ctx.L_prime > L) {
                int *indices = NULL;
                st = scl_prune_main(&ctx, L, &indices);
                if (st != POLAR_STATUS_OK) {
                    scl_ctx_free(&ctx);
                    free(y);
                    free(y_tmp);
                    return st;
                }
                prune_u8_paths(y, y_tmp, pc_buf_len, indices, L);
                free(indices);
            }
        }
    }

    if (with_crc) {
        int K = 0;
        for (i = 0; i < N; i++) {
            if (info_bit_pattern[i] && !pc_bit_pattern[i]) {
                K++;
            }
        }
        st = get_crc_generator_matrix(K, crc_poly, crc_poly_len, &G_crc);
        if (st != POLAR_STATUS_OK) {
            scl_ctx_free(&ctx);
            free(y);
            free(y_tmp);
            return st;
        }
    }

    st = sorted_path_indices(ctx.PM, ctx.L_prime, &sorted);
    if (st != POLAR_STATUS_OK) {
        scl_ctx_free(&ctx);
        free(y);
        free(y_tmp);
        free_u8_matrix(&G_crc);
        return st;
    }

    max_candidates = with_crc ? L : 1;
    if (with_crc && (1 << P2) < max_candidates) {
        max_candidates = (1 << P2);
    }
    if (max_candidates > ctx.L_prime) {
        max_candidates = ctx.L_prime;
    }

    for (i = 0; i < max_candidates; i++) {
        uint8_t *b_hat = NULL;
        int b_len = 0;

        st = extract_info_not_pc(
            &ctx,
            sorted[i],
            info_bit_pattern,
            pc_bit_pattern,
            &b_hat,
            &b_len);
        if (st != POLAR_STATUS_OK) {
            continue;
        }

        if (!with_crc) {
            *a_hat = b_hat;
            *a_len = b_len;
            break;
        }

        if (crc_pass_with_matrix(b_hat, b_len, &G_crc)) {
            uint8_t *out = NULL;
            int P = crc_poly_len - 1;
            st = alloc_u8((size_t)(b_len - P), &out);
            if (st == POLAR_STATUS_OK) {
                memcpy(out, b_hat, (size_t)(b_len - P));
                *a_hat = out;
                *a_len = b_len - P;
                free(b_hat);
                break;
            }
        }

        free(b_hat);
    }

    free(sorted);
    free_u8_matrix(&G_crc);
    scl_ctx_free(&ctx);
    free(y);
    free(y_tmp);
    return POLAR_STATUS_OK;
}

void polar_free(void *ptr) {
    free(ptr);
}

polar_status_t polar_get_3gpp_n(int K, int E, int n_max, int *N_out) {
    return get_3gpp_n_impl(K, E, n_max, N_out);
}

polar_status_t polar_get_3gpp_rate_matching_pattern(
    int K,
    int N,
    int E,
    polar_i32_vec_t *pattern,
    int *mode_out
) {
    int *p = NULL;
    polar_status_t st;

    if (pattern == NULL || mode_out == NULL) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    st = get_3gpp_rate_matching_pattern_impl(K, N, E, &p, mode_out);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    pattern->data = p;
    pattern->len = E;
    return POLAR_STATUS_OK;
}

polar_status_t polar_get_3gpp_channel_interleaver_pattern(int E, polar_i32_vec_t *pattern) {
    int *p = NULL;
    polar_status_t st;

    if (pattern == NULL) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    st = get_3gpp_channel_interleaver_pattern_impl(E, &p);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    pattern->data = p;
    pattern->len = E;
    return POLAR_STATUS_OK;
}

polar_status_t polar_pbch_encode(const uint8_t *a, int A, int E, polar_u8_vec_t *f) {
    int K;
    int N;
    int *crc_interleaver = NULL;
    int *rate_pattern = NULL;
    int mode;
    int *Q_N = NULL;
    uint8_t *info_pattern = NULL;
    uint8_t *out = NULL;
    polar_status_t st;

    if (a == NULL || f == NULL) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    if (A != 32) {
        return POLAR_STATUS_UNSUPPORTED_BLOCK_LENGTH;
    }
    if (E <= 0) {
        E = 864;
    }
    if (E != 864) {
        return POLAR_STATUS_UNSUPPORTED_BLOCK_LENGTH;
    }

    K = A + ((int)sizeof(CRC24_POLY) - 1);

    st = get_3gpp_n_impl(K, E, 9, &N);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    st = get_3gpp_crc_interleaver_pattern_impl(K, &crc_interleaver);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    st = get_3gpp_rate_matching_pattern_impl(K, N, E, &rate_pattern, &mode);
    if (st != POLAR_STATUS_OK) {
        free(crc_interleaver);
        return st;
    }
    st = get_3gpp_sequence_pattern_impl(N, &Q_N);
    if (st != POLAR_STATUS_OK) {
        free(crc_interleaver);
        free(rate_pattern);
        return st;
    }
    st = get_3gpp_info_bit_pattern_impl(K, Q_N, N, rate_pattern, E, mode, &info_pattern);
    if (st != POLAR_STATUS_OK) {
        free(crc_interleaver);
        free(rate_pattern);
        free(Q_N);
        return st;
    }

    st = dca_polar_encoder_impl(
        a,
        A,
        CRC24_POLY,
        (int)sizeof(CRC24_POLY),
        crc_interleaver,
        K,
        info_pattern,
        N,
        rate_pattern,
        E,
        &out);

    free(crc_interleaver);
    free(rate_pattern);
    free(Q_N);
    free(info_pattern);

    if (st != POLAR_STATUS_OK) {
        return st;
    }

    f->data = out;
    f->len = E;
    return POLAR_STATUS_OK;
}

polar_status_t polar_pbch_decode(
    const double *f_tilde,
    int E,
    int A,
    int L,
    int min_sum,
    const int8_t *a_tilde,
    polar_u8_vec_t *a_hat
) {
    int K;
    int N;
    int *crc_interleaver = NULL;
    int *rate_pattern = NULL;
    int mode;
    int *Q_N = NULL;
    uint8_t *info_pattern = NULL;
    int8_t *known = NULL;
    uint8_t *out = NULL;
    int out_len = 0;
    polar_status_t st;

    if (f_tilde == NULL || a_hat == NULL || L <= 0) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    if (A != 32 || E != 864) {
        return POLAR_STATUS_UNSUPPORTED_BLOCK_LENGTH;
    }

    K = A + ((int)sizeof(CRC24_POLY) - 1);

    st = get_3gpp_n_impl(K, E, 9, &N);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    st = get_3gpp_crc_interleaver_pattern_impl(K, &crc_interleaver);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    st = get_3gpp_rate_matching_pattern_impl(K, N, E, &rate_pattern, &mode);
    if (st != POLAR_STATUS_OK) {
        free(crc_interleaver);
        return st;
    }
    st = get_3gpp_sequence_pattern_impl(N, &Q_N);
    if (st != POLAR_STATUS_OK) {
        free(crc_interleaver);
        free(rate_pattern);
        return st;
    }
    st = get_3gpp_info_bit_pattern_impl(K, Q_N, N, rate_pattern, E, mode, &info_pattern);
    if (st != POLAR_STATUS_OK) {
        free(crc_interleaver);
        free(rate_pattern);
        free(Q_N);
        return st;
    }

    if (a_tilde == NULL) {
        st = alloc_u8((size_t)A, (uint8_t **)&known);
        if (st != POLAR_STATUS_OK) {
            free(crc_interleaver);
            free(rate_pattern);
            free(Q_N);
            free(info_pattern);
            return st;
        }
        memset(known, (uint8_t)0xFF, (size_t)A);
    } else {
        known = (int8_t *)a_tilde;
    }

    st = dca_like_decoder_impl(
        f_tilde,
        E,
        CRC24_POLY,
        (int)sizeof(CRC24_POLY),
        crc_interleaver,
        K,
        info_pattern,
        N,
        rate_pattern,
        mode,
        L,
        min_sum,
        3,
        1,
        known,
        0,
        NULL,
        0,
        &out,
        &out_len);

    if (a_tilde == NULL) {
        free(known);
    }

    free(crc_interleaver);
    free(rate_pattern);
    free(Q_N);
    free(info_pattern);

    if (st != POLAR_STATUS_OK) {
        return st;
    }

    a_hat->data = out;
    a_hat->len = out_len;
    return POLAR_STATUS_OK;
}

polar_status_t polar_pdcch_encode(
    const uint8_t *a,
    int A,
    int E,
    const uint8_t *rnti,
    polar_u8_vec_t *f
) {
    int K;
    int N;
    int *crc_interleaver = NULL;
    int *rate_pattern = NULL;
    int mode;
    int *Q_N = NULL;
    uint8_t *info_pattern = NULL;
    uint8_t rnti_default[16];
    uint8_t *a_pad = NULL;
    uint8_t *out = NULL;
    polar_status_t st;
    int i;

    if (a == NULL || f == NULL) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }
    if (A < 1 || A > 140 || E > 8192 || E <= 0) {
        return POLAR_STATUS_UNSUPPORTED_BLOCK_LENGTH;
    }

    if (rnti == NULL) {
        for (i = 0; i < 16; i++) {
            rnti_default[i] = 1;
        }
        rnti = rnti_default;
    }

    if (A < 12) {
        st = alloc_u8(12, &a_pad);
        if (st != POLAR_STATUS_OK) {
            return st;
        }
        memcpy(a_pad, a, (size_t)A);
        memset(a_pad + A, 0, (size_t)(12 - A));
        K = 12 + ((int)sizeof(CRC24_POLY) - 1);
    } else {
        K = A + ((int)sizeof(CRC24_POLY) - 1);
    }

    st = get_3gpp_n_impl(K, E, 9, &N);
    if (st != POLAR_STATUS_OK) {
        free(a_pad);
        return st;
    }
    st = get_3gpp_crc_interleaver_pattern_impl(K, &crc_interleaver);
    if (st != POLAR_STATUS_OK) {
        free(a_pad);
        return st;
    }
    st = get_3gpp_rate_matching_pattern_impl(K, N, E, &rate_pattern, &mode);
    if (st != POLAR_STATUS_OK) {
        free(a_pad);
        free(crc_interleaver);
        return st;
    }
    st = get_3gpp_sequence_pattern_impl(N, &Q_N);
    if (st != POLAR_STATUS_OK) {
        free(a_pad);
        free(crc_interleaver);
        free(rate_pattern);
        return st;
    }
    st = get_3gpp_info_bit_pattern_impl(K, Q_N, N, rate_pattern, E, mode, &info_pattern);
    if (st != POLAR_STATUS_OK) {
        free(a_pad);
        free(crc_interleaver);
        free(rate_pattern);
        free(Q_N);
        return st;
    }

    st = ds1ca_polar_encoder_impl(
        (A < 12) ? a_pad : a,
        (A < 12) ? 12 : A,
        CRC24_POLY,
        (int)sizeof(CRC24_POLY),
        rnti,
        16,
        crc_interleaver,
        K,
        info_pattern,
        N,
        rate_pattern,
        E,
        &out);

    free(a_pad);
    free(crc_interleaver);
    free(rate_pattern);
    free(Q_N);
    free(info_pattern);

    if (st != POLAR_STATUS_OK) {
        return st;
    }

    f->data = out;
    f->len = E;
    return POLAR_STATUS_OK;
}

polar_status_t polar_pdcch_decode(
    const double *f_tilde,
    int E,
    int A,
    int L,
    int min_sum,
    const uint8_t *rnti,
    polar_u8_vec_t *a_hat
) {
    int K;
    int N;
    int *crc_interleaver = NULL;
    int *rate_pattern = NULL;
    int mode;
    int *Q_N = NULL;
    uint8_t *info_pattern = NULL;
    uint8_t rnti_default[16];
    int8_t *known = NULL;
    uint8_t *out = NULL;
    int out_len = 0;
    polar_status_t st;
    int i;

    if (f_tilde == NULL || a_hat == NULL || L <= 0) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }
    if (A < 1 || A > 140 || E > 8192 || E <= 0) {
        return POLAR_STATUS_UNSUPPORTED_BLOCK_LENGTH;
    }

    if (rnti == NULL) {
        for (i = 0; i < 16; i++) {
            rnti_default[i] = 1;
        }
        rnti = rnti_default;
    }

    if (A < 12) {
        K = 12 + ((int)sizeof(CRC24_POLY) - 1);
    } else {
        K = A + ((int)sizeof(CRC24_POLY) - 1);
    }

    st = get_3gpp_n_impl(K, E, 9, &N);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    st = get_3gpp_crc_interleaver_pattern_impl(K, &crc_interleaver);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    st = get_3gpp_rate_matching_pattern_impl(K, N, E, &rate_pattern, &mode);
    if (st != POLAR_STATUS_OK) {
        free(crc_interleaver);
        return st;
    }
    st = get_3gpp_sequence_pattern_impl(N, &Q_N);
    if (st != POLAR_STATUS_OK) {
        free(crc_interleaver);
        free(rate_pattern);
        return st;
    }
    st = get_3gpp_info_bit_pattern_impl(K, Q_N, N, rate_pattern, E, mode, &info_pattern);
    if (st != POLAR_STATUS_OK) {
        free(crc_interleaver);
        free(rate_pattern);
        free(Q_N);
        return st;
    }

    if (A < 12) {
        st = alloc_u8(12, (uint8_t **)&known);
        if (st != POLAR_STATUS_OK) {
            free(crc_interleaver);
            free(rate_pattern);
            free(Q_N);
            free(info_pattern);
            return st;
        }
        for (i = 0; i < A; i++) {
            known[i] = -1;
        }
        for (i = A; i < 12; i++) {
            known[i] = 0;
        }

        st = dca_like_decoder_impl(
            f_tilde,
            E,
            CRC24_POLY,
            (int)sizeof(CRC24_POLY),
            crc_interleaver,
            K,
            info_pattern,
            N,
            rate_pattern,
            mode,
            L,
            min_sum,
            3,
            1,
            known,
            1,
            rnti,
            16,
            &out,
            &out_len);

        free(known);

        if (st == POLAR_STATUS_OK && out != NULL && out_len > A) {
            uint8_t *trim = NULL;
            st = alloc_u8((size_t)A, &trim);
            if (st == POLAR_STATUS_OK) {
                memcpy(trim, out, (size_t)A);
                free(out);
                out = trim;
                out_len = A;
            }
        }
    } else {
        st = dca_like_decoder_impl(
            f_tilde,
            E,
            CRC24_POLY,
            (int)sizeof(CRC24_POLY),
            crc_interleaver,
            K,
            info_pattern,
            N,
            rate_pattern,
            mode,
            L,
            min_sum,
            3,
            0,
            NULL,
            1,
            rnti,
            16,
            &out,
            &out_len);
    }

    free(crc_interleaver);
    free(rate_pattern);
    free(Q_N);
    free(info_pattern);

    if (st != POLAR_STATUS_OK) {
        return st;
    }

    a_hat->data = out;
    a_hat->len = out_len;
    return POLAR_STATUS_OK;
}

polar_status_t polar_pucch_encode(const uint8_t *a, int A, int G, polar_u8_vec_t *f) {
    const uint8_t *crc_poly;
    int crc_poly_len;
    int C;
    int P;
    int K;
    int E_r;
    int N;
    int *rate_pattern = NULL;
    int mode;
    int *Q_N = NULL;
    int *channel_pattern = NULL;
    uint8_t *info_pattern = NULL;
    uint8_t *info_pattern2 = NULL;
    uint8_t *pc_pattern = NULL;
    uint8_t *e = NULL;
    uint8_t *out = NULL;
    polar_status_t st;

    if (a == NULL || f == NULL || G <= 0) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }
    if (A < 12 || A > 1706) {
        return POLAR_STATUS_UNSUPPORTED_BLOCK_LENGTH;
    }

    if (A <= 19) {
        crc_poly = CRC6_POLY;
        crc_poly_len = (int)sizeof(CRC6_POLY);
        C = 1;
    } else {
        crc_poly = CRC11_POLY;
        crc_poly_len = (int)sizeof(CRC11_POLY);
        if ((A >= 360 && G >= 1088) || A >= 1013) {
            C = 2;
        } else {
            C = 1;
        }
    }

    P = crc_poly_len - 1;
    K = (A + C - 1) / C + P;
    E_r = G / C;

    if (E_r > 8192) {
        return POLAR_STATUS_UNSUPPORTED_BLOCK_LENGTH;
    }

    st = get_3gpp_n_impl(K, E_r, 10, &N);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    st = get_3gpp_rate_matching_pattern_impl(K, N, E_r, &rate_pattern, &mode);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    st = get_3gpp_sequence_pattern_impl(N, &Q_N);
    if (st != POLAR_STATUS_OK) {
        free(rate_pattern);
        return st;
    }
    st = get_3gpp_channel_interleaver_pattern_impl(E_r, &channel_pattern);
    if (st != POLAR_STATUS_OK) {
        free(rate_pattern);
        free(Q_N);
        return st;
    }

    if (A <= 19) {
        int n_PC = 3;
        int n_PC_wm = (G - K + 3 > 192) ? 1 : 0;

        st = get_3gpp_info_bit_pattern_impl(K + n_PC, Q_N, N, rate_pattern, E_r, mode, &info_pattern);
        if (st != POLAR_STATUS_OK) {
            free(rate_pattern);
            free(Q_N);
            free(channel_pattern);
            return st;
        }
        st = get_pc_bit_pattern_impl(info_pattern, Q_N, N, n_PC, n_PC_wm, &pc_pattern);
        if (st != POLAR_STATUS_OK) {
            free(rate_pattern);
            free(Q_N);
            free(channel_pattern);
            free(info_pattern);
            return st;
        }

        st = pcca_polar_encoder_impl(
            a,
            A,
            crc_poly,
            crc_poly_len,
            info_pattern,
            pc_pattern,
            N,
            5,
            rate_pattern,
            E_r,
            &e);
        if (st != POLAR_STATUS_OK) {
            free(rate_pattern);
            free(Q_N);
            free(channel_pattern);
            free(info_pattern);
            free(pc_pattern);
            return st;
        }

        st = alloc_u8((size_t)G, &out);
        if (st == POLAR_STATUS_OK) {
            for (int i = 0; i < E_r; i++) {
                out[i] = e[channel_pattern[i]];
            }
        }
    } else {
        st = get_3gpp_info_bit_pattern_impl(K, Q_N, N, rate_pattern, E_r, mode, &info_pattern);
        if (st != POLAR_STATUS_OK) {
            free(rate_pattern);
            free(Q_N);
            free(channel_pattern);
            return st;
        }

        if (C == 2) {
            int A1 = A / 2;
            int A2 = A - A1;

            st = alloc_u8((size_t)N, &info_pattern2);
            if (st != POLAR_STATUS_OK) {
                free(rate_pattern);
                free(Q_N);
                free(channel_pattern);
                free(info_pattern);
                return st;
            }
            memcpy(info_pattern2, info_pattern, (size_t)N);
            if (A % 2 != 0) {
                for (int i = 0; i < N; i++) {
                    if (info_pattern2[i]) {
                        info_pattern2[i] = 0;
                        break;
                    }
                }
            }

            st = ca_polar_encoder_impl(a, A1, crc_poly, crc_poly_len, info_pattern2, N, rate_pattern, E_r, &e);
            if (st != POLAR_STATUS_OK) {
                free(rate_pattern);
                free(Q_N);
                free(channel_pattern);
                free(info_pattern);
                free(info_pattern2);
                return st;
            }

            st = calloc_u8((size_t)G, &out);
            if (st != POLAR_STATUS_OK) {
                free(rate_pattern);
                free(Q_N);
                free(channel_pattern);
                free(info_pattern);
                free(info_pattern2);
                free(e);
                return st;
            }
            for (int i = 0; i < E_r; i++) {
                out[i] = e[channel_pattern[i]];
            }
            free(e);
            e = NULL;

            st = ca_polar_encoder_impl(a + A1, A2, crc_poly, crc_poly_len, info_pattern, N, rate_pattern, E_r, &e);
            if (st != POLAR_STATUS_OK) {
                free(rate_pattern);
                free(Q_N);
                free(channel_pattern);
                free(info_pattern);
                free(info_pattern2);
                free(out);
                return st;
            }
            for (int i = 0; i < E_r; i++) {
                out[E_r + i] = e[channel_pattern[i]];
            }
        } else {
            st = ca_polar_encoder_impl(a, A, crc_poly, crc_poly_len, info_pattern, N, rate_pattern, E_r, &e);
            if (st != POLAR_STATUS_OK) {
                free(rate_pattern);
                free(Q_N);
                free(channel_pattern);
                free(info_pattern);
                return st;
            }

            st = alloc_u8((size_t)G, &out);
            if (st == POLAR_STATUS_OK) {
                for (int i = 0; i < E_r; i++) {
                    out[i] = e[channel_pattern[i]];
                }
            }
        }
    }

    free(rate_pattern);
    free(Q_N);
    free(channel_pattern);
    free(info_pattern);
    free(info_pattern2);
    free(pc_pattern);
    free(e);

    if (st != POLAR_STATUS_OK) {
        return st;
    }

    f->data = out;
    f->len = G;
    return POLAR_STATUS_OK;
}

polar_status_t polar_pucch_decode(
    const double *f_tilde,
    int G,
    int A,
    int L,
    int min_sum,
    polar_u8_vec_t *a_hat
) {
    const uint8_t *crc_poly;
    int crc_poly_len;
    int C;
    int P;
    int K;
    int E_r;
    int N;
    int *rate_pattern = NULL;
    int mode;
    int *Q_N = NULL;
    int *channel_pattern = NULL;
    uint8_t *info_pattern = NULL;
    uint8_t *info_pattern2 = NULL;
    uint8_t *pc_pattern = NULL;
    double *e_tilde = NULL;
    uint8_t *out = NULL;
    int out_len = 0;
    polar_status_t st;

    if (f_tilde == NULL || a_hat == NULL || L <= 0 || G <= 0) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }
    if (A < 12 || A > 1706) {
        return POLAR_STATUS_UNSUPPORTED_BLOCK_LENGTH;
    }

    if (A <= 19) {
        crc_poly = CRC6_POLY;
        crc_poly_len = (int)sizeof(CRC6_POLY);
        C = 1;
    } else {
        crc_poly = CRC11_POLY;
        crc_poly_len = (int)sizeof(CRC11_POLY);
        if ((A >= 360 && G >= 1088) || A >= 1013) {
            C = 2;
        } else {
            C = 1;
        }
    }

    P = crc_poly_len - 1;
    K = (A + C - 1) / C + P;
    E_r = G / C;

    if (E_r > 8192) {
        return POLAR_STATUS_UNSUPPORTED_BLOCK_LENGTH;
    }

    st = get_3gpp_n_impl(K, E_r, 10, &N);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    st = get_3gpp_rate_matching_pattern_impl(K, N, E_r, &rate_pattern, &mode);
    if (st != POLAR_STATUS_OK) {
        return st;
    }
    st = get_3gpp_sequence_pattern_impl(N, &Q_N);
    if (st != POLAR_STATUS_OK) {
        free(rate_pattern);
        return st;
    }
    st = get_3gpp_channel_interleaver_pattern_impl(E_r, &channel_pattern);
    if (st != POLAR_STATUS_OK) {
        free(rate_pattern);
        free(Q_N);
        return st;
    }

    if (A <= 19) {
        int n_PC = 3;
        int n_PC_wm = (G - K + 3 > 192) ? 1 : 0;

        st = get_3gpp_info_bit_pattern_impl(K + n_PC, Q_N, N, rate_pattern, E_r, mode, &info_pattern);
        if (st != POLAR_STATUS_OK) {
            free(rate_pattern);
            free(Q_N);
            free(channel_pattern);
            return st;
        }
        st = get_pc_bit_pattern_impl(info_pattern, Q_N, N, n_PC, n_PC_wm, &pc_pattern);
        if (st != POLAR_STATUS_OK) {
            free(rate_pattern);
            free(Q_N);
            free(channel_pattern);
            free(info_pattern);
            return st;
        }

        st = alloc_f64((size_t)G, &e_tilde);
        if (st != POLAR_STATUS_OK) {
            free(rate_pattern);
            free(Q_N);
            free(channel_pattern);
            free(info_pattern);
            free(pc_pattern);
            return st;
        }
        for (int i = 0; i < E_r; i++) {
            e_tilde[channel_pattern[i]] = f_tilde[i];
        }

        st = pc_family_decoder_impl(
            e_tilde,
            E_r,
            crc_poly,
            crc_poly_len,
            1,
            info_pattern,
            pc_pattern,
            N,
            5,
            rate_pattern,
            mode,
            L,
            min_sum,
            3,
            &out,
            &out_len);
    } else {
        st = get_3gpp_info_bit_pattern_impl(K, Q_N, N, rate_pattern, E_r, mode, &info_pattern);
        if (st != POLAR_STATUS_OK) {
            free(rate_pattern);
            free(Q_N);
            free(channel_pattern);
            return st;
        }

        if (C == 2) {
            int A1 = A / 2;
            int A2 = A - A1;
            uint8_t *out1 = NULL;
            int out1_len = 0;
            uint8_t *out2 = NULL;
            int out2_len = 0;

            st = alloc_u8((size_t)N, &info_pattern2);
            if (st != POLAR_STATUS_OK) {
                free(rate_pattern);
                free(Q_N);
                free(channel_pattern);
                free(info_pattern);
                return st;
            }
            memcpy(info_pattern2, info_pattern, (size_t)N);
            if (A % 2 != 0) {
                for (int i = 0; i < N; i++) {
                    if (info_pattern2[i]) {
                        info_pattern2[i] = 0;
                        break;
                    }
                }
            }

            st = alloc_f64((size_t)E_r, &e_tilde);
            if (st != POLAR_STATUS_OK) {
                free(rate_pattern);
                free(Q_N);
                free(channel_pattern);
                free(info_pattern);
                free(info_pattern2);
                return st;
            }
            for (int i = 0; i < E_r; i++) {
                e_tilde[channel_pattern[i]] = f_tilde[i];
            }

            st = ca_polar_decoder_impl(
                e_tilde,
                E_r,
                crc_poly,
                crc_poly_len,
                info_pattern2,
                N,
                rate_pattern,
                mode,
                L,
                min_sum,
                3,
                &out1,
                &out1_len);
            if (st != POLAR_STATUS_OK) {
                free(rate_pattern);
                free(Q_N);
                free(channel_pattern);
                free(info_pattern);
                free(info_pattern2);
                free(e_tilde);
                return st;
            }

            if (out1_len == A1) {
                for (int i = 0; i < E_r; i++) {
                    e_tilde[channel_pattern[i]] = f_tilde[E_r + i];
                }
                st = ca_polar_decoder_impl(
                    e_tilde,
                    E_r,
                    crc_poly,
                    crc_poly_len,
                    info_pattern,
                    N,
                    rate_pattern,
                    mode,
                    L,
                    min_sum,
                    3,
                    &out2,
                    &out2_len);
            }

            if (st == POLAR_STATUS_OK && out1_len == A1 && out2_len == A2) {
                st = alloc_u8((size_t)A, &out);
                if (st == POLAR_STATUS_OK) {
                    memcpy(out, out1, (size_t)A1);
                    memcpy(out + A1, out2, (size_t)A2);
                    out_len = A;
                }
            }

            free(out1);
            free(out2);
        } else {
            st = alloc_f64((size_t)G, &e_tilde);
            if (st != POLAR_STATUS_OK) {
                free(rate_pattern);
                free(Q_N);
                free(channel_pattern);
                free(info_pattern);
                return st;
            }
            for (int i = 0; i < E_r; i++) {
                e_tilde[channel_pattern[i]] = f_tilde[i];
            }

            st = ca_polar_decoder_impl(
                e_tilde,
                E_r,
                crc_poly,
                crc_poly_len,
                info_pattern,
                N,
                rate_pattern,
                mode,
                L,
                min_sum,
                3,
                &out,
                &out_len);
        }
    }

    free(rate_pattern);
    free(Q_N);
    free(channel_pattern);
    free(info_pattern);
    free(info_pattern2);
    free(pc_pattern);
    free(e_tilde);

    if (st != POLAR_STATUS_OK) {
        return st;
    }

    a_hat->data = out;
    a_hat->len = out_len;
    return POLAR_STATUS_OK;
}

polar_status_t polar_custom1_encode(const uint8_t *a, int A, int E, polar_u8_vec_t *f) {
    int K;
    int N;
    int *rate_pattern = NULL;
    int mode;
    int *Q_N = NULL;
    uint8_t *info_pattern = NULL;
    int *channel_pattern = NULL;
    uint8_t *e = NULL;
    uint8_t *out = NULL;
    polar_status_t st;

    if (a == NULL || f == NULL || A <= 0 || E <= 0) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    K = A + ((int)sizeof(CRC11_POLY) - 1);

    st = get_3gpp_n_impl(K, E, 1024, &N);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    st = get_3gpp_rate_matching_pattern_impl(K, N, E, &rate_pattern, &mode);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    st = get_pw_sequence_pattern_impl(N, &Q_N);
    if (st != POLAR_STATUS_OK) {
        free(rate_pattern);
        return st;
    }

    st = get_3gpp_info_bit_pattern_impl(K, Q_N, N, rate_pattern, E, mode, &info_pattern);
    if (st != POLAR_STATUS_OK) {
        free(rate_pattern);
        free(Q_N);
        return st;
    }

    st = ca_polar_encoder_impl(
        a,
        A,
        CRC11_POLY,
        (int)sizeof(CRC11_POLY),
        info_pattern,
        N,
        rate_pattern,
        E,
        &e);
    if (st != POLAR_STATUS_OK) {
        free(rate_pattern);
        free(Q_N);
        free(info_pattern);
        return st;
    }

    st = get_3gpp_channel_interleaver_pattern_impl(E, &channel_pattern);
    if (st != POLAR_STATUS_OK) {
        free(rate_pattern);
        free(Q_N);
        free(info_pattern);
        free(e);
        return st;
    }

    st = alloc_u8((size_t)E, &out);
    if (st == POLAR_STATUS_OK) {
        for (int i = 0; i < E; i++) {
            out[i] = e[channel_pattern[i]];
        }
    }

    free(rate_pattern);
    free(Q_N);
    free(info_pattern);
    free(e);
    free(channel_pattern);

    if (st != POLAR_STATUS_OK) {
        return st;
    }

    f->data = out;
    f->len = E;
    return POLAR_STATUS_OK;
}

polar_status_t polar_custom1_decode(
    const double *f_tilde,
    int E,
    int A,
    int L,
    int min_sum,
    polar_u8_vec_t *a_hat
) {
    int K;
    int N;
    int *rate_pattern = NULL;
    int mode;
    int *Q_N = NULL;
    uint8_t *info_pattern = NULL;
    int *channel_pattern = NULL;
    double *e_tilde = NULL;
    uint8_t *out = NULL;
    int out_len = 0;
    polar_status_t st;

    if (f_tilde == NULL || a_hat == NULL || A <= 0 || E <= 0 || L <= 0) {
        return POLAR_STATUS_INVALID_ARGUMENT;
    }

    K = A + ((int)sizeof(CRC11_POLY) - 1);

    st = get_3gpp_n_impl(K, E, 1024, &N);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    st = get_3gpp_rate_matching_pattern_impl(K, N, E, &rate_pattern, &mode);
    if (st != POLAR_STATUS_OK) {
        return st;
    }

    st = get_pw_sequence_pattern_impl(N, &Q_N);
    if (st != POLAR_STATUS_OK) {
        free(rate_pattern);
        return st;
    }

    st = get_3gpp_info_bit_pattern_impl(K, Q_N, N, rate_pattern, E, mode, &info_pattern);
    if (st != POLAR_STATUS_OK) {
        free(rate_pattern);
        free(Q_N);
        return st;
    }

    st = get_3gpp_channel_interleaver_pattern_impl(E, &channel_pattern);
    if (st != POLAR_STATUS_OK) {
        free(rate_pattern);
        free(Q_N);
        free(info_pattern);
        return st;
    }

    st = alloc_f64((size_t)E, &e_tilde);
    if (st != POLAR_STATUS_OK) {
        free(rate_pattern);
        free(Q_N);
        free(info_pattern);
        free(channel_pattern);
        return st;
    }

    for (int i = 0; i < E; i++) {
        e_tilde[channel_pattern[i]] = f_tilde[i];
    }

    st = ca_polar_decoder_impl(
        e_tilde,
        E,
        CRC11_POLY,
        (int)sizeof(CRC11_POLY),
        info_pattern,
        N,
        rate_pattern,
        mode,
        L,
        min_sum,
        3,
        &out,
        &out_len);

    free(rate_pattern);
    free(Q_N);
    free(info_pattern);
    free(channel_pattern);
    free(e_tilde);

    if (st != POLAR_STATUS_OK) {
        return st;
    }

    a_hat->data = out;
    a_hat->len = out_len;
    return POLAR_STATUS_OK;
}

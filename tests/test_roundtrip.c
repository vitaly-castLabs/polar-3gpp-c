#include "polar_3gpp.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static uint32_t g_seed = 0x12345678u;

static uint32_t next_u32(void) {
    g_seed = g_seed * 1664525u + 1013904223u;
    return g_seed;
}

static void fill_random_bits(uint8_t *bits, int len) {
    int i;
    for (i = 0; i < len; i++) {
        bits[i] = (uint8_t)(next_u32() & 1u);
    }
}

static void bits_to_llr(const uint8_t *bits, int len, double *llr) {
    int i;
    for (i = 0; i < len; i++) {
        llr[i] = bits[i] ? -64.0 : 64.0;
    }
}

static int check_equal_bits(const uint8_t *a, const uint8_t *b, int len) {
    int i;
    for (i = 0; i < len; i++) {
        if (a[i] != b[i]) {
            return 0;
        }
    }
    return 1;
}

static int test_pbch(void) {
    uint8_t a[32];
    double llr[864];
    polar_u8_vec_t enc = {0};
    polar_u8_vec_t dec = {0};
    polar_status_t st;

    fill_random_bits(a, 32);

    st = polar_pbch_encode(a, 32, 864, &enc);
    if (st != POLAR_STATUS_OK) {
        fprintf(stderr, "PBCH encode failed: %d\n", st);
        return 1;
    }
    bits_to_llr(enc.data, enc.len, llr);
    st = polar_pbch_decode(llr, 864, 32, 8, 1, NULL, &dec);
    if (st != POLAR_STATUS_OK) {
        fprintf(stderr, "PBCH decode failed: %d\n", st);
        polar_free(enc.data);
        return 1;
    }
    if (dec.len != 32 || !check_equal_bits(a, dec.data, 32)) {
        fprintf(stderr, "PBCH mismatch\n");
        polar_free(enc.data);
        polar_free(dec.data);
        return 1;
    }

    polar_free(enc.data);
    polar_free(dec.data);
    return 0;
}

static int test_pdcch(int A, int E) {
    uint8_t *a = NULL;
    uint8_t rnti[16];
    double *llr = NULL;
    polar_u8_vec_t enc = {0};
    polar_u8_vec_t dec = {0};
    polar_status_t st;
    int fail = 0;

    a = (uint8_t *)malloc((size_t)A);
    llr = (double *)malloc((size_t)E * sizeof(double));
    if (a == NULL || llr == NULL) {
        fprintf(stderr, "PDCCH alloc failed\n");
        free(a);
        free(llr);
        return 1;
    }

    fill_random_bits(a, A);
    fill_random_bits(rnti, 16);

    st = polar_pdcch_encode(a, A, E, rnti, &enc);
    if (st != POLAR_STATUS_OK) {
        fprintf(stderr, "PDCCH encode failed: A=%d E=%d st=%d\n", A, E, st);
        fail = 1;
        goto done;
    }
    bits_to_llr(enc.data, enc.len, llr);
    st = polar_pdcch_decode(llr, E, A, 8, 1, rnti, &dec);
    if (st != POLAR_STATUS_OK) {
        fprintf(stderr, "PDCCH decode failed: A=%d E=%d st=%d\n", A, E, st);
        fail = 1;
        goto done;
    }
    if (dec.len != A || !check_equal_bits(a, dec.data, A)) {
        fprintf(stderr, "PDCCH mismatch: A=%d E=%d dec_len=%d\n", A, E, dec.len);
        fail = 1;
        goto done;
    }

done:
    free(a);
    free(llr);
    polar_free(enc.data);
    polar_free(dec.data);
    return fail;
}

static int test_pucch(int A, int G) {
    uint8_t *a = NULL;
    double *llr = NULL;
    polar_u8_vec_t enc = {0};
    polar_u8_vec_t dec = {0};
    polar_status_t st;
    int fail = 0;

    a = (uint8_t *)malloc((size_t)A);
    llr = (double *)malloc((size_t)G * sizeof(double));
    if (a == NULL || llr == NULL) {
        fprintf(stderr, "PUCCH alloc failed\n");
        free(a);
        free(llr);
        return 1;
    }
    fill_random_bits(a, A);

    st = polar_pucch_encode(a, A, G, &enc);
    if (st != POLAR_STATUS_OK) {
        fprintf(stderr, "PUCCH encode failed: A=%d G=%d st=%d\n", A, G, st);
        fail = 1;
        goto done;
    }
    bits_to_llr(enc.data, enc.len, llr);
    st = polar_pucch_decode(llr, G, A, 8, 1, &dec);
    if (st != POLAR_STATUS_OK) {
        fprintf(stderr, "PUCCH decode failed: A=%d G=%d st=%d\n", A, G, st);
        fail = 1;
        goto done;
    }
    if (dec.len != A || !check_equal_bits(a, dec.data, A)) {
        fprintf(stderr, "PUCCH mismatch: A=%d G=%d dec_len=%d\n", A, G, dec.len);
        fail = 1;
        goto done;
    }

done:
    free(a);
    free(llr);
    polar_free(enc.data);
    polar_free(dec.data);
    return fail;
}

static int test_custom1(void) {
    enum { A = 96, E = 256 };
    uint8_t a[A];
    double llr[E];
    polar_u8_vec_t enc = {0};
    polar_u8_vec_t dec = {0};
    polar_status_t st;

    fill_random_bits(a, A);

    st = polar_custom1_encode(a, A, E, &enc);
    if (st != POLAR_STATUS_OK) {
        fprintf(stderr, "custom1 encode failed: %d\n", st);
        return 1;
    }
    bits_to_llr(enc.data, enc.len, llr);
    st = polar_custom1_decode(llr, E, A, 8, 1, &dec);
    if (st != POLAR_STATUS_OK) {
        fprintf(stderr, "custom1 decode failed: %d\n", st);
        polar_free(enc.data);
        return 1;
    }
    if (dec.len != A || !check_equal_bits(a, dec.data, A)) {
        fprintf(stderr, "custom1 mismatch\n");
        polar_free(enc.data);
        polar_free(dec.data);
        return 1;
    }

    polar_free(enc.data);
    polar_free(dec.data);
    return 0;
}

int main(void) {
    int failures = 0;

    failures += test_pbch();
    failures += test_pdcch(20, 432);
    failures += test_pdcch(8, 216);
    failures += test_pucch(16, 108);
    failures += test_pucch(64, 256);
    failures += test_pucch(361, 1088);
    failures += test_custom1();

    if (failures == 0) {
        puts("All roundtrip tests passed");
        return 0;
    }
    fprintf(stderr, "Roundtrip tests failed: %d\n", failures);
    return 1;
}

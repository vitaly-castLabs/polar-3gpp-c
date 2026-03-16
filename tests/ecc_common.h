#ifndef ECC_COMMON_H
#define ECC_COMMON_H

#include <stdint.h>
#include <stdlib.h>

static uint32_t ecc_seed = 0xC0DEC0DEu;

static uint32_t ecc_next_u32(void) {
    ecc_seed = ecc_seed * 1664525u + 1013904223u;
    return ecc_seed;
}

static int ecc_rand_range(int min_inclusive, int max_inclusive) {
    uint32_t span = (uint32_t)(max_inclusive - min_inclusive + 1);
    return min_inclusive + (int)(ecc_next_u32() % span);
}

static void ecc_fill_random_bits(uint8_t *bits, int len) {
    int i;
    for (i = 0; i < len; i++) {
        bits[i] = (uint8_t)(ecc_next_u32() & 1u);
    }
}

static void ecc_bits_to_llr(const uint8_t *bits, int len, double *llr) {
    int i;
    for (i = 0; i < len; i++) {
        llr[i] = bits[i] ? -64.0 : 64.0;
    }
}

static int ecc_equal_bits(const uint8_t *a, const uint8_t *b, int len) {
    int i;
    for (i = 0; i < len; i++) {
        if (a[i] != b[i]) {
            return 0;
        }
    }
    return 1;
}

static void ecc_shuffle_i32(int *arr, int len) {
    int i;
    for (i = len - 1; i > 0; i--) {
        int j = (int)(ecc_next_u32() % (uint32_t)(i + 1));
        int tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }
}

#endif

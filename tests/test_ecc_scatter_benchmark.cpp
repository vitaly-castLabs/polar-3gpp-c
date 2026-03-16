#include "polar_3gpp.h"

#include <benchmark/benchmark.h>

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <random>
#include <string>
#include <vector>

namespace {

constexpr int kA = 13;
constexpr int kE = 8040;
constexpr int kL = 8;
constexpr int kMinSum = 1;

void FillPseudoRandomBits(std::vector<uint8_t> &bits, uint32_t seed) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> bit_dist(0, 1);
    for (auto &b : bits) {
        b = static_cast<uint8_t>(bit_dist(rng));
    }
}

void BitsToLlr(const uint8_t *bits, int len, std::vector<double> &llr) {
    llr.resize(static_cast<size_t>(len));
    for (int i = 0; i < len; ++i) {
        llr[static_cast<size_t>(i)] = bits[i] ? -64.0 : 64.0;
    }
}

static void BM_PucchEncodeScatterCase(benchmark::State &state) {
    std::vector<uint8_t> a(static_cast<size_t>(kA));
    FillPseudoRandomBits(a, 1234u);

    for (auto _ : state) {
        polar_u8_vec_t enc{nullptr, 0};
        const polar_status_t st = polar_pucch_encode(a.data(), kA, kE, &enc);
        if (st != POLAR_STATUS_OK) {
            state.SkipWithError("polar_pucch_encode failed");
            return;
        }

        benchmark::DoNotOptimize(enc.data);
        benchmark::DoNotOptimize(enc.len);
        polar_free(enc.data);
    }

    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(kE));
}

static void BM_PucchDecodeScatterCase(benchmark::State &state) {
    const int preserved_bits = static_cast<int>(state.range(0));

    std::vector<uint8_t> a(static_cast<size_t>(kA));
    FillPseudoRandomBits(a, 5678u + static_cast<uint32_t>(preserved_bits));

    polar_u8_vec_t enc{nullptr, 0};
    const polar_status_t enc_st = polar_pucch_encode(a.data(), kA, kE, &enc);
    if (enc_st != POLAR_STATUS_OK) {
        state.SkipWithError("polar_pucch_encode setup failed");
        return;
    }

    std::vector<double> llr_full;
    BitsToLlr(enc.data, enc.len, llr_full);

    std::vector<double> llr_work = llr_full;
    std::vector<int> indices(static_cast<size_t>(kE));
    std::iota(indices.begin(), indices.end(), 0);
    std::mt19937 rng(9001u + static_cast<uint32_t>(preserved_bits));

    for (auto _ : state) {
        state.PauseTiming();
        llr_work = llr_full;
        std::shuffle(indices.begin(), indices.end(), rng);
        for (int j = preserved_bits; j < kE; ++j) {
            llr_work[static_cast<size_t>(indices[static_cast<size_t>(j)])] = 0.0;
        }
        state.ResumeTiming();

        polar_u8_vec_t dec{nullptr, 0};
        const polar_status_t dec_st = polar_pucch_decode(llr_work.data(), kE, kA, kL, kMinSum, &dec);
        if (dec_st != POLAR_STATUS_OK) {
            polar_free(enc.data);
            state.SkipWithError("polar_pucch_decode failed");
            return;
        }

        benchmark::DoNotOptimize(dec.data);
        benchmark::DoNotOptimize(dec.len);
        polar_free(dec.data);
    }

    polar_free(enc.data);
    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(kE));
}

}  // namespace

BENCHMARK(BM_PucchEncodeScatterCase);
BENCHMARK(BM_PucchDecodeScatterCase)->DenseRange(15, 30, 1);

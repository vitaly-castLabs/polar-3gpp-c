# polar-3gpp-matlab (C port)

C implementation of the Polar encoder/decoder flow from the original MATLAB project:
https://github.com/robmaunder/polar-3gpp-matlab

This project provides a reusable static library (`polar_3gpp`) plus small test programs for PBCH/PDCCH/PUCCH and custom ECC experiments.

## Build

Requirements:
- CMake >= 3.16
- C compiler with C99 support
- C++ compiler with C++17 support (for optional benchmark target)

Configure and build:

```bash
cmake -S . -B build
cmake --build build -j
```

## Run tests

```bash
ctest --test-dir build --output-on-failure
```

## Run simulations

Examples:

```bash
./build/test_roundtrip
./build/test_ecc
./build/test_ecc_crop
./build/test_ecc_scatter
./build/test_ecc_scatter_custom1
./build/test_ecc_scatter_pdcch
```

Optional arguments for `test_ecc*` binaries:
- arg1: number of iterations
- arg2: random seed

Example with optional arguments:

```bash
./build/test_ecc 5 1234
```

## Benchmarking

`test_ecc_scatter_bench` requires Google Benchmark and is disabled by default.

```bash
cmake -S . -B build -DPOLAR_BUILD_BENCHMARK=ON
cmake --build build -j
./build/test_ecc_scatter_bench
```

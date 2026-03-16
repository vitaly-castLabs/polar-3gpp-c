cmake_minimum_required(VERSION 3.14)

# Clang FTW
set(CMAKE_C_COMPILER_NAMES clang gcc icc cc)
set(CMAKE_CXX_COMPILER_NAMES clang++ g++ icpc c++ cxx)

project(polar_3gpp_c C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# this is only required for the benchmark target
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(polar_3gpp STATIC
    src/polar_3gpp.c
)

target_include_directories(polar_3gpp
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(polar_3gpp PUBLIC m)

add_executable(test_roundtrip
    tests/test_roundtrip.c
)

target_link_libraries(test_roundtrip PRIVATE polar_3gpp m)

add_executable(test_ecc
    tests/test_ecc.c
)
target_link_libraries(test_ecc PRIVATE polar_3gpp m)

add_executable(test_ecc_crop
    tests/test_ecc_crop.c
)
target_link_libraries(test_ecc_crop PRIVATE polar_3gpp m)

add_executable(test_ecc_scatter
    tests/test_ecc_scatter.c
)
target_link_libraries(test_ecc_scatter PRIVATE polar_3gpp m)

add_executable(test_ecc_scatter_custom1
    tests/test_ecc_scatter_custom1.c
)
target_link_libraries(test_ecc_scatter_custom1 PRIVATE polar_3gpp m)

add_executable(test_ecc_scatter_pdcch
    tests/test_ecc_scatter_pdcch.c
)
target_link_libraries(test_ecc_scatter_pdcch PRIVATE polar_3gpp m)

enable_testing()
add_test(NAME roundtrip COMMAND test_roundtrip)
add_test(NAME ecc_smoke COMMAND test_ecc 2 1234)
add_test(NAME ecc_crop_smoke COMMAND test_ecc_crop 5 1234)
add_test(NAME ecc_scatter_smoke COMMAND test_ecc_scatter 5 1234)
add_test(NAME ecc_scatter_custom1_smoke COMMAND test_ecc_scatter_custom1 5 1234)
add_test(NAME ecc_scatter_pdcch_smoke COMMAND test_ecc_scatter_pdcch 5 1234)

option(POLAR_BUILD_BENCHMARK "Build benchmark target (requires Google Benchmark)" OFF)

if(POLAR_BUILD_BENCHMARK)
    find_package(benchmark CONFIG REQUIRED)

    add_executable(test_ecc_scatter_bench
        tests/test_ecc_scatter_benchmark.cpp
    )
    target_link_libraries(test_ecc_scatter_bench PRIVATE polar_3gpp benchmark::benchmark benchmark::benchmark_main m)
endif()

#pragma once

#include <stdint.h>
#include <stddef.h>

#define FFT_USE_LINEARITY_TRICK() 1

#define FFT_USE_BUILTIN_UNPREDICTABLE() 1

#if UINTPTR_MAX == 0xFFFFFFFFFFFFFFFFul && !defined(FFT_FORCE_32_BIT)
#   define FFT_LIMB_BITS    64
typedef uint64_t            FFT_ULIMB;
typedef int64_t             FFT_SLIMB;
typedef unsigned __int128   FFT_DOUBLE_ULIMB;
typedef __int128            FFT_DOUBLE_SLIMB;
#else
#   define FFT_LIMB_BITS    32
typedef uint32_t            FFT_ULIMB;
typedef int32_t             FFT_SLIMB;
typedef uint64_t            FFT_DOUBLE_ULIMB;
typedef int64_t             FFT_DOUBLE_SLIMB;
#endif

typedef struct {
    size_t ndigits;
    int base_log;
} FFT_Threshold;

extern const FFT_Threshold FFT_THRESHOLDS[];

void fft(
        FFT_ULIMB *a1, FFT_ULIMB *b1,
        FFT_ULIMB *a2, FFT_ULIMB *b2,
        FFT_ULIMB *scratch,
        size_t n);

void fft_fourstep(
        FFT_ULIMB *a1, FFT_ULIMB *b1,
        FFT_ULIMB *a2, FFT_ULIMB *b2,
        FFT_ULIMB *scratch,
        size_t n);

void fft_sixstep(
        FFT_ULIMB *a1, FFT_ULIMB *b1,
        FFT_ULIMB *a2, FFT_ULIMB *b2,
        FFT_ULIMB *scratch,
        size_t n);

void fft_recover_answer(
        FFT_ULIMB *out, size_t nout,
        FFT_ULIMB *a1, FFT_ULIMB *a2,
        int base_log);

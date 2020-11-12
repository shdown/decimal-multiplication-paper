#pragma once

#include <stddef.h>

#define FFT_UNUSED          __attribute__((unused))
#define FFT_FORCE_INLINE    __attribute__((always_inline))
#define FFT_UNREACHABLE()   __builtin_unreachable()

static inline FFT_UNUSED FFT_FORCE_INLINE
unsigned fft_zu_counttz(size_t x)
{
# if __SIZEOF_INT__ == __SIZEOF_SIZE_T__
    return __builtin_ctz(x);
# elif __SIZEOF_LONG__ == __SIZEOF_SIZE_T__
    return __builtin_ctzl(x);
# elif __SIZEOF_LONG_LONG__ == __SIZEOF_SIZE_T__
    return __builtin_ctzll(x);
# else
#  error "Unsupported platform."
# endif
}

static inline FFT_UNUSED FFT_FORCE_INLINE
unsigned fft_zu_countlz(size_t x)
{
# if __SIZEOF_INT__ == __SIZEOF_SIZE_T__
    return __builtin_clz(x);
# elif __SIZEOF_LONG__ == __SIZEOF_SIZE_T__
    return __builtin_clzl(x);
# elif __SIZEOF_LONG_LONG__ == __SIZEOF_SIZE_T__
    return __builtin_clzll(x);
# else
#  error "Unsupported platform."
# endif
}

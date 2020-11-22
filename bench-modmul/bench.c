// (c) 2020 shdown
// This code is licensed under MIT license (see LICENSE.MIT for details)

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define ULIMB uint64_t
#define SLIMB int64_t
#define DOUBLE_ULIMB unsigned __int128
#define LIMB_BITS 64
#define FORCE_INLINE __attribute__((always_inline))

static inline FORCE_INLINE
DOUBLE_ULIMB full_prod(ULIMB x, ULIMB y)
{
    return ((DOUBLE_ULIMB) x) * y;
}

#if defined(USE_NAIVE)

#define MODULO 18446744069414584321ull
static inline FORCE_INLINE
ULIMB mul(ULIMB x, ULIMB y)
{
    DOUBLE_ULIMB t = full_prod(x, y);
    ULIMB t_hi = t >> LIMB_BITS;
    ULIMB t_lo = t;

    ULIMB q, r;

    asm inline (
        "divq %4\n"

        : /*0*/ "=a" (q)
        , /*1*/ "=d" (r)

        : /*2*/ "0" (t_lo)
        , /*3*/ "1" (t_hi)
        , /*4*/ "r" (MODULO)

        : "cc");

    return r;
}

#elif defined(USE_SOLINAS)

#define MODULO 18446744069414584321ull
static inline FORCE_INLINE
ULIMB mul(ULIMB x, ULIMB y)
{
    DOUBLE_ULIMB t = full_prod(x, y);

#define ROUND() \
    do { \
        ULIMB hi = t >> LIMB_BITS; \
        ULIMB lo = t; \
        t = (((DOUBLE_ULIMB) hi) << 32) - hi + lo; \
    } while (0)

    ROUND();
    ROUND();

#undef ROUND

    // no adjustment
    return t;
}

#elif defined(USE_MONT)

#define MODULO    9203387313508319233ull
#define MODULO_P1 9203387313508319231ull
static inline FORCE_INLINE
ULIMB mul(ULIMB x, ULIMB y)
{
    DOUBLE_ULIMB t = full_prod(x, y);
    ULIMB m = t * MODULO_P1;
    DOUBLE_ULIMB mp = full_prod(m, MODULO);
    ULIMB r = (t + mp) >> LIMB_BITS;

    // no adjustment
    return r;
}

#else
#   error "Please define either of: USE_NAIVE, USE_SOLINAS, USE_MONT."
#endif

static int parse_int_or_die(const char *s)
{
    int r;
    if (sscanf(s, "%d", &r) != 1) {
        fprintf(stderr, "Cannot parse '%s' as integer.\n", s);
        abort();
    }
    return r;
}

int main(int argc, char **argv)
{
    if (argc != 3) {
        fprintf(stderr, "USAGE: %s seed nrepeat\n", argv[0]);
        abort();
    }
    ULIMB x = parse_int_or_die(argv[1]);
    int n = parse_int_or_die(argv[2]);
    for (int i = 0; i < n; ++i) {
        x = mul(x, x);
        x = mul(x, x);
        x = mul(x, x);
        x = mul(x, x);
        x = mul(x, x);
        x = mul(x, x);
        x = mul(x, x);
        x = mul(x, x);
        x = mul(x, x);
        x = mul(x, x);
    }
    printf("%u\n", (unsigned) x);
}

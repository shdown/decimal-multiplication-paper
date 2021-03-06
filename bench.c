// (c) 2020 shdown
// This code is licensed under MIT license (see LICENSE.MIT for details)

#include "fft.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

static void die_oom(void)
{
    fprintf(stderr, "Out of memory.\n");
    abort();
}

static inline size_t x_add_zu(size_t a, size_t b)
{
    size_t r = a + b;
    if (r < a)
        die_oom();
    return r;
}

static inline size_t x_mul_zu(size_t a, size_t b)
{
    if (b && a > SIZE_MAX / b)
        die_oom();
    return a * b;
}

static void *xcalloc(size_t n, size_t m)
{
    void *r = calloc(n, m);
    if (!r && n && m)
        die_oom();
    return r;
}

static size_t div_ceil_zu(size_t x, size_t y)
{
    return x / y + !!(x % y);
}

static void write_arg(FFT_ULIMB *p, size_t n, size_t ndigits, unsigned base_log)
{
    FFT_ULIMB base = 1;
    for (unsigned i = 0; i < base_log; ++i)
        base *= 10;

    FFT_ULIMB max_limb = base - 1;

    size_t nlimbs = div_ceil_zu(ndigits, base_log);
    for (size_t i = 0; i < nlimbs; ++i)
        p[i] = max_limb;
    for (size_t i = nlimbs; i < n; ++i)
        p[i] = 0;
}

static size_t x_fft_size_aux(size_t init_factor, size_t m)
{
    size_t n = init_factor;
    while (n < m)
        n = x_mul_zu(n, 2);
    return n;
}

static size_t x_get_fft_size(unsigned *method, size_t m)
{
    switch (*method) {
    case 1:
    case 6:
        return x_fft_size_aux(1, m);
    case 4:
        return x_fft_size_aux(3, m);
    case 0:
        {
            size_t n_simple = x_fft_size_aux(1, m);
            size_t n_fourstep = x_fft_size_aux(3, m);
            if (n_simple < n_fourstep) {
                *method = 1;
                return n_simple;
            } else {
                *method = 4;
                return n_fourstep;
            }
        }
    default:
        fprintf(stderr, "Unknown method: %u (expected either of: 0 1 4 6)\n", *method);
        abort();
    }
}

static int max_base_log(size_t ndigits)
{
    for (const FFT_Threshold *p = FFT_THRESHOLDS; ; ++p) {
        if (!p->ndigits)
            return -1;
        if (ndigits <= p->ndigits)
            return p->base_log;
    }
}

static size_t x_parse_zu(const char *s)
{
    size_t r;
    if (sscanf(s, "%zu", &r) != 1) {
        fprintf(stderr, "Cannot parse integer: '%s'.\n", s);
        abort();
    }
    return r;
}

int main(int argc, char **argv)
{
    if (argc != 4) {
        fprintf(stderr, "USAGE: bench METHOD NDIGITS NREPEAT\n");
        return 2;
    }
    unsigned method = x_parse_zu(argv[1]);
    size_t ndigits = x_parse_zu(argv[2]);
    size_t nrepeat = x_parse_zu(argv[3]);

    int base_log = max_base_log(ndigits);
    if (base_log < 0) {
        fprintf(stderr, "Transform of such large size is not supported.\n");
        abort();
    }

    size_t nlimbs_arg = div_ceil_zu(ndigits, base_log);
    size_t nlimbs_r = x_mul_zu(nlimbs_arg, 2);
    size_t n = x_get_fft_size(&method, nlimbs_r);

    fprintf(stderr, "method=%u, ndigits=%zu, nrepeat=%zu | base=10^%d, fft_size=%zu\n",
            method, ndigits, nrepeat, base_log, n);

    FFT_ULIMB *mem = xcalloc(sizeof(FFT_ULIMB), x_add_zu(4, x_mul_zu(n, 5)));

    FFT_ULIMB *a1      = mem + 0 * n + 0;
    FFT_ULIMB *b1      = mem + 1 * n + 1;
    FFT_ULIMB *a2      = mem + 2 * n + 2;
    FFT_ULIMB *b2      = mem + 3 * n + 3;
    FFT_ULIMB *scratch = mem + 4 * n + 4;

    clock_t t = clock();

    for (size_t i = 0; i < nrepeat; ++i) {
        write_arg(a1, n, ndigits, base_log);
        write_arg(b1, n, ndigits, base_log);
        memcpy(a2, a1, sizeof(FFT_ULIMB) * n);
        memcpy(b2, b1, sizeof(FFT_ULIMB) * n);

        switch (method) {
        case 1:
            fft(a1, b1, a2, b2, scratch, n);
            break;
        case 4:
            fft_fourstep(a1, b1, a2, b2, scratch, n);
            break;
        case 6:
            fft_sixstep(a1, b1, a2, b2, scratch, n);
            break;
        default:
            abort();
            break;
        }

        fft_recover_answer(a1, nlimbs_r, a1, a2, base_log);
    }

    t = clock() - t;
    printf("%.4f\n", ((double) t) / CLOCKS_PER_SEC);

    free(mem);
    return 0;
}

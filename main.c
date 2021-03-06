// (c) 2020 shdown
// This code is licensed under MIT license (see LICENSE.MIT for details)

#include "fft.h"
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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

static void *x_calloc(size_t n, size_t m)
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

static size_t x_calc_fft_size(unsigned method, size_t m)
{
    size_t n;
    switch (method) {
    case 1:
    case 6:
        n = 1;
        break;
    case 4:
        n = 3;
        break;
    default:
        fprintf(stderr, "Unknown method, expected either of: 1, 4, 6.\n");
        abort();
    }

    while (n < m)
        n = x_mul_zu(n, 2);
    return n;
}

static char *read_number_str(void)
{
    char *buf = NULL;
    size_t nbuf = 0;
    ssize_t r = getline(&buf, &nbuf, stdin);
    if (r < 0) {
        if (ferror(stdin)) {
            perror("getline");
        } else {
            fprintf(stderr, "Got EOF.\n");
        }
        abort();
    }

    if (r && buf[r - 1] == '\n') {
        --r;
        buf[r] = '\0';
    }

    if (!r) {
        fprintf(stderr, "Got empty line.\n");
        abort();
    }
    for (size_t i = 0; i < (size_t) r; ++i) {
        unsigned char c = buf[i];
        if (c < '0' || c > '9') {
            fprintf(stderr, "Got non-numeric data.\n");
            abort();
        }
    }

    return buf;
}

static FFT_ULIMB parse_limb(const char *s, size_t ns)
{
    FFT_ULIMB r = 0;
    for (size_t i = 0; i < ns; ++i) {
        r *= 10;
        r += s[i] - '0';
    }
    return r;
}

static void parse_str(FFT_ULIMB *out, const char *s, size_t ns, unsigned base_log)
{
    size_t i = ns;
    while (i >= base_log) {
        i -= base_log;
        *out++ = parse_limb(s + i, base_log);
    }
    if (i) {
        *out++ = parse_limb(s, i);
    }
}

#if FFT_LIMB_BITS == 32
#   define FFT_ULIMB_FMT PRIu32
#elif FFT_LIMB_BITS == 64
#   define FFT_ULIMB_FMT PRIu64
#else
#   error "Unsupported FFT_LIMB_BITS."
#endif

static void pretty_print(const FFT_ULIMB *w, size_t nw, int base_log)
{
    // normalize
    while (nw && w[nw - 1] == 0)
        --nw;

    if (!nw) {
        printf("0\n");
        return;
    }

    size_t i = nw - 1;
    printf("%" FFT_ULIMB_FMT, w[i]);
    while (i) {
        --i;
        printf("%0*" FFT_ULIMB_FMT, base_log, w[i]);
    }
    printf("\n");
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

int main(int argc, char **argv)
{
    if (argc != 2) {
        fprintf(stderr, "USAGE: main METHOD\n");
        fprintf(stderr, "  where METHOD must be either of: 1 4 6.\n");
        fprintf(stderr, "METHOD of 1 means use the 'straight' FFT (function 'fft()').\n");
        fprintf(stderr, "METHOD of 4 means use the four-step FFT (function 'fft_fourstep()').\n");
        fprintf(stderr, "METHOD of 6 means use the six-step FFT (function 'fft_sixstep()').\n");
        return 2;
    }
    unsigned method;
    if (sscanf(argv[1], "%u", &method) != 1) {
        fprintf(stderr, "Cannot parse '%s' as unsigned integer (METHOD).\n", argv[1]);
        return 2;
    }

    char *a = read_number_str();
    size_t na = strlen(a);

    char *b = read_number_str();
    size_t nb = strlen(b);

    int base_log = max_base_log(na < nb ? na : nb);
    if (base_log < 0) {
        fprintf(stderr, "Transform of such large size is not supported.\n");
        abort();
    }

    size_t nlimbs_a = div_ceil_zu(na, base_log);
    size_t nlimbs_b = div_ceil_zu(nb, base_log);
    size_t nlimbs_r = x_add_zu(nlimbs_a, nlimbs_b);
    size_t n = x_calc_fft_size(method, nlimbs_r);

    FFT_ULIMB *mem = x_calloc(sizeof(FFT_ULIMB), x_add_zu(4, x_mul_zu(n, 5)));

    FFT_ULIMB *a1      = mem + 0 * n + 0;
    FFT_ULIMB *b1      = mem + 1 * n + 1;
    FFT_ULIMB *a2      = mem + 2 * n + 2;
    FFT_ULIMB *b2      = mem + 3 * n + 3;
    FFT_ULIMB *scratch = mem + 4 * n + 4;

    parse_str(a1, a, na, base_log);
    memcpy(a2, a1, sizeof(FFT_ULIMB) * n);

    parse_str(b1, b, nb, base_log);
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
    }

    fft_recover_answer(a1, nlimbs_r, a1, a2, base_log);

    pretty_print(a1, nlimbs_r, base_log);

    free(mem);
    free(a);
    free(b);
    return 0;
}

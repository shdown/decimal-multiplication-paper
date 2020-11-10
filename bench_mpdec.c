#include <mpdecimal.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static mpd_t *xalloc_mpd_t(void)
{
    mpd_t *r = mpd_qnew();
    if (!r) {
        fprintf(stderr, "Out of memory.\n");
        abort();
    }
    return r;
}

static int x_parse_int(const char *s)
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
        fprintf(stderr, "USAGE: bench_mpdec NDIGITS NREPEAT\n");
        return 2;
    }
    int ndigits = x_parse_int(argv[1]);
    int nrepeat = x_parse_int(argv[2]);

    mpd_t *a = xalloc_mpd_t();
    mpd_t *b = xalloc_mpd_t();
    mpd_t *c = xalloc_mpd_t();

    mpd_context_t ctx;
    mpd_maxcontext(&ctx);

    mpd_set_uint(a, 10, &ctx);
    mpd_set_uint(b, ndigits, &ctx);

    mpd_pow(a, a, b, &ctx);

    mpd_sub_uint(a, a, 1, &ctx);

    mpd_sub_uint(b, a, 1, &ctx);

    fprintf(stderr, "ndigits=%d, nrepeat=%d\n", ndigits, nrepeat);
    clock_t t = clock();

    for (int i = 0; i < nrepeat; ++i) {
        mpd_mul(c, a, b, &ctx);
    }

    t = clock() - t;
    printf("%.4f\n", ((double) t) / CLOCKS_PER_SEC);

    mpd_del(a);
    mpd_del(b);
    mpd_del(c);
}

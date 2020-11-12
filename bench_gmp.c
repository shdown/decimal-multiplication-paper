#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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
    if (argc != 4) {
        fprintf(stderr, "USAGE: bench_gmp BASE NDIGITS NREPEAT\n");
        return 2;
    }
    int base = x_parse_int(argv[1]);
    int ndigits = x_parse_int(argv[2]);
    int nrepeat = x_parse_int(argv[3]);

    mpz_t x;
    mpz_init(x);
    mpz_ui_pow_ui(x, base, ndigits);
    mpz_sub_ui(x, x, 1);

    mpz_t y;
    mpz_init_set(y, x);
    mpz_sub_ui(y, y, 1);

    mpz_t z;
    mpz_init(z);

    fprintf(stderr, "base=%d, ndigits=%d, nrepeat=%d\n", base, ndigits, nrepeat);
    clock_t t = clock();

    for (int i = 0; i < nrepeat; ++i) {
        mpz_mul(z, x, y);
    }

    t = clock() - t;
    printf("%.4f\n", ((double) t) / CLOCKS_PER_SEC);

    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(z);
    return 0;
}

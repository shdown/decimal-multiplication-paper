#include <mpdecimal.h>
#include <stdio.h>

enum {
	NDIGITS = 15000000,
	NREPEAT = 5,
};

static mpd_t *xalloc_mpd_t(void)
{
	mpd_t *r = mpd_qnew();
	if (!r) {
		fprintf(stderr, "Out of memory.\n");
		abort();
	}
	return r;
}

int main()
{
	mpd_t *a = xalloc_mpd_t();
	mpd_t *b = xalloc_mpd_t();
	mpd_t *c = xalloc_mpd_t();

	mpd_context_t ctx;
	mpd_maxcontext(&ctx);

	mpd_set_uint(a, 10, &ctx);
	mpd_set_uint(b, NDIGITS, &ctx);

	mpd_pow(a, a, b, &ctx);

	mpd_sub_uint(a, a, 1, &ctx);

	mpd_sub_uint(b, a, 1, &ctx);

	for (int i = 0; i < NREPEAT; ++i)
		mpd_mul(c, a, b, &ctx);

	mpd_del(a);
	mpd_del(b);
	mpd_del(c);
}

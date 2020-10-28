#include "fft_utils.h"
#include "fft_compdep.h"

#define SWAP(Type_, X_, Y_) \
    do { \
        Type_ swap_tmp__ = (X_); \
        (X_) = (Y_); \
        (Y_) = swap_tmp__; \
    } while (0)

enum { SIDE = 64 };

void fftu_permute(FFT_ULIMB *x, size_t n)
{
    unsigned order = fft_zu_counttz(n);

    if (!order)
        return;

    unsigned shift_r = sizeof(size_t) * 8 - order;

    size_t i = 0;
    size_t rev_i = 0;
    while (i < n) {
        // swap if needed
        if (i < rev_i) {
            SWAP(FFT_ULIMB, x[i], x[rev_i]);
        }
        // advance 'i' and update 'rev_i'
        size_t tail = i ^ (i + 1);
        tail <<= fft_zu_countlz(tail);
        tail >>= shift_r;
        rev_i ^= tail;

        ++i;
    }
}

static inline void transpose_square_simple(FFT_ULIMB *m, size_t side, size_t actual_width)
{
    FFT_ULIMB *row_i = m;
    for (size_t i = 0; i < side; ++i) {
        FFT_ULIMB *row_j = row_i + actual_width;
        for (size_t j = i + 1; j < side; ++j) {
            SWAP(FFT_ULIMB, row_i[j], row_j[i]);
            row_j += actual_width;
        }
        row_i += actual_width;
    }
}

static inline FFT_FORCE_INLINE
void copy_in(FFT_ULIMB *buf, const FFT_ULIMB *m, size_t actual_width)
{
    for (int i = 0; i < SIDE; ++i) {
        for (int j = 0; j < SIDE; ++j)
            buf[j] = m[j];
        m += actual_width;
        buf += SIDE;
    }
}

static inline FFT_FORCE_INLINE
void copy_out(const FFT_ULIMB *buf, FFT_ULIMB *m, size_t actual_width)
{
    for (int i = 0; i < SIDE; ++i) {
        for (int j = 0; j < SIDE; ++j)
            m[j] = buf[j];
        m += actual_width;
        buf += SIDE;
    }
}

static void transpose_square(FFT_ULIMB *m, size_t side, size_t actual_width, FFT_ULIMB *scratch)
{
    if (side <= SIDE) {
        transpose_square_simple(m, side, actual_width);
        return;
    }

    FFT_ULIMB *scratch1 = scratch;
    FFT_ULIMB *scratch2 = scratch + SIDE * SIDE;

    size_t stride = actual_width * SIDE;

    FFT_ULIMB *row_r = m;
    for (size_t r = 0; r < side; r += SIDE) {

        copy_in(scratch1, row_r + r, actual_width);
        transpose_square_simple(scratch1, SIDE, SIDE);
        copy_out(scratch1, row_r + r, actual_width);

        FFT_ULIMB *row_c = row_r + stride;
        for (size_t c = r + SIDE; c < side; c += SIDE) {
            copy_in(scratch1, row_r + c, actual_width);
            transpose_square_simple(scratch1, SIDE, SIDE);

            copy_in(scratch2, row_c + r, actual_width);
            transpose_square_simple(scratch2, SIDE, SIDE);

            copy_out(scratch1, row_c + r, actual_width);

            copy_out(scratch2, row_r + c, actual_width);

            row_c += stride;
        }
        row_r += stride;
    }
}

static inline void special_permute_row(FFT_ULIMB *x, size_t n, FFT_ULIMB *scratch)
{
    for (size_t i = 0; i < n; ++i)
        scratch[i] = x[i];

    FFT_ULIMB *y = x + n / 2;
    for (size_t i = 0; i < n; i += 2) {
        x[i / 2] = scratch[i    ];
        y[i / 2] = scratch[i + 1];
    }
}

static inline void special_unpermute_row(FFT_ULIMB *x, size_t n, FFT_ULIMB *scratch)
{
    for (size_t i = 0; i < n; ++i)
        scratch[i] = x[i];

    FFT_ULIMB *scratch_hi = scratch + n / 2;
    for (size_t i = 0; i < n; i += 2) {
        x[i    ] = scratch   [i / 2];
        x[i + 1] = scratch_hi[i / 2];
    }
}

void fftu_transpose(FFT_ULIMB *x, size_t r, size_t c, FFT_ULIMB *scratch)
{
    if (r == c) {
        // This is 2^k by 2^k matrix.
        transpose_square(x, r, c, scratch);

    } else if (r < c) {
        // This is 2^k by 2^(k+1) matrix.
        size_t n = r * c;
        for (size_t i = 0; i < n; i += c)
            special_permute_row(x + i, c, scratch);

        transpose_square(x + 0, r, c, scratch);

        transpose_square(x + r, r, c, scratch);

    } else {
        // This is 2^(k+1) by 2^k matrix.
        transpose_square(x + 0, c, r, scratch);

        transpose_square(x + c, c, r, scratch);

        size_t n = c * r;
        for (size_t i = 0; i < n; i += r)
            special_unpermute_row(x + i, r, scratch);
    }
}

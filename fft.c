#include "fft.h"
#include "fft_utils.h"
#include "fft_compdep.h"

// R = 2^FFT_LIMB_BITS.
//
// M(x) is Montgomery representation of x: R*x % p.
//
// All primitive roots are taken in F_p, x^-1 means inverse of x in F_p.
typedef struct {
    // Prime of form c*3*2^k+1.
    FFT_ULIMB p;

    // The 'k' from the formula above.
    unsigned k;

    // The unique integer in [1; R) such that mont_p1 * p = -1 (mod R).
    FFT_ULIMB mont_p1;

    // Value of R modulo p.
    FFT_ULIMB mont_rmodp;

    // Value of R^2 modulo p.
    FFT_ULIMB mont_r2modp;

    // Let 'rK' be a primitive root of order 2^k.
    // Then 'root_M' is M(rK), 'inv_root_M' is M(rK^-1).
    FFT_ULIMB root_M;
    FFT_ULIMB inv_root_M;

    // Let 'rT' be a primitive root of order 3.
    // Then 'root3_M' is M(rT), 'inv_root3_M' is M(rT^-1).
    FFT_ULIMB root3_M;
    FFT_ULIMB inv_root3_M;

    // Let 'rA' be a primitive root of order 3*2^k.
    // Then 'rootA_M' is M(rA), 'inv_rootA_M' is M(rA^-1).
    FFT_ULIMB rootA_M;
    FFT_ULIMB inv_rootA_M;

    // Value of 3^-1.
    FFT_ULIMB inv_3;

    // inv_2[i] = (2^i)^-1.
    //
    // Points to a static array of length k, so 0 <= i < k.
    const FFT_ULIMB *inv_2;
} Field;

#if FFT_LIMB_BITS == 64
#   include "consts_64.inc"
#elif FFT_LIMB_BITS == 32
#   include "consts_32.inc"
#else
#   error "Unknown value of FFT_LIMB_BITS."
#endif

static inline FFT_FORCE_INLINE
FFT_DOUBLE_ULIMB full_prod(FFT_ULIMB x, FFT_ULIMB y)
{
    return ((FFT_DOUBLE_ULIMB) x) * y;
}

// Performs an "adjustment" modulo field->p: the behaviour is equivalent to
//   if (x >= field->p) x -= field->p; return x;
//
// Assumes: x < 2*field->p.
//
// Returns: y such that y < field->p and y = x (mod field->p).
static inline FFT_FORCE_INLINE
FFT_ULIMB modp_adjust(FFT_ULIMB x, const Field *field)
{
    FFT_SLIMB d = x - field->p;
    FFT_SLIMB cf = d >> (FFT_LIMB_BITS - 1);
    return d + (field->p & cf);
}

// Performs addition modulo field->p.
//
// Assumes: x < field->p, y < field->p.
//
// Returns: z such that z < field->p and z = x + y (mod field->p).
static inline FFT_FORCE_INLINE
FFT_ULIMB modp_add(FFT_ULIMB x, FFT_ULIMB y, const Field *field)
{
    return modp_adjust(x + y, field);
}

static inline FFT_FORCE_INLINE
FFT_ULIMB modp_add3(FFT_ULIMB x, FFT_ULIMB y, FFT_ULIMB z, const Field *field)
{
    return modp_add(x, modp_add(y, z, field), field);
}

// Performs subtraction modulo field->p.
//
// Assumes: x < field->p, y < field->p.
//
// Returns: z such that z < field->p and z = x - y (mod field->p).
static inline FFT_FORCE_INLINE
FFT_ULIMB modp_sub(FFT_ULIMB x, FFT_ULIMB y, const Field *field)
{
    FFT_SLIMB d = x - y;
    FFT_SLIMB cf = d >> (FFT_LIMB_BITS - 1);
    return d + (field->p & cf);
}

// Performs Montgomery reduction on (x * y).
//
// Assumes: x * y < field->p * R, where:
//   * R = 2^FFT_LIMB_BITS;
//   * the products are in mathematical sense, not truncated to FFT_ULIMB.
//
// Since p < R, it works when both x and y are "adjusted" (less than field->p); but since 2*p < R,
// it can also be used when exactly one of {x, y} is in "unadjusted" form (less than 2*field->p),
// and the other one is "adjusted".
//
// Returns: z such that z < field->p and z = x * y * R^-1 (mod field->p).
static inline FFT_FORCE_INLINE
FFT_ULIMB redc(FFT_ULIMB x, FFT_ULIMB y, const Field *field)
{
    FFT_DOUBLE_ULIMB t = full_prod(x, y);

    FFT_ULIMB m = ((FFT_ULIMB) t) * field->mont_p1;

    FFT_DOUBLE_ULIMB mp = full_prod(m, field->p);

    FFT_ULIMB r = (t + mp) >> FFT_LIMB_BITS;

    return modp_adjust(r, field);
}

// Decimation-in-time butterfly.
static inline FFT_FORCE_INLINE
void butterfly_dit(FFT_ULIMB *x, FFT_ULIMB *y, FFT_ULIMB w, const Field *field)
{
    FFT_ULIMB u = *x;
    FFT_ULIMB v = redc(*y, w, field);
    *x = modp_add(u, v, field);
    *y = modp_sub(u, v, field);
}

// Decimation-in-time butterfly for index zero.
static inline FFT_FORCE_INLINE
void butterfly_dit_0(FFT_ULIMB *x, FFT_ULIMB *y, const Field *field)
{
    FFT_ULIMB u = *x;
    FFT_ULIMB v = *y;
    *x = modp_add(u, v, field);
    *y = modp_sub(u, v, field);
}

// Decimation-in-frequency butterfly.
static inline FFT_FORCE_INLINE
void butterfly_dif(FFT_ULIMB *x, FFT_ULIMB *y, FFT_ULIMB w, const Field *field)
{
    FFT_ULIMB u = *x;
    FFT_ULIMB v = *y;
    *x = modp_add(u, v, field);
    // 'u - v + field->p' is unadjusted, this is fine; see the comment for "redc".
    *y = redc(u - v + field->p, w, field);
}

// Decimation-in-frequency butterfly for index zero.
static inline FFT_FORCE_INLINE
void butterfly_dif_0(FFT_ULIMB *x, FFT_ULIMB *y, const Field *field)
{
    FFT_ULIMB u = *x;
    FFT_ULIMB v = *y;
    *x = modp_add(u, v, field);
    *y = modp_sub(u, v, field);
}

// Assumes:
//   * w < field->p;
//   * n is a power of two;
//   * n <= (1 << field->k).
//
// Returns: w^(n >> field->k).
static inline FFT_FORCE_INLINE
FFT_ULIMB rescale_root(FFT_ULIMB w, size_t n, const Field *field)
{
    for (FFT_ULIMB x = ((FFT_ULIMB) 1) << field->k; x > n; x /= 2)
        w = redc(w, w, field);
    return w;
}

// Decimation in time.
static inline FFT_FORCE_INLINE
void dit_full(
        FFT_ULIMB *w, size_t n,
        FFT_ULIMB *root_pows,
        const Field *field)
{
    FFT_ULIMB *w_end = w + n;

    for (size_t half_len = 1; half_len < n; half_len *= 2) {
        FFT_ULIMB *x = w;
        do {
            FFT_ULIMB *y = x + half_len;

            butterfly_dit_0(x, y, field);
            for (size_t i = 1; i < half_len; ++i)
                butterfly_dit(x + i, y + i, root_pows[i], field);

            x = y + half_len;
        } while (x != w_end);

        root_pows += half_len;
    }
}

// Decimation in frequency.
static inline FFT_FORCE_INLINE
void dif_full(
        FFT_ULIMB *w, size_t n,
        FFT_ULIMB *root_pows,
        const Field *field)
{
    FFT_ULIMB *w_end = w + n;
    root_pows += n - 1;

    for (size_t half_len = n / 2; half_len; half_len /= 2) {
        root_pows -= half_len;

        FFT_ULIMB *x = w;
        do {
            FFT_ULIMB *y = x + half_len;

            butterfly_dif_0(x, y, field);
            for (size_t i = 1; i < half_len; ++i)
                butterfly_dif(x + i, y + i, root_pows[i], field);

            x = y + half_len;
        } while (x != w_end);
    }
}

// Make a "normal" power table: write w^i to out[i], 0 <= i < n.
static inline FFT_FORCE_INLINE
void mk_pow_table(FFT_ULIMB w, FFT_ULIMB *out, size_t n, const Field *field)
{
    FFT_ULIMB x = field->mont_rmodp;
    for (size_t i = 0; i < n; ++i) {
        out[i] = x;
        x = redc(x, w, field);
    }
}

// Make a "special" power table.
// Assumes n is a power of two.
static inline FFT_FORCE_INLINE
void mk_special_pow_table(FFT_ULIMB w, FFT_ULIMB *out, size_t n, const Field *field)
{
    if (n == 1)
        return;

    size_t half_n = n / 2;
    FFT_ULIMB *normal_tab_start = out + half_n - 1;
    w = rescale_root(w, n, field);
    mk_pow_table(w, normal_tab_start, half_n, field);

    FFT_ULIMB *p = out;
    for (size_t step = half_n; step > 1; step /= 2)
        for (size_t i = 0; i < half_n; i += step)
            *p++ = normal_tab_start[i];
}

// Assumes n is a power of two.
static inline FFT_FORCE_INLINE
void convolve(
        FFT_ULIMB *a, FFT_ULIMB *b, size_t n,
        FFT_ULIMB *scratch,
        const Field *field)
{
    unsigned order = fft_zu_counttz(n);

    // make table of powers for root (in Montgomery representation)
    mk_special_pow_table(field->root_M, scratch, n, field);

#if ! FFT_USE_LINEARITY_TRICK()
    // convert a into Montgomery representation
    for (size_t i = 0; i < n; ++i)
        a[i] = redc(a[i], field->mont_r2modp, field);
#endif

    // forward transform on a
    dif_full(a, n, scratch, field);

    if (b != a) {
#if ! FFT_USE_LINEARITY_TRICK()
        // convert b into Montgomery representation
        for (size_t i = 0; i < n; ++i)
            b[i] = redc(b[i], field->mont_r2modp, field);
#endif
        // forward transform on b
        dif_full(b, n, scratch, field);
    }

    // pointwise multiply: a *= b
    for (size_t i = 0; i < n; ++i)
        a[i] = redc(a[i], b[i], field);

    // make table of powers for root^-1 (in Montgomery representation)
    mk_special_pow_table(field->inv_root_M, scratch, n, field);

    // inverse transform on a (without division by n)
    dit_full(a, n, scratch, field);

#if ! FFT_USE_LINEARITY_TRICK()
    // divide by n and convert out of Montgomery representation
    FFT_ULIMB factor = field->inv_2[order];
#else
    // divide by n and account for Montgomery stuff
    FFT_ULIMB r3 = redc(field->mont_r2modp, field->mont_r2modp, field);
    FFT_ULIMB factor = redc(r3, field->inv_2[order], field);
#endif
    for (size_t i = 0; i < n; ++i)
        a[i] = redc(a[i], factor, field);
}

static inline FFT_FORCE_INLINE
void fourstep_forward(
        FFT_ULIMB *w0, size_t c,
        FFT_ULIMB factor,
        FFT_ULIMB *root_pows,
        const Field *field)
{
    // 1. Multiply everything by 'factor'.
    // 2. Length 3 transform on columns.
    // 3. Multiply each matrix element with index i,j by root^(i*j).

    FFT_ULIMB *w1 = w0 + c;
    FFT_ULIMB *w2 = w1 + c;

    FFT_ULIMB aroot1 = rescale_root(field->rootA_M, c, field);
    FFT_ULIMB aroot2 = redc(aroot1, aroot1, field);

    FFT_ULIMB root1 = field->root3_M;
    FFT_ULIMB root2 = redc(root1, root1, field);

    FFT_ULIMB factor0 = factor;
    FFT_ULIMB factor1 = factor;
    FFT_ULIMB factor2 = factor;

    for (size_t i = 0; i < c; ++i) {
        FFT_ULIMB x0 = w0[i];
        FFT_ULIMB x1 = w1[i];
        FFT_ULIMB x2 = w2[i];

        // 'v0', 'v1', 'v2' are unadjusted, this is fine; see the comment for "redc".
        FFT_ULIMB v0 = x0 + modp_add(x1, x2, field);
        FFT_ULIMB v1 = x0 + modp_add(redc(x1, root1, field), redc(x2, root2, field), field);
        FFT_ULIMB v2 = x0 + modp_add(redc(x1, root2, field), redc(x2, root1, field), field);

        w0[i] = redc(v0, factor0, field);
        w1[i] = redc(v1, factor1, field);
        w2[i] = redc(v2, factor2, field);

        factor1 = redc(factor1, aroot1, field);
        factor2 = redc(factor2, aroot2, field);
    }

    // 4. Length C transform on rows.
    FFT_ULIMB *w_end = w0 + c * 3;
    for (FFT_ULIMB *w = w0; w != w_end; w += c)
        dif_full(w, c, root_pows, field);
}

static inline FFT_FORCE_INLINE
void fourstep_inverse(
        FFT_ULIMB *w0, size_t c,
        FFT_ULIMB factor,
        FFT_ULIMB *root_pows,
        const Field *field)
{
    // 1. Length C transform on rows.
    FFT_ULIMB *w_end = w0 + c * 3;
    for (FFT_ULIMB *w = w0; w != w_end; w += c)
        dit_full(w, c, root_pows, field);

    // 2. Multiply each matrix element with index i,j by root^(i*j).
    // 3. Length 3 transform on columns.
    // 4. Multiply everything by 'factor'.

    FFT_ULIMB *w1 = w0 + c;
    FFT_ULIMB *w2 = w1 + c;

    FFT_ULIMB aroot1 = rescale_root(field->inv_rootA_M, c, field);
    FFT_ULIMB aroot2 = redc(aroot1, aroot1, field);

    FFT_ULIMB root1 = field->inv_root3_M;
    FFT_ULIMB root2 = redc(root1, root1, field);

    FFT_ULIMB factor0 = factor;
    FFT_ULIMB factor1 = factor;
    FFT_ULIMB factor2 = factor;

    for (size_t i = 0; i < c; ++i) {
        FFT_ULIMB x0 = w0[i];
        FFT_ULIMB x1 = w1[i];
        FFT_ULIMB x2 = w2[i];

        x0 = redc(x0, factor0, field);
        x1 = redc(x1, factor1, field);
        x2 = redc(x2, factor2, field);

        w0[i] = modp_add3(x0, x1, x2, field);
        w1[i] = modp_add3(x0, redc(x1, root1, field), redc(x2, root2, field), field);
        w2[i] = modp_add3(x0, redc(x1, root2, field), redc(x2, root1, field), field);

        factor1 = redc(factor1, aroot1, field);
        factor2 = redc(factor2, aroot2, field);
    }
}

// Assumes n = 3*2^k.
static inline FFT_FORCE_INLINE
void convolve_fourstep(
        FFT_ULIMB *a, FFT_ULIMB *b, size_t n,
        FFT_ULIMB *scratch,
        const Field *field)
{
    size_t c = n / 3;
    unsigned order = fft_zu_counttz(c);

    // make table of powers for root (in Montgomery representation)
    mk_special_pow_table(field->root_M, scratch, c, field);

    // forward transform on a (also convert it into Montgomery representation: thus factor is
    // 'field->mont_r2modp')
    fourstep_forward(a, c, field->mont_r2modp, scratch, field);

    if (b != a) {
        // Forward transform on b (also convert it into Montgomery representation: thus factor is
        // 'field->mont_r2modp')
        fourstep_forward(b, c, field->mont_r2modp, scratch, field);
    }

    // pointwise multiply: a *= b
    for (size_t i = 0; i < n; ++i)
        a[i] = redc(a[i], b[i], field);

    // make table of powers for root^-1 (in Montgomery representation)
    mk_special_pow_table(field->inv_root_M, scratch, c, field);

    // inverse transform on a (also divide it by n and convert out of Montgomery representation)
    FFT_ULIMB factor = redc(
        field->mont_r2modp,
        redc(field->inv_2[order], field->inv_3, field),
        field);
    fourstep_inverse(a, c, factor, scratch, field);
}

// Multiply each matrix element with index i,j by factor * root^(i*j).
static inline FFT_FORCE_INLINE
void sixstep_mul_matrix(
        FFT_ULIMB *x, size_t r, size_t c,
        FFT_ULIMB factor,
        FFT_ULIMB root,
        const Field *field)
{
    for (size_t j = 0; j < c; ++j)
        x[j] = redc(x[j], factor, field);
    x += c;

    FFT_ULIMB cur_root = root;
    for (size_t i = 1; i < r; ++i) {
        FFT_ULIMB cur_factor = factor;
        for (size_t j = 0; j < c; ++j) {
            x[j] = redc(x[j], cur_factor, field);
            cur_factor = redc(cur_factor, cur_root, field);
        }
        cur_root = redc(cur_root, root, field);
        x += c;
    }
}

static inline FFT_FORCE_INLINE
void sixstep_forward(
        FFT_ULIMB *x, size_t r, size_t c,
        FFT_ULIMB factor,
        FFT_ULIMB *root_pows,
        FFT_ULIMB *transpose_scratch,
        const Field *field)
{
    size_t n = r * c;

    // 1. Transpose the matrix.
    fftu_transpose(x, r, c, transpose_scratch);

    // 2. Length R transform on rows.
    for (size_t i = 0; i < n; i += r) {
        fftu_permute(x + i, r);
        dit_full(x + i, r, root_pows, field);
    }

    // 3. Transpose the matrix.
    fftu_transpose(x, c, r, transpose_scratch);

    // 4. Multiply each matrix element with index i,j by factor * root^(i*j).
    FFT_ULIMB root = rescale_root(field->root_M, n, field);
    sixstep_mul_matrix(x, r, c, factor, root, field);

    // 5. Length C transform on rows.
    for (size_t i = 0; i < n; i += c)
        dif_full(x + i, c, root_pows, field);
}

static inline FFT_FORCE_INLINE
void sixstep_inverse(
        FFT_ULIMB *x, size_t r, size_t c,
        FFT_ULIMB factor,
        FFT_ULIMB *root_pows,
        FFT_ULIMB *transpose_scratch,
        const Field *field)
{
    size_t n = r * c;

    // 1. Length C transform on rows.
    for (size_t i = 0; i < n; i += c)
        dit_full(x + i, c, root_pows, field);

    // 2. Multiply each matrix element with index i,j by factor * root^(i*j).
    FFT_ULIMB root = rescale_root(field->inv_root_M, n, field);
    sixstep_mul_matrix(x, r, c, factor, root, field);

    // 3. Transpose the matrix.
    fftu_transpose(x, r, c, transpose_scratch);

    // 4. Length R transform on rows.
    for (size_t i = 0; i < n; i += r) {
        fftu_permute(x + i, r);
        dit_full(x + i, r, root_pows, field);
    }

    // 5. Transpose the matrix.
    fftu_transpose(x, c, r, transpose_scratch);
}

// Assumes n is a power of two.
static inline FFT_FORCE_INLINE
void convolve_sixstep(
        FFT_ULIMB *a, FFT_ULIMB *b, size_t n,
        FFT_ULIMB *scratch,
        const Field *field)
{
    unsigned order = fft_zu_counttz(n);

    unsigned r_order = order / 2;
    unsigned c_order = order - r_order;

    size_t r = ((size_t) 1) << r_order;
    size_t c = ((size_t) 1) << c_order;

    FFT_ULIMB *transpose_scratch = scratch + c - 1;

    // make table of powers for root (in Montgomery representation)
    mk_special_pow_table(field->root_M, scratch, c, field);

    // forward transform on a (also convert it into Montgomery representation: thus factor is
    // 'field->mont_r2modp')
    sixstep_forward(a, r, c, field->mont_r2modp, scratch, transpose_scratch, field);

    if (b != a) {
        // Forward transform on b (also convert it into Montgomery representation: thus factor is
        // 'field->mont_r2modp')
        sixstep_forward(b, r, c, field->mont_r2modp, scratch, transpose_scratch, field);
    }

    // pointwise multiply: a *= b
    for (size_t i = 0; i < n; ++i)
        a[i] = redc(a[i], b[i], field);

    // make table of powers for root^-1 (in Montgomery representation)
    mk_special_pow_table(field->inv_root_M, scratch, c, field);

    // inverse transform on a (also divide it by n and convert out of Montgomery representation)
    FFT_ULIMB factor = field->inv_2[order];
    sixstep_inverse(a, r, c, factor, scratch, transpose_scratch, field);
}

static inline FFT_FORCE_INLINE
FFT_DOUBLE_ULIMB crt2(FFT_ULIMB x1, FFT_ULIMB x2)
{
    FFT_ULIMB a1 = x1;

    FFT_ULIMB a2 = redc(
        FFT_CRT_R12_FACTOR,
        // 'x2 - a1 + FFT_field_2.p' is unadjusted, this is fine; see the comment for "redc".
        x2 - a1 + FFT_field_2.p,
        &FFT_field_2);

    return a1 + full_prod(a2, FFT_field_1.p);
}

static inline FFT_FORCE_INLINE
FFT_DOUBLE_ULIMB big_add3(
        FFT_DOUBLE_ULIMB a,
        FFT_DOUBLE_ULIMB b,
        FFT_DOUBLE_ULIMB c,
        FFT_DOUBLE_ULIMB *add_carry_to)
{
    FFT_DOUBLE_ULIMB s1 = a + b;
    FFT_DOUBLE_ULIMB s2 = s1 + c;
    *add_carry_to += (s1 < a) + (s2 < c);
    return s2;
}

static inline FFT_FORCE_INLINE
FFT_DOUBLE_ULIMB big_mulh(FFT_DOUBLE_ULIMB a, FFT_ULIMB b0, FFT_ULIMB b1)
{
    FFT_ULIMB a0 = (FFT_ULIMB) a;
    FFT_ULIMB a1 = a >> FFT_LIMB_BITS;

    FFT_DOUBLE_ULIMB x00 = full_prod(a0, b0);
    FFT_DOUBLE_ULIMB x01 = full_prod(a0, b1);
    FFT_DOUBLE_ULIMB x10 = full_prod(a1, b0);
    FFT_DOUBLE_ULIMB x11 = full_prod(a1, b1);

    FFT_DOUBLE_ULIMB r = x11;
    FFT_DOUBLE_ULIMB mid = big_add3(x01, x10, x00 >> FFT_LIMB_BITS, &r);
    return r + (mid >> FFT_LIMB_BITS);
}

void fft(
        FFT_ULIMB *a1, FFT_ULIMB *b1,
        FFT_ULIMB *a2, FFT_ULIMB *b2,
        FFT_ULIMB *scratch,
        size_t n)
{
    convolve(a1, b1, n, scratch, &FFT_field_1);
    convolve(a2, b2, n, scratch, &FFT_field_2);
}

void fft_fourstep(
        FFT_ULIMB *a1, FFT_ULIMB *b1,
        FFT_ULIMB *a2, FFT_ULIMB *b2,
        FFT_ULIMB *scratch,
        size_t n)
{
    convolve_fourstep(a1, b1, n, scratch, &FFT_field_1);
    convolve_fourstep(a2, b2, n, scratch, &FFT_field_2);
}

void fft_sixstep(
        FFT_ULIMB *a1, FFT_ULIMB *b1,
        FFT_ULIMB *a2, FFT_ULIMB *b2,
        FFT_ULIMB *scratch,
        size_t n)
{
    convolve_sixstep(a1, b1, n, scratch, &FFT_field_1);
    convolve_sixstep(a2, b2, n, scratch, &FFT_field_2);
}

#if FFT_LIMB_BITS == 64
#   include "recover_answer_64.inc"
#   include "thresholds_64.inc"
#elif FFT_LIMB_BITS == 32
#   include "recover_answer_32.inc"
#   include "thresholds_32.inc"
#else
#   error "Unknown value of FFT_LIMB_BITS."
#endif

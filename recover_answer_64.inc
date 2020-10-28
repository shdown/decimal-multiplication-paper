// Auto-generated; do not edit.

static void recover_answer_14(
    FFT_ULIMB *out, size_t nout, FFT_ULIMB *a1, FFT_ULIMB *a2)
{
 FFT_DOUBLE_ULIMB carry = 0;
 for (size_t i = 0; i < nout; ++i) {
  carry += crt2(a1[i], a2[i]);
  FFT_DOUBLE_ULIMB q = (big_mulh((carry), UINT64_C(6018080969204141205), UINT64_C(12980742146337069071))) >> 46;
  out[i] = carry - q * UINT64_C(100000000000000);
  carry = q;
 }
}

static void recover_answer_15(
    FFT_ULIMB *out, size_t nout, FFT_ULIMB *a1, FFT_ULIMB *a2)
{
 FFT_DOUBLE_ULIMB carry = 0;
 for (size_t i = 0; i < nout; ++i) {
  carry += crt2(a1[i], a2[i]);
  FFT_DOUBLE_ULIMB q = (big_mulh((carry) >> 15, UINT64_C(12926525502088622041), UINT64_C(633825300114114))) >> 20;
  out[i] = carry - q * UINT64_C(1000000000000000);
  carry = q;
 }
}

static void recover_answer_16(
    FFT_ULIMB *out, size_t nout, FFT_ULIMB *a1, FFT_ULIMB *a2)
{
 FFT_DOUBLE_ULIMB carry = 0;
 for (size_t i = 0; i < nout; ++i) {
  carry += crt2(a1[i], a2[i]);
  FFT_DOUBLE_ULIMB q = (big_mulh((carry), UINT64_C(15207441643216202349), UINT64_C(4153837486827862102))) >> 51;
  out[i] = carry - q * UINT64_C(10000000000000000);
  carry = q;
 }
}

static void recover_answer_17(
    FFT_ULIMB *out, size_t nout, FFT_ULIMB *a1, FFT_ULIMB *a2)
{
 FFT_DOUBLE_ULIMB carry = 0;
 for (size_t i = 0; i < nout; ++i) {
  carry += crt2(a1[i], a2[i]);
  FFT_DOUBLE_ULIMB q = (big_mulh((carry), UINT64_C(702278666647013315), UINT64_C(13292279957849158729))) >> 56;
  out[i] = carry - q * UINT64_C(100000000000000000);
  carry = q;
 }
}

void fft_recover_answer(
        FFT_ULIMB *out, size_t nout, FFT_ULIMB *a1, FFT_ULIMB *a2, int base_log)
{
 switch (base_log) {
 case 14: recover_answer_14(out, nout, a1, a2); break;
 case 15: recover_answer_15(out, nout, a1, a2); break;
 case 16: recover_answer_16(out, nout, a1, a2); break;
 case 17: recover_answer_17(out, nout, a1, a2); break;
 default: FFT_UNREACHABLE();
 }
}

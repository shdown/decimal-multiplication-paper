// Auto-generated; do not edit.

static void recover_answer_5(
    FFT_ULIMB *out, size_t nout, FFT_ULIMB *a1, FFT_ULIMB *a2)
{
 FFT_DOUBLE_ULIMB carry = 0;
 for (size_t i = 0; i < nout; ++i) {
  carry += crt2(a1[i], a2[i]);
  FFT_DOUBLE_ULIMB q = (big_mulh((carry) >> 5, UINT32_C(1907652675), UINT32_C(175921860))) >> 7;
  out[i] = carry - q * UINT32_C(100000);
  carry = q;
 }
}

static void recover_answer_6(
    FFT_ULIMB *out, size_t nout, FFT_ULIMB *a1, FFT_ULIMB *a2)
{
 FFT_DOUBLE_ULIMB carry = 0;
 for (size_t i = 0; i < nout; ++i) {
  carry += crt2(a1[i], a2[i]);
  FFT_DOUBLE_ULIMB q = (big_mulh((carry), UINT32_C(3619042523), UINT32_C(1125899906))) >> 18;
  out[i] = carry - q * UINT32_C(1000000);
  carry = q;
 }
}

static void recover_answer_7(
    FFT_ULIMB *out, size_t nout, FFT_ULIMB *a1, FFT_ULIMB *a2)
{
 FFT_DOUBLE_ULIMB carry = 0;
 for (size_t i = 0; i < nout; ++i) {
  carry += crt2(a1[i], a2[i]);
  FFT_DOUBLE_ULIMB q = (big_mulh((carry), UINT32_C(3849994941), UINT32_C(3602879701))) >> 23;
  out[i] = carry - q * UINT32_C(10000000);
  carry = q;
 }
}

void fft_recover_answer(
        FFT_ULIMB *out, size_t nout, FFT_ULIMB *a1, FFT_ULIMB *a2, int base_log)
{
 switch (base_log) {
 case 5: recover_answer_5(out, nout, a1, a2); break;
 case 6: recover_answer_6(out, nout, a1, a2); break;
 case 7: recover_answer_7(out, nout, a1, a2); break;
 default: FFT_UNREACHABLE();
 }
}

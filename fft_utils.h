// (c) 2020 shdown
// This code is licensed under MIT license (see LICENSE.MIT for details)

#pragma once

#include "fft.h"
#include <stddef.h>

// Performs bit-reversal permutation of size n, n=2^k.
void fftu_permute(FFT_ULIMB *x, size_t n);

// Transposes an r-by-c matrix in row-major order.
// Assumes either of:
//   * r = c = 2^k (in this case, scratch must have capacity of r>1 ? r*r/2 : 0);
//   * r = 2*c = 2^k (in this case, scratch must have capacity of max(c*c/2, 2*c));
//   * c = 2*r = 2^k (in this case, scratch must have capacity of max(r*r/2, 2*r)).
void fftu_transpose(FFT_ULIMB *x, size_t r, size_t c, FFT_ULIMB *scratch);

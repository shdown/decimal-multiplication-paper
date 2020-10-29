#!/usr/bin/env python3
import argparse

from isprime import is_prime
from utils import egcd, canonical_mod, prim_root, mod_pow, modp_inv
import fastdiv


class Platform:
    def __init__(self, bitwidth, k_min, p_min, base_logs):
        self.bitwidth = bitwidth
        self.k_min = k_min
        self.p_min = p_min
        self.base_logs = base_logs


PLATFORM_32 = Platform(
    bitwidth=32,
    k_min=25,
    p_min=int(2e9),
    base_logs=[5, 6, 7])


PLATFORM_64 = Platform(
    bitwidth=64,
    k_min=47,
    p_min=int(9.2e18),
    base_logs=[14, 15, 16, 17])


class Prime:
    def __init__(self, c, k):
        self.c = c
        self.k = k
        self.p = (c << k) + 1


def generate_primes(p_min, p_max, k_min, k_max):

    def _div_ceil(a, b):
        q, r = divmod(a, b)
        return q + (1 if r else 0)

    def _round_up_c(n, k):
        q = _div_ceil(_div_ceil(n, 1 << k), 3)
        if q % 2 == 0:
            q += 1
        return 3 * q

    def _round_down_c(n, k):
        q = (n >> k) // 3
        if q % 2 == 0:
            q -= 1
        return 3 * q

    for k in range(k_max, k_min - 1, -1):
        c_min = _round_up_c(p_min, k)
        c_max = _round_down_c(p_max, k)
        for c in range(c_max, c_min - 1, -6):
            p = (c << k) + 1
            assert p >= p_min
            assert p <= p_max
            if not is_prime(p):
                continue
            yield Prime(c, k)


def find_primes_for_platform(platform):
    g = generate_primes(
        p_min=platform.p_min,
        p_max=1 << (platform.bitwidth - 1),
        k_min=platform.k_min,
        k_max=platform.bitwidth - 2)
    prime1 = next(g)
    prime2 = next(g)
    return list(sorted([prime1, prime2], key=lambda prime: prime.p))


def print_consts_for_platform(platform, prime1, prime2):
    print('// Auto-generated; do not edit.')

    MONT_R = 1 << platform.bitwidth

    for number, prime in [(1, prime1), (2, prime2)]:
        neg_mont_p1, _ = egcd(prime.p, MONT_R)
        mont_p1 = canonical_mod(-neg_mont_p1, MONT_R)
        x = prim_root(prime.p)

        print(f'static const FFT_ULIMB FFT_ARRAY_{number}[] =')
        print('{')
        for i in range(prime.k):
            value = (1 << i) % prime.p
            print(f' {modp_inv(value, prime.p)},')
        print('};')

        print(f'static const Field FFT_field_{number} =')
        print('{')

        print(f' /*p=*/ {prime.p},')
        print(f' /*k=*/ {prime.k},')
        print(f' /*mont_p1=*/ {mont_p1},')
        print(f' /*mont_rmodp=*/ {MONT_R % prime.p},')
        print(f' /*mont_r2modp=*/ {MONT_R ** 2 % prime.p},')

        def to_mont(a):
            return a * MONT_R % prime.p

        def print_root_data(name, order):
            root = mod_pow(x, (prime.p - 1) // order, prime.p)
            inv_root = modp_inv(root, prime.p)
            print(f' /*{name}_M=*/ {to_mont(root)},')
            print(f' /*inv_{name}_M=*/ {to_mont(inv_root)},')

        print_root_data('root',  order=1 << prime.k)
        print_root_data('root3', order=3)
        print_root_data('rootA', order=3 << prime.k)

        print(f' /*inv_3=*/ {modp_inv(3, prime.p)},')
        print(f' /*inv_2=*/ FFT_ARRAY_{number},')
        print('};')

    crt_r12 = modp_inv(prime1.p, prime2.p) * MONT_R % prime2.p
    print(f'static const FFT_ULIMB FFT_CRT_R12_FACTOR = {crt_r12};')

    min_k = min(prime1.k, prime2.k)
    print('enum {')
    print(f' FFT_MIN_K = {min_k},')
    print('};')


def print_thresholds_for_platform(platform, prime1, prime2):
    print('// Auto-generated; do not edit.')

    product = prime1.p * prime2.p

    print('const FFT_Threshold FFT_THRESHOLDS[] =')
    print('{')
    print('#define NDIGITS(X_) ((X_) < SIZE_MAX ? (X_) : SIZE_MAX)')

    def make_const(n):
        if n <= 0xFFFFFFFF:
            return f'UINT32_C({n})'
        return f'UINT64_C({n})'

    for base_log in sorted(platform.base_logs, reverse=True):
        base = 10 ** base_log
        max_n = (product - 1) // ((base - 1) ** 2)
        ndigits = max_n * base_log * 2
        print(' {NDIGITS(%s), %s},' % (make_const(ndigits), base_log))

    print('#undef NDIGITS')

    print(' {0, 0},')
    print('};')


def print_fastdiv_for_platform(platform, prime1, prime2):
    print('// Auto-generated; do not edit.')

    def issue_srl(expr, nbits):
        if nbits:
            return f'({expr}) >> {nbits}'
        else:
            return expr

    def issue_const(x):
        return f'UINT{platform.bitwidth}_C({x})'

    def get_lo_hi(x):
        lo = x & ((1 << platform.bitwidth) - 1)
        hi = x >> platform.bitwidth
        return lo, hi

    for base_log in platform.base_logs:
        base = 10 ** base_log
        c = fastdiv.fastdiv(d=base, N=platform.bitwidth * 2)
        if not isinstance(c, fastdiv.CaseEasy):
            raise NotImplementedError()

        print(f'''
static void recover_answer_{base_log}(
    FFT_ULIMB *out, size_t nout, FFT_ULIMB *a1, FFT_ULIMB *a2)''')
        print('{')
        print(' FFT_DOUBLE_ULIMB carry = 0;')
        print(' for (size_t i = 0; i < nout; ++i) {')
        print('  carry += crt2(a1[i], a2[i]);')

        f0, f1 = get_lo_hi(c.factor)
        q = issue_srl('carry', c.sh_pre)
        q = f'big_mulh({q}, {issue_const(f0)}, {issue_const(f1)})'
        q = issue_srl(q, c.sh_post)

        print(f'  FFT_DOUBLE_ULIMB q = {q};')
        print(f'  out[i] = carry - q * {issue_const(base)};')
        print(f'  carry = q;')

        print(' }')
        print('}')

    print('''
void fft_recover_answer(
        FFT_ULIMB *out, size_t nout, FFT_ULIMB *a1, FFT_ULIMB *a2, int base_log)''')
    print('{')
    print(' switch (base_log) {')

    for base_log in platform.base_logs:
        print(f' case {base_log}: recover_answer_{base_log}(out, nout, a1, a2); break;')

    print(' default: FFT_UNREACHABLE();')

    print(' }')
    print('}')


def main():
    platforms = {
        '32': PLATFORM_32,
        '64': PLATFORM_64,
    }
    actions = {
        'print_fastdiv': print_fastdiv_for_platform,
        'print_consts': print_consts_for_platform,
        'print_thresholds': print_thresholds_for_platform,
    }

    parser = argparse.ArgumentParser()

    parser.add_argument(
        'platform',
        help='platform word size',
        type=lambda x: platforms[x])

    parser.add_argument(
        'action',
        help='action to perform',
        type=lambda x: actions[x])

    args = parser.parse_args()

    platform, action = args.platform, args.action
    p1, p2 = find_primes_for_platform(platform)
    action(platform, p1, p2)


if __name__ == '__main__':
    main()

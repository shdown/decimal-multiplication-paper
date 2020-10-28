def egcd(a, b):
    s = 0
    old_s = 1
    r = b
    old_r = a

    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s

    if b:
        bz_t = (old_r - old_s * a) // b
    else:
        bz_t = 0

    return old_s, bz_t


def canonical_mod(a, b):
    a %= b
    if a < 0:
        a += b
    return a


def factorize(n):
    result = []
    i = 2
    while i * i <= n:
        if n % i == 0:
            result.append(i)
            while n % i == 0:
                n //= i
        i += 1

    if n > 1:
        result.append(i)

    return result


def mod_pow(b, e, m):
    r = 1 % m
    while e:
        if e & 1:
            r = r * b % m
        b = b * b % m
        e >>= 1
    return r


def modp_inv(a, p):
    return mod_pow(a, p - 2, p)


def prim_root(p):
    phi = p - 1
    factors = factorize(phi)

    def check_candidate(g):
        for factor in factors:
            if mod_pow(g, phi // factor, p) == 1:
                return False
        return True

    for g in range(2, p):
        if check_candidate(g):
            return g

    raise ValueError(f'cannot find primitive root modulo {p}')


def decompose_pow2(n):
    assert n > 0
    e = 0
    while (n & 1) == 0:
        n >>= 1
        e += 1
    return e, n


def ceil_log2(n):
    assert n > 0
    r = 0
    n -= 1
    while n:
        n >>= 1
        r += 1
    return r

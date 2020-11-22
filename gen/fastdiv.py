# (c) 2020 shdown
# This code is licensed under MIT license (see LICENSE.MIT for details)

from utils import ceil_log2, decompose_pow2


# See: https://gmplib.org/~tege/divcnst-pldi94.pdf


class CasePowerOfTwo:
    """
    q = SRL(n, <self.ell>)
    """
    def __init__(self, ell):
        self.ell = ell


class CaseEasy:
    """
    q = SRL(MULUH(<self.factor>, SRL(n, <self.sh_pre>)), <self.sh_post>)
    """
    def __init__(self, factor, sh_pre, sh_post):
        self.factor = factor
        self.sh_pre = sh_pre
        self.sh_post = sh_post


class CaseHard:
    """
    t1 = MULUH(<self.factor>, n)
    q = SRL(t1 + SRL(n - t1, 1), <self.sh>)
    """
    def __init__(self, factor, sh):
        self.factor = factor
        self.sh = sh


def choose_multiplier(d, prec, N):
    ell = ceil_log2(d)
    sh_post = ell

    m_lo = (1 << (N + ell)) // d

    m_hi = ((1 << (N + ell)) + (1 << (N + ell - prec))) // d

    while (m_lo // 2) < (m_hi // 2) and sh_post:
        m_lo //= 2
        m_hi //= 2
        sh_post -= 1

    return m_hi, sh_post, ell


# N is the bitwidth of numbers that we are going to divide by 'd', e.g. 32 or 64.
def fastdiv(d, N):
    assert d > 0

    m, sh_post, ell = choose_multiplier(d, N, N)
    if m >= (1 << N) and d % 2 == 0:
        e, d_odd = decompose_pow2(d)
        sh_pre = e
        m, sh_post, _ = choose_multiplier(d_odd, N - e, N)
    else:
        sh_pre = 0

    if d == (1 << ell):
        return CasePowerOfTwo(ell=ell)
    elif m >= (1 << N):
        assert sh_pre == 0
        return CaseHard(factor=m - (1 << N), sh=sh_post - 1)
    else:
        return CaseEasy(factor=m, sh_pre=sh_pre, sh_post=sh_post)

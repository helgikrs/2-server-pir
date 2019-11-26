#!/usr/bin/env python3
"""Faster, pure Python MVF generator."""

import functools
import math

from tqdm import tqdm


def nCr(n, r):
    if r > n:
        return 0
    return math.factorial(n)//math.factorial(r)//math.factorial(n-r)


def prod(l):
    r = 1
    for p in l:
        r *= p
    return r


def crt(res):
    # XXX: hardcoded CRT over Z6
    m2 = res[0] % 2
    m3 = res[1] % 3
    if (m2, m3) == (0, 0):
        return 0
    if (m2, m3) == (1, 0):
        return 3
    if (m2, m3) == (0, 1):
        return 4
    if (m2, m3) == (1, 1):
        return 1
    if (m2, m3) == (0, 2):
        return 2
    if (m2, m3) == (1, 2):
        return 5


def fp(x_ham, p, r):
    return 1 - prod(1 - nCr(x_ham, p**i)**(p-1) for i in range(r))


@functools.lru_cache(maxsize=None)
def f(x_ham, x_len, ps):
    n = x_len
    t = len(ps)
    m = prod(ps)

    res = []

    nt = n**(1/t)

    for p in ps:
        r = 0
        while p**r <= nt:
            r += 1
        res.append(fp(x_ham, p, r))

    return crt(res)


@functools.lru_cache(maxsize=None)
def ham(x):
    return bin(x).count('1')


def make_mvf_M(n):
    # writes directly to file instead of keeping in-memory
    with open('mvf-large.dat', 'w') as file:
        file.write('[')
        for y in tqdm(range(2**n)):
            if y > 0:
                file.write(',\n')
            file.write('[')
            for x in range(2**n):
                if x > 0:
                    file.write(', ')
                mxy = f(ham(x ^ y), n, (2, 3))
                file.write(str(mxy))
            file.write(']')
        file.write(']')

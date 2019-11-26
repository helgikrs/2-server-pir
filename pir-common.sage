# -*- encoding=utf8
from collections import namedtuple

load("mvf.sage")

MVF = namedtuple('MVF', ['v', 'u'])

def gammify(v):
    return str(v).replace('gamma', '𝛾')

def R(m, r):
    """Returns the ring of polynomials of degree < r over integers mod m."""
    Zm = Integers(m)
    R.<gamma> = PolynomialRing(Zm)
    Rmr.<gamma> = R.quotient(gamma ^ r - 1)
    return Rmr, gamma

def random_vector(group, k):
    return vector([group.random_element() for _ in range(k)])

def vecexp(gamma, vec):
    return vector([gamma ^ i for i in vec])

def vecvecexp(xs, zs):
    k = len(xs)
    assert len(zs) == k
    return prod([xs[i] ^ zs[i] for i in range(k)])

def embed(PR, a, mvf):
    n = len(a)
    k = len(mvf.u[0])
    xs = PR.gens()
    assert len(xs) == k
    assert len(mvf.u) == n

    r = [a[i] * vecvecexp(xs, mvf.u[i]) for i in range(n)]
    d = [a[i] * mvf.u[i] for i in range(n)]

    return sum(r), lambda *xs: sum(d[i] * r[i](*xs) for i in range(len(r)))

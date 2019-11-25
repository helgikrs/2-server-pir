# -*- encoding=utf8
from collections import namedtuple
import sys
import time

load("mvf.sage")

MVF = namedtuple('MVF', ['v', 'u'])

def gammify(v):
    return str(v).replace('gamma', '𝛾')


def R(m, r):
    Zm = Integers(m)
    R.<gamma> = PolynomialRing(Zm)
    R66.<gamma> = R.quotient(gamma ^ r - 1)
    return R66, gamma

def random_vector(group, k):
    return vector([group.random_element() for _ in range(k)])


def is_smvf(s, v, u):
    n = len(v)
    assert len(u) == n

    for i in range(n):
        for j in range(n):
            dot = u[i].dot_product(v[j])
            if i == j and dot != 0:
                return False
            if i != j and dot not in s:
                return False
    return True


def MVF_gen(S, group, k, n):
    while True:
        print('Generating a random family')
        u = [random_vector(group, k) for _ in range(n)]
        v = [random_vector(group, k) for _ in range(n)]

        if is_smvf(S, u, v):
            print(u)
            print(v)
            break


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


def check(db, w):
    """
    Embed db as polynomial, query point i and assert that the results are
    correct.
    """

    n = ceil(log(len(db), 2))

    print("Computing MVF(n={})".format(n))
    u, v = make_mvf(n)
    # print("U = {}".format(u))
    # print("V = {}".format(v))

    # print('MVF degree:', len(u[0]))
    # print('MVF size:', len(u))

    mvf = MVF(u=u, v=v)

    n = len(u)
    k = len(u[0])

    Z = Integers(6)
    R66, gamma = R(6, 6)
    PR = PolynomialRing(R66, k, 'x')
    # PR = PolynomialRing(Z, k, 'x')
    # gamma = Z(-1)

    print("Computing db polynomial F")
    edb, deriv = embed(PR, db, mvf)
    print("F{} = {}".format(PR.gens(), edb))

    print("==============")

    print("User querying letter w = {}".format(w))

    letter = []
    for j in range(8):
        i = j + w * 8

        # protocol, for quering point i
        z = random_vector(Z, k)
        q1 = vecexp(gamma, z)
        q2 = vecexp(gamma, z + mvf.v[i])

        # U sends q_i to server S_i
        # S_i responds with r_i and d_i
        r1, d1 = edb(*q1), deriv(*q1)
        r2, d2 = edb(*q2), deriv(*q2)

        # User recovers bit of a_i
        gp1 = d1.dot_product(mvf.v[i])
        gp2 = d2.dot_product(mvf.v[i])

        M = Matrix([
            [1, 1, 1, 1],
            [0, 1, 3, 4],
            [1, gamma, gamma^3, gamma^4],
            [0, gamma, 3*gamma^3, 4*gamma^4]
        ])

        resv = (M.adjoint() * vector([r1, gp1, r2, gp2]))

        res = int(resv[0] != 0)

        print('Got {}{} bit, {}'.format(j, (['th', 'st', 'nd', 'rd'] + ['th'] * 10)[j], res))

        letter.append(res)

        assert res == (db[i] == 1)

    print("Received letter a_w = '{}'".format(chr(int(''.join(map(str, letter)), 2))))


def main():
    print("Started")
    dbt = "Hello, World!"
    print('db = "{}"'.format(dbt))

    n = N(2 ^ ceil(log(len(dbt) * 8, 2)))
    w = randint(0, len(dbt))

    db = []
    for i in dbt:
        for j in bin(ord(i))[2:].rjust(8, '0'):
            db.append(int(j))

    while len(db) < n:
        db.append(0)


    print('binary db = {}'.format(db))

    time.sleep(3)

    check(db, w)

    print("Actual a_w = '{}'".format(dbt[w]))


if __name__ == '__main__':
    main()

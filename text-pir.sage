# -*- encoding=utf8
import time

load("pir-common.sage")


def check(db, w):
    n = ceil(log(len(db), 2))

    print("Computing MVF(n={})".format(n))
    u, v = make_mvf(n)
    mvf = MVF(u=u, v=v)

    n = len(u)
    k = len(u[0])

    Z = Integers(6)
    R66, gamma = R(6, 6)
    PR = PolynomialRing(R66, k, 'x')

    print("Computing db polynomial F")
    edb, deriv = embed(PR, db, mvf)
    print("F{} = {}".format(PR.gens(), edb))

    print("==============")

    print("User querying letter w = {}".format(w))

    letter = []
    for j in range(8):
        i = j + w * 8

        z = random_vector(Z, k)
        q1 = vecexp(gamma, z)
        q2 = vecexp(gamma, z + mvf.v[i])

        r1, d1 = edb(*q1), deriv(*q1)
        r2, d2 = edb(*q2), deriv(*q2)

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

# -*- encoding=utf8

load("pir-common.sage")


def check(db, i):
    """
    Embed db as polynomial, query point i and assert that the results are
    correct.
    """
    print("db = {}\n".format(db))

    n = ceil(log(len(db), 2))

    print("Computing MVF(n={})".format(n))
    u, v = make_mvf(n)
    print("U = {}".format(u))
    print("V = {}".format(v))

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

    print("User querying point w = {}".format(i))

    # protocol, for quering point i
    z = random_vector(Z, k)
    print("User generates random vector z={}".format(z))

    # U sends q_i to server S_i
    q1 = vecexp(gamma, z)
    q2 = vecexp(gamma, z + mvf.v[i])
    print("User constructs two queries to send to servers S1 and S2")
    print(" q1 = {}".format(str(q1).replace('gamma', '𝛾')))
    print(" q2 = {}".format(str(q2).replace('gamma', '𝛾')))


    # S_i responds with r_i and d_i
    r1, d1 = edb(*q1), deriv(*q1)
    print("< Server 1 received q1")
    print("< Server 1 responds with [{}, {}]".format(gammify(r1), gammify(d1)))

    r2, d2 = edb(*q2), deriv(*q2)
    print("< Server 2 received q2")
    print("< Server 2 responds with [{}, {}]".format(gammify(r2), gammify(d2)))

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

    print('User computes:')

    print("adj(M) * {}\n\t= {}".format(gammify(vector([r1, gp1, r2, gp2])),
                                       gammify(resv)))

    res = resv[0] != 0

    print("User recovery of a_w = {}".format(int(res)))
    print("db[w] = {}".format(db[i]))

    assert res == (db[i] == 1)

def main():
    n = 4
    db = random_vector(GF(2), n).change_ring(Integers())
    w = randint(0, n-1)
    check(db, w)

if __name__ == '__main__':
    main()

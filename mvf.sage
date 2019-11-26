import os


def fp(xs, p, r):
    return 1 - prod(1 - Combinations(range(sum(xs)), p^i).cardinality()^(p-1) for i in range(r))

def f(xs, ps):
    n = len(xs)
    t = len(ps)
    m = prod(ps)

    res = []

    nt = n^(1/t)

    for p in ps:
        r = 0
        while p^r <= nt:
            r += 1
        res.append(fp(xs, p, r))

    return crt(res, ps)


def make_mvf(n):
    datname = '{}-mvf.dat'.format(n)

    if os.path.isfile(datname):
        with open(datname) as fd:
            return loads(fd.read())

    Zm = Integers(6)

    M = Matrix([[ f((vector(GF(2), x) + vector(GF(2), y)).change_ring(Integers()), [2, 3]) for x in Words([0, 1], n) ] for y in Words([0, 1], n) ], ring=Zm)

    A = M.change_ring(QQ)
    B = A.echelon_form()
    F = B[:B.rank()]

    A.pivots()
    C = Matrix([ A.column(j) for j in A.pivots() ]).T
    C = C.change_ring(Zm)
    F = F.change_ring(Zm)

    U, V = list(C), list(F.T)

    with open(datname, 'w') as fd:
        fd.write(dumps((U, V)))

    return U, V

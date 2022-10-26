#########################################
########## toy parameters ###############
# modulus
q = 8380417
# dimension of the polynomial
n = 16

# norm bound for the noise
beta = 1
# number of polynomials in a vector
l = 4

# number of rows
k = 3

#########################################
########## setup ########################
Rcal = Zmod(q)
P.<x> = PolynomialRing(Rcal)
F = P(x^n+1)
PQ.<x> = PolynomialQuotientRing(P, F)
#########################################
# scheme
# following the definition of `streaming` paper
# moving to a polynomial ring for efficiency


def parameters():
    '''
    generating the parameters L and R for the hash function
    note that we are using ring-SIS instead of SIS, i.e
        - L and R are k by l matrix, each entry is a polynomial;
        - x and y are l dimensional vectors of polynomials
    '''
    L = []
    R = []
    for _ in range (k):
        L_row = []
        R_row = []
        for _ in range (l):
            f = PQ([Rcal.random_element() for _ in range (n)])
            L_row.append(f)
            f = PQ([Rcal.random_element() for _ in range (n)])
            R_row.append(f)
        L.append(L_row)
        R.append(R_row)
    return (L, R)


def hash(x, y, L, R):
    '''
    digest = L x  + R y
    '''
    res = [0 for _ in range(l)]

    for i in range (k):
        for j in range (l):
            res[i] += (L[i][j] * x[j] + R[i][j] *  y[j])

    return res


#########################################
# testing

(L, R) = parameters()
x1 = [ PQ([ZZ.random_element(0, beta+1) for _ in range (n)]) for _ in range (l)]
y1 = [ PQ([ZZ.random_element(0, beta+1) for _ in range (n)]) for _ in range (l)]
digest1 = hash(x1, y1, L, R)


x2 = [ PQ([ZZ.random_element(0, beta+1) for _ in range (n)]) for _ in range (l)]
y2 = [ PQ([ZZ.random_element(0, beta+1) for _ in range (n)]) for _ in range (l)]
digest2 = hash(x2, y2, L, R)

# check that the hash function is additive homomorphic
x12 = [ x1[i] + x2[i] for i in range(l) ]
y12 = [ y1[i] + y2[i] for i in range(l) ]
digest12 = hash(x12, y12, L, R)

digest12_hom = [ digest1[i] + digest2[i] for i in range(l) ]


assert(digest12 == digest12_hom)

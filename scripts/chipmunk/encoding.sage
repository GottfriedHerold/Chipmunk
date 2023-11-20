# N: dimension
n = 16
# q: modulus
q = 202753
# q + 1 over 2
q_over_2 = 101377
# eta: log base
eta = 29
# kappa is log(q, 2eta+1) = 3
kappa = 3
# polynomial ring over Z
P.<x> = PolynomialRing(ZZ)

def poly_mul(a, b):
    return vector((P(a.list()) * P(b.list()) % P(x^n+1)).list())

# random polynomial
def random_poly(n, q):
    return vector([ZZ.random_element(q)- q_over_2 for _ in range(n)])

# random binary polynomial
def random_binary_poly(n, q):
    return vector([ZZ.random_element(2) for _ in range(n)])

# lift to [0, q)
def lift(poly):
    for e in poly:
        e = e % q
        if e < 0:
            e = e + q
    return poly

def center_mod(i, mod):
    i = i % mod
    if i > mod/2:
        i = i - mod
    return i

# center lift polynomial to [-mod/2, mod/2)
def center_lift(poly, mod):
    res = deepcopy(poly)
    for i in range(0, n):
        res[i] = center_mod(res[i], mod)
    return res

# project a vector of R elements into an R_q element
def proj_q(v):
    return center_lift(proj_eta_k(v), q)

# decompose an R_q element into a vector of R elements
def dec_q(poly):
    v1 = center_lift(poly, (2*eta+1))
    poly = (poly - v1)/(2*eta+1)
    v2 = center_lift(poly, (2*eta+1))
    v3 = (poly - v2)/(2*eta+1)
    return [v1,v2,v3]

# project a vector of R elements into an R element
def proj_eta_k(v):
    return v[0] + (2*eta+1)*v[1] + (2*eta+1)^2*v[2]

# decompose an R element into a vector of R elements
def dec_eta_k(v):
    poly = center_lift(v, q)
    print("p:", poly)
    return dec_q(poly)

# encode a vector of small elements
# figure 4 from paper
def encode(v):
    hint = proj_q(v)
    print("hint:   ", hint)
    hint_prime = lift(hint)
    print("hint':  ", hint_prime)
    proj_eta_k_v = proj_eta_k(v)
    dec_eta_k_v = dec_eta_k(proj_eta_k(v))
    alpha_star = (proj_eta_k_v - hint_prime)/q
    print("alpha*: ", alpha_star)
    print("proj_eta_k_v:    ", proj_eta_k_v)
    print("dec_eta_k_v:    ", dec_eta_k_v)
    print("v:      ", v)
    delta_v = [v[0] - dec_eta_k_v[0], v[1] - dec_eta_k_v[1], v[2] - dec_eta_k_v[2]]

    alpha_1 = -(v[0] - dec_eta_k_v[0])
    alpha_2 = alpha_1 - (v[1] - dec_eta_k_v[1])
    alpha_3 = dec_eta_k_v[2] - v[2]
    print("alpha 1:", alpha_1)
    print("alpha 2:", alpha_2)
    print("alpha 3:", alpha_3)

    print("detla_v ", delta_v)

    return delta_v

# find alpha s.t.
#  alpha1 (-2eta + 1, 1, 0) + alpha2 (0, -2eta + 1, 1) = delta_v
# we build a lattice spanned by the row vector of the following matrix
#
#  |rot(delta_v[0]) | rot(delta_v[1]) | rot(delta_v[2]) |   0    |   0    |
#  |dignal(-2eta+1) |     identity    |        0        |identity|   0    |
#  |       0        | dignal(-2eta+1) |     identity    |   0    |identity|
# ---------------------------------------------------------------------------
#  |     xxx        |       xxx       |       xxx       | alpha1 | alpha2 |
def find_alpha(delta_v):
    # m = matrix(n*3, n*5)
    # v0 = anti_cyclic_rotate(delta_v[0])
    # v1 = anti_cyclic_rotate(delta_v[1])
    # v2 = anti_cyclic_rotate(delta_v[2])
    # for i in range(n):
    #     for j in range(n):
    #         m[i,j] = v0[i,j]
    #         m[i,j+n] = v1[i,j]
    #         m[i,j+n*2] = v2[i,j]
    # for i in range(n):
    #     m[n+i,i] = -(eta*2+1)
    #     m[n+i,n+i] = 1
    #     m[n+i,3*n+i] = 1
    #     m[2*n+i,n+i] = -(eta*2+1)
    #     m[2*n+i,2*n+i] = 1
    #     m[2*n+i,4*n+i] = 1
    # print("before reducing:", m)
    #
    # return m
    return
#
# def anti_cyclic_rotate(f):
#     m = matrix(n,n)
#     for i in range(n):
#         for j in range(n):
#             if i+j >=n:
#                 m[i, (i+j)%n] = -f[j]
#             else:
#                 m[i, i+j] = f[j]
#     return m

# tests
p = random_poly(n, q)
print("p:  ", p)
# print("center lift: ", center_lift(p, 2*eta+1))

v = dec_q(p)
print("v0: ", v[0])
print("v1: ", v[1])
print("v2: ", v[2])

p2 = proj_q(v)
print("p2: ", p2)
print("p==p2: ", p == p2)

r = random_binary_poly(n, q)
vr = [poly_mul(v[0], r), poly_mul(v[1], r),poly_mul(v[2], r)]
vrp = proj_q(vr)
vrp2 = poly_mul(p, r)
print("vrp: ", vrp)
print("vrp2:", vrp2)
print("vrp - vrp2:", vector(vrp)-vector(vrp2))
print("vrp - vrp2:", (vector(vrp)-vector(vrp2))%q)

delta_v = encode(vr)
# m = find_alpha(delta_v)
# print("after reducing:", m.LLL())

from math import factorial as fac
from math import log, ceil, erf, sqrt

R = RealField(20)

def prune(d):
    r = {}
    for key in d:
        if d[key] >= 2^-128:
            r[key] = R(d[key])
    return r


def gaussian_center_weight(sigma, t):
    """ Weight of the gaussian of std deviation s, on the interval [-t, t]
    :param x: (float)
    :param y: (float)
    :returns: erf( t / (sigma*\\sqrt 2) )
    """
    return erf(t / (sigma * sqrt(2.)))


def binomial(x, y):
    """ Binomial coefficient
    :param x: (integer)
    :param y: (integer)
    :returns: y choose x
    """
    try:
        binom = fac(x) // fac(y) // fac(x - y)
    except ValueError:
        binom = 0
    return binom


def centered_binomial_pdf(k, x):
    """ Probability density function of the centered binomial law of param k at x
    :param k: (integer)
    :param x: (integer)
    :returns: p_k(x)
    """
    return binomial(2*k, x+k) / 2.**(2*k)


def build_centered_binomial_law(k):
    """ Construct the binomial law as a dictionnary
    :param k: (integer)
    :param x: (integer)
    :returns: A dictionnary {x:p_k(x) for x in {-k..k}}
    """
    D = {}
    for i in range(-k, k+1):
        D[i] = centered_binomial_pdf(k, i)
    return D


def mod_switch(x, q, rq):
    """ Modulus switching (rounding to a different discretization of the Torus)
    :param x: value to round (integer)
    :param q: input modulus (integer)
    :param rq: output modulus (integer)
    """
    return int(round(1.* rq * x / q) % rq)


def mod_centered(x, q):
    """ reduction mod q, centered (ie represented in -q/2 .. q/2)
    :param x: value to round (integer)
    :param q: input modulus (integer)
    """
    a = x % q
    if a < q/2:
        return a
    return a - q


def build_mod_switching_error_law(q, rq):
    """ Construct Error law: law of the difference introduced by switching from and back a uniform value mod q
    :param q: original modulus (integer)
    :param rq: intermediate modulus (integer)
    """
    D = {}
    V = {}
    for x in range(q):
        y = mod_switch(x, q, rq)
        z = mod_switch(y, rq, q)
        d = mod_centered(x - z, q)
        D[d] = D.get(d, 0) + 1./q
        V[y] = V.get(y, 0) + 1

    return D


def law_convolution(A, B):
    """ Construct the convolution of two laws (sum of independent variables from two input laws)
    :param A: first input law (dictionnary)
    :param B: second input law (dictionnary)
    """

    C = {}
    for a in A:
        for b in B:
            c = a+b
            C[c] = C.get(c, 0) + A[a] * B[b]
    return prune(C)


def law_product(A, B):
    """ Construct the law of the product of independent variables from two input laws
    :param A: first input law (dictionnary)
    :param B: second input law (dictionnary)
    """
    C = {}
    for a in A:
        for b in B:
            c = a*b
            C[c] = C.get(c, 0) + A[a] * B[b]
    return C


def clean_dist(A):
    """ Clean a distribution to accelerate further computation (drop element of the support with proba less than 2^-300)
    :param A: input law (dictionnary)
    """
    B = {}
    for (x, y) in A.items():
        if y>2**(-300):
            B[x] = y
    return B


def iter_law_convolution(A, i):
    """ compute the -ith forld convolution of a distribution (using double-and-add)
    :param A: first input law (dictionnary)
    :param i: (integer)
    """
    D = {0: 1.0}
    i_bin = bin(i)[2:]  # binary representation of n
    for ch in i_bin:
        D = law_convolution(D, D)
        D = clean_dist(D)
        if ch == '1':
            D = law_convolution(D, A)
            D = clean_dist(D)
    return D


def tail_probability(D, t):
    '''
    Probability that an drawn from D is strictly greater than t in absolute value
    :param D: Law (Dictionnary)
    :param t: tail parameter (integer)
    '''
    s = 0
    ma = max(D.keys())
    if t >= ma:
        return 0
    for i in reversed(range(int(ceil(t)), ma)):  # Summing in reverse for better numerical precision (assuming tails are decreasing)
        s += D.get(i, 0) + D.get(-i, 0)
    return s

def tail_param(D, l):
    '''
    Find the tail parameters t that a drawn from D is strictly greater than t with probability > 1-2^l
    :param D: Law (Dictionnary)
    :param l: log 2 of probability
    '''

    s = 0
    ma = max(D.keys())
    for i in reversed(range(ma)):  # Summing in reverse for better numerical precision (assuming tails are decreasing)
        if s > 2^-l:
            return (i, log(s,2))
        s += D.get(i, 0) + D.get(-i, 0)
    return (i, log(s,2))

#
# chi = {-5: 1/11, -4: 1/11, -3: 1/11, -2: 1/11, -1: 1/11, 0: 1/11, 1: 1/11, 2: 1/11, 3: 1/11, 4: 1/11, 5: 1/11}
# chi2 = iter_law_convolution(chi, 2)
# chi44 = iter_law_convolution(chi, 44)
# chi60 = iter_law_convolution(chi, 60)
#
#
# gamma_1 = 523776
# space = 2*gamma_1 + 1
# rho = RR(1/space)
# gamma_dist = {0: rho}
# for i in range(1, gamma_1+1):
#     gamma_dist[i] = rho
#     gamma_dist[-i]= rho
#
#
# # for n = 256, we have h = 60, the space is large enough
# # assert choose(256, 60) * 2^60 > 2^256
#
# # the distribution of c*s is an accumulation of 60 [-5, 5]
# chi60 = iter_law_convolution(chi, 60)
#
# # the distribution of y + c*s is
# final_dist = law_convolution(gamma_dist, chi60)
# print "final dist computed"
#
# beta = 275
# # unform between [-gamma, gamma]
# space = 2*(gamma_1) + 1
# rho = RR(1/space)
# uni_dist = {0: rho}
# for i in range (1, (gamma_1-beta)+1):
#     uni_dist[i] = rho
#     uni_dist[-i] = rho
#
# # now compute the SD between uni_dist and final_dist
# distance = (uni_dist[0] - final_dist[0])^2
# for i in range (1, (gamma_1-beta)+1):
#     distance += (uni_dist[i] - final_dist[i])^2
#     distance += (uni_dist[-i] - final_dist[-i])^2
#
# print "distance", sqrt(distance)
#
#
#
# max44 = 44*5
# max60 = 60*5
#
# sum44 = 0
# counter = 0
# while sum44 < 2^-128:
#     sum44 += chi44.get(max44-counter)*2
#     counter += 1
# print max44-counter, log(sum44,2)
# # 213 -127.293311878
#
#
# sum60 = 0
# counter = 0
# while sum60 < 2^-128:
#     sum60 += chi60.get(max60-counter)*2
#     counter += 1
# print max60-counter, log(sum60,2)
# # 269 -127.410773065

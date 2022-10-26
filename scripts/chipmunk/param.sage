load("dist.sage")
R = RealField(20)

# search for q
def search_q(beta, n):
    q = next_prime(beta * 4)
    while q % (2*n) != 1:
        q = next_prime(q)
    return q

def search_trinary_hamming_weight(n):
    for i in range(n):
        order = factorial(n)/factorial(n-i)/factorial(i/2)/factorial(i/2)
        if order > 2^128:
            return i

# dimension of the polynomial
n = 1024

# we set q to be an 14 bits prime integer
hvc_q = 12289 # q = 2^13 + 2^12 + 1, same q used by Falcon
hvc_log_q = 14


# base-alpha decomposition: an F_q element is decomposed into
# width number integers with [-base/2, base/2]
base = 9
width = ZZ(round(log(hvc_q, base)))

# the # non-zero entries in a randomizer polynomial
alpha = search_trinary_hamming_weight(n)

# each non-aggregated polynomial, denoted by `a` has coefficients between
# [0, base); after applying randomizer, denoted by `r`, has alpha/2 number of
# 1s and alpha/2 number of -1s, the result, denoted by `b := a r` has
# coefficients which are the accumulation of alpha numbers between -base/2 and
# base/2.

chi = {-4: 1/base, -3: 1/base, -2: 1/base, -1: 1/base, 0: 1/base, 1: 1/base, 2: 1/base, 3: 1/base, 4: 1/base}
chi_alpha = iter_law_convolution(chi, alpha)

# now we try to get an l2 norm of b whose coefficients follow chi_alpha
# first we get the quadratic distribution of chi_alpha
chi_alpha_square = {}
for key, value in chi_alpha.items():
    if key > 0:
        chi_alpha_square[key^2] = R(value)*2
    if key == 0:
        chi_alpha_square[key] = R(value)
        
# then we accumulate chi_alpha_square for n times to get the norm
# this part of the code runs forever...
l2_norm_dist = iter_law_convolution(chi_alpha_square, n)

for i in range(100, 100000):
    print(i, tail_probability(l2_norm_dist, i))

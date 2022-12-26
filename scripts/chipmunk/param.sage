reset()

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
        order = factorial(n)/factorial(n-i*2)/factorial(i)/factorial(i)
        if order > 2^128:
            return i*2

# =================================
# system setting
## number of users
rho = 1024
## height of the tree
tau = 21

# ring setting
## dimension of the polynomial
n = 512

## we set q to be an 14 bits prime integer
hvc_q = 12289 # q = 2^13 + 2^12 + 1, same q used by Falcon
hvc_log_q = 14

## base-alpha decomposition: an F_q element is decomposed into
## width number integers with [-base/2, base/2]
base = 9
base_over_2 = (base -1)/2
width = ZZ(ceil(log(hvc_q, base)))

# =================================
# computing the rest of the parameters

# the # non-zero entries in a randomizer polynomial
# we assume equal number of -1s and 1s
alpha = search_trinary_hamming_weight(n)

# now we want to expect the norm bound for the a * r where
# - each non-aggregated polynomial, denoted by `a` has coefficients between
# [-base/2, base/2];
# - each randomizer has alpha/2 number of 1s and alpha/2 number of -1s
#
# after applying randomizer, the result, denoted by `b := a r` has
# coefficients which are the accumulation of alpha numbers between -base/2 and
# base/2.
#
# at last we are adding rho different signatures together,
# the final bound, beta_agg, is a binomial distribution of alph * rho elements,
# scaled by base/2

chi = {-1: 1/2, 1: 1/2}
chi_alpha = iter_law_convolution(chi, alpha*rho)
(beta_agg, tail_prob_bits) = tail_param(chi_alpha, 40)  # use 40 here because beta_agg will be close to a power-of-2
beta_agg = beta_agg * base_over_2
log_beta_agg = ZZ(ceil(log(beta_agg, 2)))

# a single coefficient has 2^-40 probablity to be greater than beta_agg.
# here an union bound for n * rho * tau gives roughly 2^-16 fail rate
# which is probably fine? or maybe we want even larger failure rate?
fail_prob = ceil(log(n * rho * tau, 2) + tail_prob_bits)

# =================================
# estimate the security of the scheme

## using same root hermite factor as before
delta = 1.004

## best m for attacker
m = ceil(n * log(hvc_q,2) * log(delta, 2))

## infinity norm of lambda_1 * lambda^dim
## page 6 of "pqc paper" MR08
capable_to_find = 2^(2*sqrt(n * log(hvc_q, 2) * log(delta, 2)))
print (capable_to_find, beta_agg, capable_to_find > beta_agg)

# =================================
# now let's try to get the sizes
hvc_ring_element_size = hvc_log_q * n

fresh_sig_size = hvc_log_q * n * tau
agg_sig_size = (log_beta_agg + 1) * n * width * tau * 2

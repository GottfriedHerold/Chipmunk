# decompose a ring element a with coefficients uniformly in 0 and q
# into a list of ring elements with coefficients in {0,base} for base >= 1

# toy example parameters

# modulus
# we use falcon modulus 2^13 + 2^12 + 1
q = 12289
base = 8
half_base = 4
width = ZZ(round(log(q, base)))

# dimension
n = 512


P.<x> = PolynomialRing(Zmod(q))

# input f in R
# output (f0, ... f15) s.t. f0 + f1 * 2 + f2 * 4 + ... + f15 * 2^15 = f
def decompose(input):
    '''
    input: [u16; n], an array of u16 for dimension 512
    output: [[u16; n]; width], width such arrays; i.e., width ring elements
    '''
    assert(len(input)==n)
    input = vector(input)
    res = []
    for i in range (width):
        current = input % base
        for j in range(n):
            if current[j] >= half_base:
                current[j] = current[j] - base

        input = (input - current)/base
        res.append(current)
    return res


# input (f0, ... f15)
# output f = f0 + f1 * 2 + f2 * 4 + ... + f15 * 2^15
def project(input):
    '''
    input: [[u16; n]; width] width arrays of u16 for dimension n
    output: [u16; n], an array representation of a ring element
    '''
    assert(len(input)==width)
    res = input[0]
    for i in range(1,width):
        res = res + input[i] * base^i

    return vector(res)

def aggregate(input1, input2):
    '''
    input1: [[u16; n]; width], width arrays of u16 for dimension n
    input2: [[u16; n]; width], width arrays of u16 for dimension n
    output: [[u16; n]; width], width arrays of u16 for dimension n
    '''
    assert(len(input1)==width)
    assert(len(input2)==width)
    res = []
    for i in range(width):
        cur = input1[i]+ input2[i]
        res.append(cur)
    return res

# ring multiplication for R_q elements
def ring_mul(input1, input2):
    '''
    input1: [u16; n], representing a ring element
    input2: [u16; n], representing a ring element
    output: [u16; n], representing a ring element
    '''
    assert(len(input1)==n)
    assert(len(input2)==n)
    a = P(input1)
    b = P(input2)
    c = a*b % P(x^n +1)
    c = c.list()
    for i in range(len(c), n):
        c.append(0)
    return vector(c)

# ring multiplication for Decompose(R_q) * Rq element
#
# here we will need a new method to handle this multiplication
# if we view input1 as a list of 512 * 16 integers that is trunked into
# 512 blocks, each block is of 16 bits, denoted by
#   [trunk[0], ... trunk[512]]
# then ring multiplication becomes treating each trunk as a coefficients
# of the polynomial, and replacing integer x integer multiplications with
# trunk x integer multiplication
#
def ring_mul_decomposed(input1, input2):
    '''
    input1: [[u16; n]; width], representing a decomposed ring element
    input2: [u16; n], representing a ring element
    output: [[u16; n]; width], representing a decomposed ring element
    '''
    assert(len(input1)==width)
    assert(len(input2)==n)

    res = []
    for e in input1:
        res.append(ring_mul(e.list(), input2))

    return res


def hash(param_a, input_x):
    '''
    param_a: a list of log(q) = witdh ring elements
    input_x: a list of log(q) = witdh ring elements of small norm
    return \\sum a_i * x_i
    '''
    assert(len(param_a) == len(input_x))
    res = 0
    for i in range (len(param_a)):
        res += P(param_a[i].list()) * P(input_x[i].list()) % P(x^n + 1)
    res = vector(res.list()) % q
    return res

#
#
#
test_repeat = 10
#
# =======================================
# test for decomposition correctness
# =======================================
sample1 = vector([ZZ(i) for i in range(n)])
tmp = decompose(sample1)
sample1_rec = project(tmp)
assert(sample1 == sample1_rec)

# print (tmp, sample1, sample1_rec)


# random testing
for _ in range (test_repeat):
    sample2 = vector([ZZ.random_element(0, q) for i in range(n)])
    tmp = decompose(sample2)
    sample2_rec = project(tmp)
    assert(sample2 == sample2_rec)

# =======================================
# test for decomposed multiplication correctness
# =======================================

# first input is a small element decomposed into binary form
sample1 = [ZZ.random_element(0, q) for i in range(n)]
sample1_dec = decompose(sample1)

# second input is the randomizer which is __ternary__
sample2 = [ZZ.random_element(0, 3) - 1 for i in range(n)]
sample_1_mul_2 = ring_mul(sample1, sample2)

# use the `decomposed multiplication`
sample_1_dec_mul_2 = ring_mul_decomposed(sample1_dec, sample2)
sample_1_mul_2_rec = project(sample_1_dec_mul_2)

# check if sample_1_mul_2 and sample_1_mul_2_rec represent a
# same ring element
sample_1_mul_2_rec = sample_1_mul_2_rec % q
assert(sample_1_mul_2 == sample_1_mul_2_rec)


# =======================================
# test for decomposition homomorphism
# =======================================
# test for homomorphic additions for 2 elements
for _ in range(test_repeat):
    sample3 = vector([ZZ.random_element(0, q) for i in range(n)])
    sample3_dec = decompose(sample3)
    sample3_rec = project(sample3_dec)
    assert(sample3 == sample3_rec)

    sample4 = vector([ZZ.random_element(0, q) for i in range(n)])
    sample4_dec = decompose(sample4)
    sample4_rec = project(sample4_dec)
    assert(sample4 == sample4_rec)

    sum_sample_3_4 = vector([(sample3[i] + sample4[i])%q for i in range(n)])
    sum_sample_3_4_dec = aggregate(sample3_dec, sample4_dec)
    sum_sample_3_4_rec = project(sum_sample_3_4_dec)

    # check if sum_sample_3_4 and sum_sample_3_4_rec represent a
    # same ring element
    sum_sample_3_4_rec = sum_sample_3_4_rec%q
    assert(sum_sample_3_4 == sum_sample_3_4_rec)



# =======================================
# test for hash function
# =======================================
input_x = []
param_a = []
for _ in range(width):
    a_i = vector([ZZ.random_element(0, q) for i in range(n)])
    x_i = vector([ZZ.random_element(0, 2) for i in range(n)])
    param_a.append(a_i)
    input_x.append(x_i)

y = hash(param_a, input_x)


# =======================================
# test for hash function with ternary polynomial randomizer
# =======================================

## user 1's input
a = vector([ZZ.random_element(0, q) for i in range(n)])
a_decom = decompose(a)


## user 2's input
b = vector([ZZ.random_element(0, q) for i in range(n)])
b_decom = decompose(b)

## public parameter to the hash
param_a = []
for _ in range(width):
    a_i = vector([ZZ.random_element(0, q) for i in range(n)])
    param_a.append(a_i)

##
## user 1's randomizer
w1 = vector([ZZ.random_element(0, 3) - 1 for i in range(n)])
## user 2's randomizer
w2 = vector([ZZ.random_element(0, 3) - 1 for i in range(n)])

## first case:
## - apply the randomizer with "decomposite form multiplication" first
## - then hash
hash_input_decomp = []
user1 = ring_mul_decomposed(a_decom, w1.list())
user2 = ring_mul_decomposed(b_decom, w2.list())
for i in range(width):
    hash_input_decomp.append( user1[i] + user2[i] )

hash_output = hash(param_a, hash_input_decomp)

## second case:
## - directly do the hash
## - apply the randomizer to the digests
user1_hash = hash(param_a, a_decom)
user2_hash = hash(param_a, b_decom)
hash_output_rec = (P(w1.list()) * P(user1_hash.list()) + P(w2.list()) * P(user2_hash.list())) % P(x^n+1)
hash_output_rec = vector(hash_output_rec)

assert(hash_output == hash_output_rec)



# =======================================
# test for hash function with binary polynomial randomizer, and more users
# =======================================

## user 1's input
a = vector([ZZ.random_element(0, q) for i in range(n)])
a_decom = decompose(a)


## user 2's input
b = vector([ZZ.random_element(0, q) for i in range(n)])
b_decom = decompose(b)

## user 3's input
c = vector([ZZ.random_element(0, q) for i in range(n)])
c_decom = decompose(c)


## user 4's input
d = vector([ZZ.random_element(0, q) for i in range(n)])
d_decom = decompose(d)


## public parameter to the hash
param_a = []
for _ in range(width):
    a_i = vector([ZZ.random_element(0, q) for i in range(n)])
    param_a.append(a_i)

## in fact, you can change the norm of randomizer from
## "2" to something like 100
## and the code says homomorphism still holds
##
## user 1's randomizer
w1 = vector([ZZ.random_element(0, 3)-1 for i in range(n)])
## user 2's randomizer
w2 = vector([ZZ.random_element(0, 3)-1 for i in range(n)])
## user 3's randomizer
w3 = vector([ZZ.random_element(0, 3)-1 for i in range(n)])
## user 4's randomizer
w4 = vector([ZZ.random_element(0, 3)-1 for i in range(n)])

## first case:
## - apply the randomizer with "decomposite form multiplication" first
## - then hash
hash_input_decomp = []
user1 = ring_mul_decomposed(a_decom, w1.list())
user2 = ring_mul_decomposed(b_decom, w2.list())
user3 = ring_mul_decomposed(c_decom, w3.list())
user4 = ring_mul_decomposed(d_decom, w4.list())
for i in range(width):
    hash_input_decomp.append( user1[i] + user2[i] + user3[i] + user4[i] )

hash_output = hash(param_a, hash_input_decomp)

## second case:
## - directly do the hash
## - apply the randomizer to the digests
user1_hash = hash(param_a, a_decom)
user2_hash = hash(param_a, b_decom)
user3_hash = hash(param_a, c_decom)
user4_hash = hash(param_a, d_decom)
hash_output_rec = (P(w1.list()) * P(user1_hash.list()) + P(w2.list()) * P(user2_hash.list()) + P(w3.list()) * P(user3_hash.list()) + P(w4.list()) * P(user4_hash.list())) % P(x^n+1)
hash_output_rec = vector(hash_output_rec)

assert(hash_output == hash_output_rec)


print("finished")

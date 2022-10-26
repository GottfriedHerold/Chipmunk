# decompose a ring element a with coefficients uniformly in 0 and q
# into a list of ring elements with coefficients in {0,1}

# toy example parameters

# modulus
# let's pick a power-of-2 q so we can argue the decomposed coefficients
# are uniformly distributed
q = 2^16

# dimension
n = 512

P.<x> = PolynomialRing(Zmod(q))

def decompose(input):
    '''
    input: [u16; 512], an array of u16 for dimension 512
    output: [[u16; 512]; 16], 16 such arrays; i.e., 16 ring elements
    '''
    assert(len(input)==512)
    res = []
    for i in range (16):
        current = []
        for j in range(i*32, (i+1)*32):
            # print(i, j, input[j])
            coeff_bits = int_decom(input[j])
            current += coeff_bits
        res.append(current)

    return res


def int_decom(input):
    '''
    input: u16, an integer
    output: [u16, 16], bit decompose of the input
    '''
    res = input.bits();
    for _ in range(len(res), 16):
        res.append(0)

    return res


def int_project(input):
    '''
    input: [u16, 16], bit decompose of the input where coefficients
        are not necessarily in {0,1} after aggregations
    output: integer, may overflow 2^16 after aggregations
    '''
    assert(len(input)==16)
    res = 0
    for i in range(16):
        res += input[i]*2^i
    return res

def project(input):
    '''
    input: [[u16; 512]; 16] 16 arrays of u16 for dimension 512
    output: [u16; 512], an array representation of a ring element
    '''
    assert(len(input)==16)
    res = []
    for e in input:
        # each array is 32 coefficients
        for i in range(32):
            res.append(int_project(e[i*16:(i+1)*16]))

    return res

def aggregate(input1, input2):
    '''
    input1: [[u16; 512]; 16], 16 arrays of u16 for dimension 512
    input2: [[u16; 512]; 16], 16 arrays of u16 for dimension 512
    output: [[u16; 512]; 16], 16 arrays of u16 for dimension 512
    '''
    assert(len(input1)==16)
    assert(len(input2)==16)
    res = []
    for i in range(16):
        cur = [input1[i][j] + input2[i][j] for j in range(512)]
        res.append(cur)
    return res

# ring multiplication for R_q elements
def ring_mul(input1, input2):
    '''
    input1: [u16; 512], representing a ring element
    input2: [u16; 512], representing a ring element
    output: [u16; 512], representing a ring element
    '''
    assert(len(input1)==512)
    assert(len(input2)==512)
    a = P(input1)
    b = P(input2)
    c = a*b % P(x^512 +1)
    c = c.list()
    for i in range(len(c), 512):
        c.append(0)
    return c

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
    input1: [[u16; 512];16], representing a decomposed ring element
    input2: [u16; 512], representing a ring element
    output: [[u16; 512];16], representing a decomposed ring element
    '''
    assert(len(input1)==16)
    assert(len(input2)==512)
    input_concat = []
    for e in input1:
        input_concat += e
    res_concat = [0 for _ in range(512*16*2)]

    for i in range(512):
        for j in range(512):
            tmp = trunk_int_mul(input_concat[i*16: (i+1)*16], input2[j])
            for e in range((i+j)*16, (i+j+1)*16):
                res_concat[e] += tmp[e - (i+j)*16]

    for i in range(512*16):
        res_concat[i] -= res_concat[i+512*16]

    res = []
    for i in range(16):
        res.append(res_concat[i*512:(i+1)*512])
    return res


def trunk_int_mul(input1, input2):
    '''
    input1: [u16;16], representing a decomposed coefficient of a ring element
    input2: int, representing a coefficient of a ring element
    output: [int;16], representing a decomposed coefficient of a ring element
    '''
    assert(len(input1)==16)
    output = [e*input2 for e in input1]
    return output



test_repeat = 1

# =======================================
# test for decomposition correctness
# =======================================
sample1 = [ZZ(i) for i in range(512)]
tmp = decompose(sample1)
sample1_rec = project(tmp)
assert(sample1 == sample1_rec)

# random testing
for _ in range (test_repeat):
    sample2 = [ZZ.random_element(0, q) for i in range(512)]
    tmp = decompose(sample2)
    sample2_rec = project(tmp)
    assert(sample2 == sample2_rec)

# =======================================
# test for decomposed multiplication correctness
# =======================================

# first input is a small element decomposed into binary form
sample1 = [ZZ.random_element(0, 1000) for i in range(512)]
sample1_dec = decompose(sample1)
# second input is the randomizer which is binary
sample2 = [ZZ.random_element(0, 2) for i in range(512)]
sample_1_mul_2 = ring_mul(sample1, sample2)
# use the `decomposed multiplication`
sample_1_dec_mul_2 = ring_mul_decomposed(sample1_dec, sample2)
sample_1_mul_2_rec = project(sample_1_dec_mul_2)

# check if sample_1_mul_2 and sample_1_mul_2_rec represent a
# same ring element
sample_1_mul_2_rec = [e%q for e in sample_1_mul_2_rec]
assert(sample_1_mul_2 == sample_1_mul_2_rec)


# =======================================
# test for decomposition homomorphism
# =======================================
# test for homomorphic additions for 2 elements
for _ in range(test_repeat):
    sample3 = [ZZ.random_element(0, q) for i in range(512)]
    sample3_dec = decompose(sample3)
    sample3_rec = project(sample3_dec)
    assert(sample3 == sample3_rec)

    sample4 = [ZZ.random_element(0, q) for i in range(512)]
    sample4_dec = decompose(sample4)
    sample4_rec = project(sample4_dec)
    assert(sample4 == sample4_rec)

    sum_sample_3_4 = [(sample3[i] + sample4[i])%q for i in range(512)]
    sum_sample_3_4_dec = aggregate(sample3_dec, sample4_dec)
    sum_sample_3_4_rec = project(sum_sample_3_4_dec)

    # check if sum_sample_3_4 and sum_sample_3_4_rec represent a
    # same ring element
    sum_sample_3_4_rec = [e%q for e in sum_sample_3_4_rec]
    assert(sum_sample_3_4 == sum_sample_3_4_rec)

# test for homomorphic additions for 1000 elements
for _ in range(test_repeat):
    samples = []
    decomposes = []
    for _ in range(1000):
        sample = [ZZ.random_element(0, q) for i in range(512)]
        sample_dec = decompose(sample)
        sample_rec = project(sample_dec)
        assert(sample == sample_rec)

        samples.append(sample)
        decomposes.append(sample_dec)

    sum_samples = samples[0]
    sum_samples_dec = decomposes[0]

    for ctr in range(1,1000):
        sum_samples = [(sum_samples[i] + samples[ctr][i])%q for i in range(512)]
        sum_samples_dec = aggregate(sum_samples_dec, decomposes[ctr])

    # check if sum_samples and sum_samples_dec represent a
    # same ring element
    sum_samples_rec = project(sum_samples_dec)
    sum_samples_rec = [e%q for e in sum_samples_rec]
    assert(sum_samples == sum_samples_rec)


# test for homomorphic additions for 10 ring elements with random scalars
# where the scalars are binary polynomials
# Note that this is super slow... so we only setup for 10 users
for _ in range(test_repeat):
    samples = []
    decomposes = []
    scalars = []
    for _ in range(10):
        sample = [ZZ.random_element(0, q) for i in range(512)]
        sample_dec = decompose(sample)
        sample_rec = project(sample_dec)
        assert(sample == sample_rec)

        samples.append(sample)
        decomposes.append(sample_dec)

        scalar = [ZZ.random_element(0, 2) for i in range(512)]
        scalars.append(scalar)

    sum_samples = ring_mul(samples[0], scalars[0])
    sum_samples_dec = ring_mul_decomposed(decomposes[0], scalars[0])

    for ctr in range(1,10):
        tmp = ring_mul(samples[ctr], scalars[ctr])
        sum_samples = [(sum_samples[i] + tmp[i])%q for i in range(512)]
        tmp = ring_mul_decomposed(decomposes[ctr], scalars[ctr])
        sum_samples_dec = aggregate(sum_samples_dec, tmp)


    # check if sum_samples and sum_samples_dec represent a
    # same ring element
    sum_samples_rec = project(sum_samples_dec)
    sum_samples_rec = [e%q for e in sum_samples_rec]
    assert(sum_samples == sum_samples_rec)

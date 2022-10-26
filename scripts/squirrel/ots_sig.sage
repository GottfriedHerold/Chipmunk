import hashlib

#########################################
########## toy parameters ###############
# modulus
q = 8380417
# dimension of the polynomial
n = 16
# norm bound for the noise
beta = 1
# number of polynomials in a vector
l = 3
#########################################
########## setup ########################
Rcal = Zmod(q)
P.<x> = PolynomialRing(Rcal)
F = P(x^n+1)
PQ.<x> = PolynomialQuotientRing(P, F)

#########################################
# scheme
# following the definition of `agg_ots` paper

def PrmsGen():
    a = []
    for _ in range (l):
        f = PQ([Rcal.random_element() for _ in range (n)])
        a.append(f)
    return a

def KeyGen(pp):
    s0 = []
    s1 = []
    v0 = 0
    v1 = 0
    for i in range (l):
        f1 = PQ([ZZ.random_element(0, beta+1) for _ in range (n)])
        f2 = PQ([ZZ.random_element(0, beta+1) for _ in range (n)])
        v0 += f1 * pp[i]
        v1 += f2 * pp[i]
        s0.append(f1)
        s1.append(f2)
    return ((s0,s1), (v0, v1))

def Sign(sk, m):
    c = Hash0(m)
    sigma = [ c*sk[0][i] + sk[1][i] for i in range(l)]
    return sigma

def SigAgg(pk_list, msg_list, sig_list):
    t = Hash1(pk_list, msg_list)

    sig_agg = [0 for _ in range(l)]
    for e in range(len(sig_list)):
        for i in range (l):
            sig_agg[i] += t[e] * sig_list[e][i]

    return sig_agg

def Verify(pp, pk_list, msg_list, sig):

    c_list = []
    for m in msg_list:
        c = Hash0(m)
        c_list.append(c)
    t = Hash1(pk_list, msg_list)

    # verify a^t sigma = \sum (v_{i,0} c_i + v_{i,1}) * t_i
    left = 0
    for e in range (l):
        left += pp[e] * sig[e]

    right = 0
    for e in range (len(pk_list)):
        right += (pk_list[e][0]*c_list[e] + pk_list[e][1])*t[e]

    return left == right

# additional functions
def Hash0(msg):
    '''
    hash into a polynomial of small norm, rather than a constant
    __this is different from the paper__
    '''
    hasher = hashlib.sha256()
    hasher.update(msg)
    s = ZZ("0x" + hasher.hexdigest()).binary()
    res = []
    for i in range(n):
        res.append(ZZ(s[i]))
    return PQ(res)

def Hash1(pk_list, msg_list):
    '''
    hash into a list of polynomials of small norm, rather than a
    list of constants
    __this is different from the paper__
    '''
    if (len(pk_list) == 1):
        return PQ(1)


    assert(len(pk_list) == len(msg_list))

    hash_input = b''
    for e in range(len(pk_list)):
        for pk in pk_list:
            for v_0_i in pk[0]:
                hash_input += bytes((ZZ(v_0_i)%256).str(), 'utf-8')
            for v_0_i in pk[1]:
                hash_input += bytes((ZZ(v_0_i)%256).str(), 'utf-8')
        hash_input += msg_list[e]

    hasher = hashlib.sha256()
    hasher.update(hash_input)
    s =  ZZ("0x" + hasher.hexdigest())
    set_random_seed(s)
    t = []
    for e in range(len(pk_list)):
        tmp = PQ([ZZ.random_element(0, beta+1) for _ in range(n)])
        t.append(tmp)
    return t


#########################################
# testing
pp = PrmsGen()
(sk1, pk1) = KeyGen(pp)
msg1 = b"42"
sig1 = Sign(sk1, msg1)

(sk2, pk2) = KeyGen(pp)
msg2 = b"43"
sig2 = Sign(sk2, msg2)

print(Verify(pp, [pk1], [msg1], sig1))
print(Verify(pp, [pk2], [msg2], sig2))

sig_agg = SigAgg([pk1, pk2], [msg1, msg2], [sig1, sig2])
print(Verify(pp, [pk1, pk2], [msg1, msg2], sig_agg))

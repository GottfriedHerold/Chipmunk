load("hom_hash.sage")



#########################################
########## toy parameters ###############
# height of the tree
#           rt
#        m10, m11
#   l20, l21, l22, l23
# lxx are the leaves and are the (hashes of) pks for ots-sig
h = 2

def new_tree():
    return [0, [0,0], [0,0,0,0]]


def gen_tree(leaves, param):
    assert(len(leaves) == 4)
    (L, R) = param
    m10 = hash(leaves[0], leaves[1], L, R)
    m11 = hash(leaves[2], leaves[3], L, R)
    m10 = bit_decompose(m10)
    m11 = bit_decompose(m11)
    rt = hash(m10, m11, L, R)

    return [rt, [m10, m11], leaves]


def bit_decompose(a):
    # todo!
    # it should return an l-dim vector of poly with small coeffs
    # for now, we simply return itself: a k-dim vector of poly
    return a



param = parameters()

# leaves is a vector of 4 elements, each element is an l-dim vector
# of polynomials of degree n
leaves1 = [[ PQ([ZZ.random_element(0, beta+1) for _ in range (n)]) for _ in range (l)] for _ in range(4)]
t1 = gen_tree(leaves1, param)

# leaves is a vector of 4 elements, each element is an l-dim vector
# of polynomials of degree n
leaves2 = [[ PQ([ZZ.random_element(0, beta+1) for _ in range (n)]) for _ in range (l)] for _ in range(4)]
t2 = gen_tree(leaves2, param)


leaves_sum = []
for i in range (4):
    leaf = [ leaves1[i][j] + leaves2[i][j] for j in range (l)]
    leaves_sum.append(leaf)
t_sum = gen_tree(leaves, param)


# this code doesnot work because bit_decompose is not implemented properly
# m10 and m11 are too large so that it is no longer additive homomorphic
assert(t1[0][0] + t2[0][0] == t_sum[0][0])



# - rev(i) is the reverse bit decomposition of i, i.e.,
#   0   ->  0
#   1   ->  100 0000
#   2   ->  010 0000
#   3   ->  110 0000   ...
def reverse_bits(i, n):
    t = i.binary()[::-1]

    while len(t) < n:
        t = t + "0"

    res = 0
    for e in t:
        res *= 2
        res += ZZ(e)
    return res

def print_hots_ntt():
    q_hots = 0x305801

    P.<x> = PolynomialRing(Zmod(q_hots))
    f = P(x^1024+1)
    r = f.roots()[0][0]
    r_inv = 1/r
    print(r)

    for i in range (1024):
        e = reverse_bits(ZZ(i), 10)
        t = ZZ(r^e)
        if t > 0x305801/2:
            t = t - 0x305801
        print(t, end = ', ')
        # print(i, e, r^e)
    print()

def print_hots_inv_ntt():
    q_hots = 0x305801

    P.<x> = PolynomialRing(Zmod(q_hots))
    f = P(x^1024+1)
    r = f.roots()[0][0]
    r_inv = 1/r
    print(r)

    for i in range (1024):
        e = reverse_bits(ZZ(i), 10)
        t = ZZ(r_inv^e)
        if t > 0x305801/2:
            t = t - 0x305801
        print(t, end = ', ')
    print()

def print_hvc_ntt():
    q_hvc = 202753

    P.<x> = PolynomialRing(Zmod(q_hvc))
    f = P(x^1024+1)
    r = f.roots()[0][0]
    r_inv = 1/r
    print(r)

    for i in range (1024):
        e = reverse_bits(ZZ(i), 10)
        t = ZZ(r^e)
        if t > 202753/2:
            t = t - 202753
        print(t, end = ', ')
        # print(i, e, r^e)
    print()

def print_hvc_inv_ntt():
    q_hvc = 202753

    P.<x> = PolynomialRing(Zmod(q_hvc))
    f = P(x^1024+1)
    r = f.roots()[0][0]
    r_inv = 1/r
    print(r)

    for i in range (1024):
        e = reverse_bits(ZZ(i), 10)
        t = ZZ(r_inv^e)
        if t > 202753/2:
            t = t - 202753
        print(t, end = ', ')
    print()


print_hots_ntt()
print()
print_hots_inv_ntt()
print()
print_hvc_ntt()
print()
print_hvc_inv_ntt()

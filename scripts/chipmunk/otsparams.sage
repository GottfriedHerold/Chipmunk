from tabulate import tabulate
load("dist.sage")

def cardinalitySetOfTernaryPoly(n,alpha):
  return binomial(n,alpha)*2^alpha

def findHammingWeight(n,l):
  for alpha in range(n+1):
    if cardinalitySetOfTernaryPoly(n,alpha) > l:
      return alpha
  raise ValueError("There does not exist a Hamming weight satisfying the specified conditions.")

#def slowFindGamma(secpar,delta,n,q,phi):
#  for gamma in range(1,101):
#    if 2^((3*secpar+delta)/(n*gamma))*q^(1/gamma) <= phi+.5:
#      return gamma

def findGamma(secpar,delta,n,q,phi):
  u = 50
  l = 1
  while u!=l:
    gamma = l+(floor((u-l)/2))
    if 2^((3*secpar+delta)/(n*gamma))*q^(1/gamma) <= phi+.5:
      if 2^((3*secpar+delta)/(n*(gamma-1)))*q^(1/(gamma-1)) > phi+.5:
        return gamma
      else:
        u = gamma
    else:
      l = gamma
  raise ValueError("Could not find gamma satisfying the specified conditions.")

def findGammaROM(secpar,n,q,phi):
  u = 100
  l = 1
  while u!=l:
    gamma = l+(floor((u-l)/2))
    if 2^((2*secpar)/(n*gamma))*q^(1/gamma) <= phi+.5:
      if 2^((2*secpar)/(n*(gamma-1)))*q^(1/(gamma-1)) > phi+.5:
        return gamma
      else:
        u = gamma
    else:
      l = gamma
  raise ValueError("Could not find gamma satisfying the specified conditions.")

def findPrime(beta):
  q = next_prime(2*beta)
  while q % (2*n) != 1:
    q = next_prime(q)
  return q

secpar = 128
n = 512
rho = 1024
hvc_arity = 23
hvc_beta_agg = 11231
c=1.005

# Hamming weight for random weights
alpha_w = findHammingWeight(n,2^secpar)

# Hamming weight for hash function outputs
alpha_H = findHammingWeight(n,2^(2*secpar))

# The range of the hash function is at most 2^delta times larger than needed. Usually delta=1.
delta = ceil(log(cardinalitySetOfTernaryPoly(n,alpha_H),2))-2*secpar

list = []

chi = {-1: 1/2, 1: 1/2}
chi_alpha = iter_law_convolution(chi, alpha_w*rho)
for phi in range(1,15):
  (beta_sig, tail_prob_bits) = tail_param(chi_alpha, 40)  # use 40 here because beta_agg will be close to a power-of-2
  beta_sig = beta_sig * 2*phi*alpha_H

#  beta_sig  = (2*rho*alpha_w*alpha_H)*phi
#  beta_KOTS = ((4*rho+4)*alpha_w*alpha_H)*phi
  beta_KOTS = 2*beta_agg+4*alpha_w*phi*alpha_H
  q = findPrime(beta_KOTS)
  gamma = findGamma(secpar,delta,n,q,phi)
#  gamma = findGammaROM(secpar,n,q,phi)
  fail_prob = ceil(log(n * rho * gamma, 2) + tail_prob_bits)
  size = gamma*n*ceil(log(2*beta_sig+1,2))/8/1024
#  sizeROM = gammaROM*n*ceil(log(2*beta_sig+1,2))/8/1024
  size_key = (ceil(log(2*hvc_beta_agg+1,2))*n*ceil(log(q,hvc_arity)))/8/1024
  
  check_one = false
  if 2*beta_KOTS < c^(n*gamma)*(q*10)^(1/(n*gamma)) -1:
    check_one = true
  
  check_two = false
  if sqrt(n*gamma)*beta_KOTS < c^(n*gamma)*sqrt(n*gamma/(2*pi*e))*q^(1/(gamma)):
    check_two = true
  
  list.append([phi,str(q),gamma,str((size.n())) + " KB",str((size_key.n())) + " KB",str((size+size_key).n()) + " KB",check_one,check_two,fail_prob])#,str(ceil(sizeROM)) + " KB"])

print("alpha_w: ", alpha_w)
print("alpha_H: ", alpha_H)
print("delta: ", delta)

print(tabulate(list,headers=["phi","q","gamma","size sig", "size key", "total size","check 1","check 2","fail_prob"],disable_numparse=True))

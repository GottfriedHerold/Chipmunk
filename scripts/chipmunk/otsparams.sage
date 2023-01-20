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

# security parameter
secpar = 128
# polynomial degree
n = 512
# number of users
rho = 1024
# height of the tree
tau = 21

# Hamming weight for random weights
alpha_w = findHammingWeight(n,2^secpar)
c=1.005

chi = {-1: 1/2, 1: 1/2}
chi_alpha = iter_law_convolution(chi, alpha_w*rho)

list_hvc = []

for hvc_arity in range(3,24,2):
  (beta_path, tail_prob_bits_path) = tail_param(chi_alpha, 32)
  beta_path = beta_path * (hvc_arity-1)/2

  beta_hvc = 4*beta_path

  q_hvc = findPrime(beta_hvc)
  width = ZZ(ceil(log(q_hvc, hvc_arity)))

  fail_prob_path = ceil(log(n * width * tau, 2) + tail_prob_bits_path)

  check_one_hvc = bool(2*beta_hvc < c^(n*width)*(q_hvc*10)^(1/(n*width)) -1)

  check_two_hvc = bool(sqrt(n*width)*beta_hvc < c^(n*width)*sqrt(n*width/(2*pi*e))*q_hvc^(1/(width)))

  size_path = (ceil(log(beta_path*2+1,2))*n*width*2*tau)/8/1024
  
  list_hvc.append([hvc_arity,q_hvc,beta_path,width,size_path,check_one_hvc,check_two_hvc,fail_prob_path])

print(tabulate(list_hvc,headers=["arity","q_hvc","beta_agg","width","size","check 1","check2","fail prob"]))

# Hamming weight for hash function outputs
alpha_H = findHammingWeight(n,2^(2*secpar))

# The range of the hash function is at most 2^delta times larger than needed. Usually delta=1.
delta = ceil(log(cardinalitySetOfTernaryPoly(n,alpha_H),2))-2*secpar

i = int(input("Choose Arity:"))

#print("arity:",list_hvc[int((i-1)/2)-1])

hvc_arity = list_hvc[int((i-1)/2)-1][0]
hvc_beta_agg = list_hvc[int((i-1)/2)-1][2]
size_path = list_hvc[int((i-1)/2)-1][4]


list = []

for phi in range(10,20):
  (beta_sig, tail_prob_bits) = tail_param(chi_alpha, 28)  # use 40 here because beta_agg will be close to a power-of-2
  beta_sig = beta_sig*2*phi*alpha_H

#  beta_sig  = (2*rho*alpha_w*alpha_H)*phi
#  beta_KOTS = ((4*rho+4)*alpha_w*alpha_H)*phi
  beta_KOTS = 2*beta_sig+4*alpha_w*phi*alpha_H
  q = findPrime(beta_KOTS)
  gamma = findGamma(secpar,delta,n,q,phi)
#  gamma = findGammaROM(secpar,n,q,phi)
  fail_prob = ceil(log(n * gamma, 2) + tail_prob_bits)
  size = gamma*n*ceil(log(2*beta_sig+1,2))/8/1024
#  sizeROM = gammaROM*n*ceil(log(2*beta_sig+1,2))/8/1024
  size_key = (ceil(log(2*hvc_beta_agg+1,2))*n*ceil(log(q,hvc_arity)))/8/1024
  
  check_one = bool(2*beta_KOTS < c^(n*gamma)*(q*10)^(1/(n*gamma)) -1)
  
  check_two = bool(sqrt(n*gamma)*beta_KOTS < c^(n*gamma)*sqrt(n*gamma/(2*pi*e))*q^(1/(gamma)))
  
  list.append([phi,str(q),gamma,beta_sig,("%.4f" % size) + " KB",("%.4f" % size_key) + " KB",("%.4f" % (size+size_key+size_path)) + " KB",check_one,check_two,fail_prob])#,str(ceil(sizeROM)) + " KB"])

print("arity: ", hvc_arity)
print("q_hvc: ", list_hvc[int((i-1)/2)-1][1])
print("beta_agg: ", hvc_beta_agg)
print("alpha_w: ", alpha_w)
print("alpha_H: ", alpha_H)
print("delta: ", delta)

print(tabulate(list,headers=["phi","q","gamma","beta_sig","size sig", "size key", "total size","check 1","check 2","fail_prob"],disable_numparse=True))

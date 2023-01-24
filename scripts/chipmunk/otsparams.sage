from tabulate import tabulate
from colorama import init, Fore, Back, Style
load("dist.sage")

def cardinalitySetOfTernaryPoly(n,alpha):
  return binomial(n,alpha)*2^alpha

def findHammingWeight(n,l):
  for alpha in range(n+1):
    if cardinalitySetOfTernaryPoly(n,alpha) > l:
      return alpha
  raise ValueError("There does not exist a Hamming weight satisfying the specified conditions.")

def findGamma(secpar,delta,n,q,phi):
  u = 40
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
  u = 40
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
# Hermite factor?
c=1.005

# Hamming weight for random weights
alpha_w = findHammingWeight(n,2^secpar)

chi = {-1: 1/2, 1: 1/2}
chi_alpha = iter_law_convolution(chi, alpha_w*rho)

hvc_params = {}
hvc_min_params = 0
hvc_min_size = oo

(hvc_tail_cut, hvc_tail_cut_prob) = tail_param(chi_alpha, 33)
for hvc_arity in range(3,40,2):
  beta_path = hvc_tail_cut * (hvc_arity-1)/2
  beta_hvc = 4*beta_path

  q_hvc = findPrime(beta_hvc)
  width = ZZ(ceil(log(q_hvc, hvc_arity)))


#  check_one_hvc = bool(2*beta_hvc < c^(n*width)*q_hvc^(1/width) -1)

#  check_two_hvc = bool(sqrt(n*width)*beta_hvc < c^(n*width)*sqrt(n*width/(2*pi*e))*q_hvc^(1/(width)))

  if bool(sqrt(n*width)*beta_hvc < c^(n*width)*sqrt(n*width/(2*pi*e))*q_hvc^(1/(width))) and bool(2*beta_hvc < c^(n*width)*q_hvc^(1/width) -1):
    fail_prob_path = floor(log(2 * n * width * tau, 2) + hvc_tail_cut_prob)
    hvc_path_size = (ceil(log(beta_path*2+1,2))*n*width*2*tau)/8/1024
    if hvc_path_size < hvc_min_size:
      hvc_min_params = hvc_arity
      hvc_min_size = hvc_path_size
    hvc_params[hvc_arity] = {"arity" : hvc_arity,"q" : q_hvc, "norm bound" : beta_path, "SIS norm bound" : beta_hvc, "SIS width" : width, "size" : hvc_path_size,"failure prob" : fail_prob_path, "tail cut" : hvc_tail_cut, "tail cut prob" : hvc_tail_cut_prob}
  else:
    print("Oh no, ", hvc_arity, " not viable")
    break
to_tabulate = []
for p in hvc_params.values():
  if p["arity"] == hvc_min_params:
    to_tabulate.append([Back.YELLOW+str(p["arity"]),p["q"],p["norm bound"],p["SIS norm bound"],p["SIS width"],("%.4f" % p["size"]) + " KB","2^"+str(p["failure prob"])+Back.RESET])
  else:
    to_tabulate.append([str(p["arity"]),p["q"],p["norm bound"],p["SIS norm bound"],p["SIS width"],("%.4f" % p["size"]) + " KB","2^"+str(p["failure prob"])])
#  p[5] = ("%.4f" % p[5]) + " KB"
print(tabulate(to_tabulate,headers=["arity","q_HVC","beta_agg","beta_HVC","width","size","fail prob"]))

# Hamming weight for hash function outputs
alpha_H = findHammingWeight(n,2^(2*secpar))

# The range of the hash function is at most 2^delta times larger than needed. Usually delta=1.
delta = ZZ(ceil(log(cardinalitySetOfTernaryPoly(n,alpha_H),2))-2*secpar)

chosen_arity = input("Choose arity ["+str(hvc_min_params) + "]:")
if chosen_arity == "":
  hvc = hvc_params[hvc_min_params]
else:
  hvc = hvc_params[int(chosen_arity)]

hvc_arity = hvc["arity"]
hvc_beta_agg = hvc["norm bound"]
bitsize_coeff_beta_agg = ceil(log(2*hvc_beta_agg+1,2))
hvc_path_size = hvc["size"]
hvc_width = hvc["SIS width"]
hvc_tail_cut_prob = hvc["tail cut prob"]


sig_params = {}
sig_min_params = 0
sig_min_size = oo

(sig_tail_cut, tail_prob_bits) = tail_param(chi_alpha, 28)
for phi in range(1,100):
  beta_sig = sig_tail_cut*2*phi*alpha_H

  beta_KOTS = 2*beta_sig+4*alpha_w*phi*alpha_H
  # Second condition comes from Lemma 19 of the squirrel paper
  q = findPrime(max(beta_KOTS,8*alpha_w*alpha_H*phi))
  gamma = findGamma(secpar,delta,n,q,phi)
  key_arity = hvc_arity
  beta_key = hvc_tail_cut*(key_arity-1)/2

  
  check_one = bool(2*beta_KOTS < c^(n*gamma)*q^(1/gamma) -1)
  check_one_key = bool(2*(4*beta_key) < c^(n*ceil(log(q,key_arity)))*hvc["q"]^(1/ceil(log(q,key_arity))) -1)
  
  check_two = bool(sqrt(n*gamma)*beta_KOTS < c^(n*gamma)*sqrt(n*gamma/(2*pi*e))*q^(1/(gamma)))
  check_two_key = bool(sqrt(n*ceil(log(q,key_arity)))*(4*beta_key) < c^(n*ceil(log(q,key_arity)))*sqrt(n*ceil(log(q,key_arity))/(2*pi*e))*hvc["q"]^(1/(ceil(log(q,key_arity)))))
  
  if check_one and check_two and check_one_key and check_two_key:
    sig_fail_prob = ceil(log(n * gamma, 2) + tail_prob_bits)
    hvc_fail_prob = ceil(log(2*n*ceil(log(q,hvc_arity))+2*n*hvc_width*tau,2)+hvc_tail_cut_prob)
    total_fail_prob = max(sig_fail_prob,hvc_fail_prob)+1
    size = gamma*n*ceil(log(2*beta_sig+1,2))/8/1024
    size_key = (ceil(log(2*beta_key+1,2))*2*n*ceil(log(q,key_arity)))/8/1024
    total_size = size + size_key + hvc_path_size
    if total_size < sig_min_size:
      sig_min_params = phi
      sig_min_size = total_size
    sig_params[phi] = {"phi" : phi, "q" : q, "norm bound" : beta_sig, "SIS norm bound" : beta_KOTS, "SIS width" : gamma, "sig size" : size, "key size" : size_key, "total size" : total_size, "sig failure prob" : sig_fail_prob, "hvc failure prob" : hvc_fail_prob, "total failure prob" : total_fail_prob}
    
    
#    list.append([phi,str(q),gamma,beta_sig,("%.4f" % size) + " KB",("%.4f" % size_key) + " KB",("%.4f" % (size+size_key+hvc_path_size)) + " KB",check_one,check_two,sig_fail_prob,hvc_fail_prob,total_fail_prob,check_one_key,check_two_key,key_arity])
  else:
    break
  
print("secpar: ", secpar)
print("n:",n)
print("rho:",rho)
print("tau:",tau)
print("arity: ", hvc_arity)
print("q_HVC: ", hvc["q"])
print("beta_agg: ", hvc_beta_agg)
print("alpha_w: ", alpha_w)
print("alpha_H: ", alpha_H)
print("delta: ", delta)

to_tabulate = []
for p in sig_params.values():
  if p["phi"] == sig_min_params:
    to_tabulate.append([Back.YELLOW+str(p["phi"]),str(p["q"]),str(p["norm bound"]),str(p["SIS norm bound"]),p["SIS width"],("%.4f" % p["sig size"]) + " KB",("%.4f" % p["key size"]) + " KB",("%.4f" % p["total size"]) + " KB" ,"2^"+str(p["total failure prob"])+Back.RESET])
  else:
    to_tabulate.append([str(p["phi"]),str(p["q"]),str(p["norm bound"]),str(p["SIS norm bound"]),p["SIS width"],("%.4f" % p["sig size"]) + " KB",("%.4f" % p["key size"]) + " KB",("%.4f" % p["total size"]) + " KB" ,"2^"+str(p["total failure prob"])])
print(tabulate(to_tabulate,headers=["phi","q_KOTS","beta_sigma","beta_KOTS","gamma","signature size", "key size", "total size","total fail prob"]))

from tabulate import tabulate
from colorama import init, Fore, Back, Style
import cProfile
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
  union_bound=3*secpar+delta
  while u!=l:
    gamma = l+(floor((u-l)/2))
    if 2^(union_bound/(n*gamma))*q^(1/gamma) <= phi+.5:
      if 2^(union_bound/(n*(gamma-1)))*q^(1/(gamma-1)) > phi+.5:
        return gamma
      else:
        u = gamma
    else:
      l = gamma
  raise ValueError("Could not find gamma satisfying the specified conditions.")

def findGamma(secpar,delta,n,q,phi):
  log_q = log(q,2)
  log_phi_plus_one_half = log(phi+.5,2)
  return ceil((((3*secpar+delta)/n)+log_q)/log_phi_plus_one_half)
#  raise ValueError("Could not find gamma satisfying the specified conditions.")

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

def rsisIsHard(beta, q, n, m, delta):
  """ The sis problem is hard for the inputs
  :param beta: infinity norm of the short solution
  :param q: modulus
  :param n: dimension of the polynomial ring
  :param m: number of ring elements
  :param delta: root Hermite factor
  """

  # check that sis problem is not trivial
  if beta >= q/2:
    return (false, "SIS is trivial")

  # find the best m for the attacker
  # if there are more m than required, the attacker can ignore some of them
  m = min(m, ceil(sqrt(n * log(q, 2) * log(delta, 2))))

  # check that it is not possible to find short solution, using infinity norm
  if 2 * beta >= delta^(n * m) * q^(1/m) - 1:
    return (false, "infinity norm fails")

  # check that it is not possible to find short solution, using l2 norm
  if sqrt(n * m) * beta >= delta^(n * m) * sqrt( n * m / 2 / pi / e) * q^(1/m):
    return (false, "l2 norm fails")

  return (true, "pass")

def getRootHermiteFactor(secpar):
  """ Return the handwaving root Hermite factor for the given security level
  """
  if secpar == 128:
    return 1.004
  else:
    if secpar == 112:
      return 1.005
    else:
      raise ValueError("Input security parameter should be either 112 or 128")

def find_param_hvc_mcdiarmid(n, secpar, rho, tau, alpha_w, fail_prob_target, verbose):

  # root hermite factor
  c = getRootHermiteFactor(secpar)

  hvc_params = {}
  hvc_min_params = 0
  hvc_min_size = oo

  hvc_arity = 1
  sis_is_hard = true
  while sis_is_hard:
    hvc_arity+=2
    
    log_width = 1
    width_too_big = True
    beta_fresh = int((hvc_arity-1)/2)
    while width_too_big:
      log_width += 1
      epsilon = (fail_prob_target + 1) + 1 + log(n*tau,2) + log_width
      beta_path = ceil(beta_fresh*sqrt((epsilon + 1) * log(2) * 2 * rho * alpha_w))
      beta_hvc = 4 * beta_path

      q_hvc = findPrime(beta_hvc)
      width = ZZ(ceil(log(q_hvc, hvc_arity)))
      width_too_big = bool(2^log_width < width)

    (sis_is_hard, error_msg) = rsisIsHard(beta_hvc, q_hvc, n, 2 * width, c)

    if sis_is_hard:
      fail_prob_path = - fail_prob_target -1 

      hvc_path_size = (ceil(log(beta_path*2+1,2))*n*width*2*tau)/8/1024
      if hvc_path_size < hvc_min_size:
        hvc_min_params = hvc_arity
        hvc_min_size = hvc_path_size
        hvc_params[hvc_arity] = {"arity" : hvc_arity,"q" : q_hvc, "norm bound" : beta_path, "SIS norm bound" : beta_hvc, "SIS width" : width, "size" : hvc_path_size,"failure prob" : fail_prob_path}

  to_tabulate = []
  for p in hvc_params.values():
    if p["arity"] == hvc_min_params:
      to_tabulate.append([Back.YELLOW+str(p["arity"]),str(p["q"]),str(p["norm bound"]),str(p["SIS norm bound"]),str(p["SIS width"]),("%.4f" % p["size"]) + " KB","2^"+str(p["failure prob"])+Back.RESET])
    else:
      to_tabulate.append([str(p["arity"]),str(p["q"]),str(p["norm bound"]),str(p["SIS norm bound"]),str(p["SIS width"]),("%.4f" % p["size"]) + " KB","2^"+str(p["failure prob"])])
  #  p[5] = ("%.4f" % p[5]) + " KB"
  if verbose:
    print(tabulate(to_tabulate,headers=["arity","q_HVC","beta_agg","beta_HVC","width","size","fail prob"]))

  f = open("summary.txt", "a")
  f.write("security param: " + str(secpar) +  ", rho: " + str(rho) + ", tau: " + str(tau) + "\n")
  f.write(str(hvc_params[hvc_min_params]) + "\n")
  f.close()

  f = open("hvc_sec" + str(secpar) + "_rho" + str(rho) + "_tau" + str(tau) + ".log", "w")
  f.write(tabulate(to_tabulate,headers=["arity","q_HVC","beta_agg","beta_HVC","width","size","fail prob"]))
  f.close()
  return hvc_params[hvc_min_params]


def find_param_hots_mcdiarmid(n, secpar, rho, tau, alpha_w, fail_prob_target, hvc, verbose):

  hvc_arity = hvc["arity"]
  hvc_beta_agg = hvc["norm bound"]
  bitsize_coeff_beta_agg = ceil(log(2*hvc_beta_agg+1,2))
  hvc_path_size = hvc["size"]
  hvc_width = hvc["SIS width"]


  # root hermite factor
  c = getRootHermiteFactor(secpar)

  # Hamming weight for hash function outputs
  alpha_H = findHammingWeight(n, 2^(2*secpar))
  print(alpha_H)

  # The range of the hash function is at most 2^delta times larger than needed. Usually delta=1.
  delta = ZZ(ceil(log(cardinalitySetOfTernaryPoly(n,alpha_H),2))-2*secpar)
  
  
  sig_params = {}
  sig_min_params = 0
  sig_min_size = oo

  sis_is_hard = True
  phi = 0
  while sis_is_hard:
    phi += 1
    beta_fresh = 2 * phi * alpha_H
    log_gamma = 2
    gamma_too_big = True
    while gamma_too_big:
      epsilon = fail_prob_target + 1 + log(n,2) + log_gamma
      beta_sigma = ceil(beta_fresh*sqrt((epsilon + 1) * log(2) * 2 * rho * alpha_w))
      beta_KOTS = 2*beta_sigma+4*alpha_w*phi*alpha_H
      q = findPrime(max(beta_KOTS,8*alpha_w*alpha_H*phi))
      gamma = findGamma(secpar,delta,n,q,phi)
      gamma_too_big = bool(2^log_gamma < gamma)
      log_gamma = floor(log(gamma,2))+1
    

    (sig_check, sig_check_msg) = rsisIsHard(beta_KOTS, q, n, gamma, c)
    (key_check, key_check_msg) = rsisIsHard(4*hvc_beta_agg, hvc["q"], n, 2*ceil(log(q,hvc_arity)), c)

    if sig_check:
      if key_check:
        sig_fail_prob = -fail_prob_target -1 
        hvc_beta_fresh = int((hvc_arity-1)/2)
        
        hvc_epsilon = fail_prob_target + 1 + ceil(log(2*n*ceil(log(q,hvc_arity))+2*n*hvc_width*tau,2))
        
        hvc_fail_prob = -(floor((hvc_beta_agg^2)/(2*log(2)*rho*alpha_w*hvc_beta_fresh^2) - 1) - ceil(log(2*n*ceil(log(q,hvc_arity))+2*n*hvc_width*tau,2)))
        total_fail_prob = max(sig_fail_prob,hvc_fail_prob)+1
        
        size = gamma*n*ceil(log(2*beta_sigma+1,2))/8/1024
        size_key = (ceil(log(2*hvc_beta_agg+1,2))*2*n*ceil(log(q,hvc_arity)))/8/1024
        total_size = size + size_key + hvc_path_size
        if total_size < sig_min_size:
          sig_min_params = phi
          sig_min_size = total_size
          sig_params[phi] = {"phi" : phi, "q" : q, "norm bound" : beta_sigma, "SIS norm bound" : beta_KOTS, "SIS width" : gamma, "sig size" : size, "key size" : size_key, "total size" : total_size, "sig failure prob" : sig_fail_prob, "hvc failure prob" : hvc_fail_prob, "total failure prob" : total_fail_prob}
    else:
      break

  to_tabulate = []
  for p in sig_params.values():
    if p["phi"] == sig_min_params:
      to_tabulate.append([Back.YELLOW+str(p["phi"]),str(p["q"]),str(p["norm bound"]),str(p["SIS norm bound"]),str(p["SIS width"]),("%.4f" % p["sig size"]) + " KB",("%.4f" % p["key size"]) + " KB",("%.4f" % p["total size"]) + " KB","2^"+str(p["sig failure prob"]),"2^"+str(p["hvc failure prob"]) ,"2^"+str(p["total failure prob"])+Back.RESET])
    else:
      to_tabulate.append([str(p["phi"]),str(p["q"]),str(p["norm bound"]),str(p["SIS norm bound"]),str(p["SIS width"]),("%.4f" % p["sig size"]) + " KB",("%.4f" % p["key size"]) + " KB",("%.4f" % p["total size"]) + " KB" ,"2^"+str(p["sig failure prob"]),"2^"+str(p["hvc failure prob"]), "2^"+str(p["total failure prob"])])

  if verbose:
    print(tabulate(to_tabulate,headers=["phi","q_KOTS","beta_sigma","beta_KOTS","gamma","signature size", "key size", "total size","sig fail prob","hvc fail prob","total fail prob"]))
    print("\n")

  f = open("summary.txt", "a")
  f.write(str(sig_params[sig_min_params]) + "\n\n")
  f.close()

  f = open("hots_sec" + str(secpar) + "_rho" + str(rho) + "_tau" + str(tau) + ".log", "w")
  f.write(tabulate(to_tabulate,headers=["phi","q_KOTS","beta_sigma","beta_KOTS","gamma","signature size", "key size", "total size","total fail prob"]))
  f.close()

def find_param(n, secpar, rho, tau, fail_prob_target, hots_key_size_hedge, verbose):
  print("finding param for sec =", secpar, " rho =", rho, " and tau =", tau)

  # Hamming weight for random weights
  alpha_w = findHammingWeight(n,2^secpar)
  print("alpha_w=", alpha_w)

  hvc = find_param_hvc_mcdiarmid(n, secpar, rho, tau, alpha_w, fail_prob_target+hots_key_size_hedge, verbose)

  find_param_hots_mcdiarmid(n, secpar, rho, tau, alpha_w, fail_prob_target, hvc, verbose)


# security parameter
secpars = [128]
# polynomial degree
n = 512
# number of users
rhos = [8192]
# height of the tree
taus = [21]

verbose = true

fail_prob_target = 14
hots_key_size_hedge = 1

for secpar in secpars:
  for rho in rhos:
    for tau in taus:
#      cProfile.run('find_param(n, secpar, rho, tau, fail_prob_target, hots_key_size_hedge, verbose)',sort='cumulative')
      find_param(n, secpar, rho, tau, fail_prob_target, hots_key_size_hedge, verbose)

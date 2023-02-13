from tabulate import tabulate
from colorama import init, Fore, Back, Style
import cProfile
load("dist.sage")

def cardinalitySetOfTernaryPoly(n,alpha):
  """ Determines the size of the set of ternary polynomials of degree n and Hamming weight alpha.
  """
  
  return binomial(n,alpha)*2^alpha

def findHammingWeight(n,l):
  """ Finds the minimal value alpha, such that the set of ternary polynomials of degree n and Hamming weight alpha has size at least l.
  """
  
  for alpha in range(n+1):
    if cardinalitySetOfTernaryPoly(n,alpha) >= l:
      return alpha
  raise ValueError("There does not exist a Hamming weight satisfying the specified conditions.")

def findGamma(secpar,delta,n,q,phi):
  """ Finds the minimal value gamma, such that the condition of Lemma TODO is satisfied.
  
  Lemma TODO requires that gamma is a positive integer such that 
    2^((3*secpar+deta)/(n*gamma))*q^(1/gamma) <= phi+1/2.
  This function finds the minimal value of gamma, such that the condition is 
  satisfied by solving the above inequality for gamma.
  """
  
  return ceil((((3*secpar+delta)/n)+log(q,2))/log(phi+.5,2))

def findNTTFriendlyPrime(n,beta):
  """ Finds the smallest prime q > beta, such that Z[x]/(x^n + 1) is NTT friendly.
  """
  
  q = next_prime(beta)
  while q % (2*n) != 1:
    q = next_prime(q)
  return q

def rsisIsHard(beta, q, n, m, c):
  """ Determines whether the Ring-SIS instance described by the inputs is hard.
  
  :param beta: infinity norm of the short solution
  :param q: modulus
  :param n: dimension of the polynomial ring
  :param m: number of ring elements
  :param c: root Hermite factor
  """

  # check that SIS problem is not trivial
  if beta >= q/2:
    return (false, "SIS is trivial")

  # find the best m for the attacker
  # if there are more m than required, the attacker can ignore some of them
  m = min(m, ceil(sqrt(n * log(q, 2) * log(c, 2))))

  # check that it is not possible to find short solution, using infinity norm
  if 2 * beta >= c^(n * m) * q^(1/m) - 1:
    return (false, "infinity norm fails")

  # check that it is not possible to find short solution, using l2 norm
  if sqrt(n * m) * beta >= c^(n * m) * sqrt( n * m / 2 / pi / e) * q^(1/m):
    return (false, "l2 norm fails")

  return (true, "pass")

def getRootHermiteFactor(secpar):
  """ Determines the handwavy root Hermite factor for the given security level
  """
  if secpar == 128:
    return 1.004
  else:
    if secpar == 112:
      return 1.005
    else:
      raise ValueError("Input security parameter should be either 112 or 128")

def find_kots_params(n, secpar, rho, alpha_w, fail_prob_target, verbose):
  """ Finds parameters for the key homomorphic one-time signature scheme compatible with the inputs.
  
  Specifically, the parameters should result in a scheme with secpar bits security that supports aggregation of up to rho signatures.
  When aggregating using uniformly random ternary polynomials with Hamming weight alpha_w, the aggregated signature will verify with probability at least 1-2^(-fail_prob_target).
  """
  
  # root hermite factor
  c = getRootHermiteFactor(secpar)

  # Hamming weight for hash function outputs
  alpha_H = findHammingWeight(n, 2^(2*secpar))

  # The range of the hash function is at most 2^delta times larger than needed. Usually delta=1.
  delta = ZZ(ceil(log(cardinalitySetOfTernaryPoly(n,alpha_H),2))-2*secpar)
  
  params = {}
  min_params = 0
  min_size = oo

  sis_is_hard = True
  phi = 1
  while sis_is_hard:
    # The norm bound for a freshly generated OTS (Lemma TODO)
    beta_fresh = 2 * phi * alpha_H
    # We have a circular dependecy between gamma, beta_sigma, and q. To solve 
    # this, we guess log_2(gamma) and check afterwards if the guess was correct.
    guessed_gamma = 1
    gamma_too_big = True
    while gamma_too_big:
      # The proof requires a union bound over gamma*n coefficients. To achieve 
      # overall failure probability below 2^(-fail_prob_target) we thus require 
      # a per-coefficient failure probability below 
      #   (2^{-fail_prob_target})/(gamma*n) 
      # = (2^{-(fail_prob_target + log_2(gamma) + log_2(n))}
      # (See Lemma TODO)
      epsilon = fail_prob_target + log(n*guessed_gamma,2)
      # Use McDiarmid bound to choose a norm bound that a single coefficient 
      # will exceed with probability at most 2^{-epsilon}
      beta_sigma = ceil(beta_fresh*sqrt((epsilon + 1) * log(2) * 2 * rho * alpha_w))
      # The norm bound of the corresponding SIS-instance (Theorem TODO)
      beta_kots = 2*beta_sigma+4*alpha_w*phi*alpha_H
      # Find a large enough NTT friendly prime q. The first bound is required 
      # for SIS-Hardness, the second one for Lemma TODO 19 of the squirrel paper
      q = findNTTFriendlyPrime(n,max(2*beta_kots,16*alpha_w*alpha_H*phi))
      # Choose gamma according to the conditions in Lemma TODO
      gamma = findGamma(secpar,delta,n,q,phi)
      # Check if we guessed gamma correctly
      gamma_too_big = bool(guessed_gamma < gamma)
      guessed_gamma +=1
    
    #Check whether the resulting SIS-Instance is hard
    (sis_is_hard, sis_check_msg) = rsisIsHard(beta_kots, q, n, gamma, c)
    if sis_is_hard:
      # Determine tighter bound on failure probability by inverting McDiarmid
      fail_prob = -(beta_sigma^2 / (beta_fresh^2 * log(2) * 2 * rho * alpha_w) - 1) + log(n*gamma,2)
      # A Signature consists of gamma ring elements with n coefficients each.
      # Each coefficients is bounded by beta_sigma, requiring 
      # log_2(2*beta_sigma+1) bits to store.
      size = (gamma * n) * ceil(log(2 * beta_sigma + 1,2))/8/1024
      if size < min_size:
        min_params = phi
        min_size = size
        params[phi] = {"alpha_H" : alpha_H, "delta" : delta, "phi" : phi, "gamma" : gamma, "beta_sigma" : beta_sigma, "q" : q, "SIS beta" : beta_kots, "SIS width" : 2*gamma, "size" : size, "failure prob" : fail_prob}
    phi += 1
    
  to_tabulate = []
  if verbose:
    for p in params.values():
      if p["phi"] == min_params:
        to_tabulate.append([Back.YELLOW+str(p["alpha_H"]),str(p["delta"]), str(p["phi"]),str(p["gamma"]),str(p["beta_sigma"]),str(p["q"]),str(p["SIS beta"]),str(p["SIS width"]),("%.4f" % p["size"]) + " KB","2^"+str(p["failure prob"])+Back.RESET])
      else:
        to_tabulate.append([str(p["alpha_H"]),str(p["delta"]), str(p["phi"]),str(p["gamma"]),str(p["beta_sigma"]),str(p["q"]),str(p["SIS beta"]),str(p["SIS width"]),("%.4f" % p["size"]) + " KB","2^"+str(p["failure prob"])])
    print(tabulate(to_tabulate,headers=["alpha_H", "delta", "phi","gamma","beta_sigma","q_kots","beta_kots","SIS width","signature size","sig fail prob"]),"\n")

#  f = open("summary.txt", "a")
#  f.write(str(params[min_params]) + "\n\n")
#  f.close()

#  f = open("kots_sec" + str(secpar) + "_rho" + str(rho) + "_tau" + str(tau) + ".log", "w")
#  f.write(tabulate(to_tabulate,headers=["alpha_w", "phi","gamma","beta_sigma","q_kots","beta_kots","SIS width","signature size","sig fail prob"]))
#  f.close()
  return params[min_params]

def find_hvc_params(n, secpar, rho, tau, alpha_w, kots, fail_prob_target, verbose):
  """ Finds parameters for the homomorphic vector commitment compatible with the inputs.
  
  Specifically, the parameters should result in a vector commitment with secpar bits security that supports vectors of length 2^tau and aggregation of up to rho openings.
  When aggregating using uniformly random ternary polynomials with Hamming weight alpha_w, the aggregated opening will verify correctly with probability at least 1-2^(-fail_prob_target).
  """
  
  # root hermite factor
  c = getRootHermiteFactor(secpar)

  hvc_params = {}
  hvc_min_params = 0
  hvc_min_size = oo

  q_kots = kots["q"]
  
  hvc_arity = 1
  sis_is_hard = True
  while sis_is_hard:
    hvc_arity+=2
    # The keys we're hashing into the leaves consist of 2 R_{q_kots} elements 
    # decompsed into ceil(log_arity(q_kots)) R_q elements.
    payload_width = 2*ceil(log(q_kots,hvc_arity))
    guessed_width = 1
    width_too_big = True
    beta_fresh = ZZ((hvc_arity-1)/2)
    while width_too_big:
      guessed_width+=1
      epsilon = fail_prob_target + log(n*(payload_width+2*guessed_width*tau),2)
      beta_path = ceil(beta_fresh*sqrt((epsilon + 1) * log(2) * 2 * rho * alpha_w))
      beta_hvc = 4 * beta_path

      q_hvc = findNTTFriendlyPrime(n,2*beta_hvc)
      width = ZZ(ceil(log(q_hvc, hvc_arity)))
      width_too_big = bool(guessed_width < width)

    (path_sis_is_hard, path_error_msg) = rsisIsHard(beta_hvc, q_hvc, n, 2 * width, c)
    (payload_sis_is_hard, payload_error_msg) = rsisIsHard(beta_hvc, q_hvc, n, payload_width, c)
    sis_is_hard = (path_sis_is_hard and payload_sis_is_hard)
    if sis_is_hard:
      fail_prob = -(beta_path^2 / (beta_fresh^2 * log(2) * 2 * rho * alpha_w) - 1) + log(n*(payload_width+2*width*tau),2)

      path_size = (ceil(log(beta_path*2+1,2))*n*width*2*tau)/8/1024
      payload_size = (ceil(log(beta_path*2+1,2))*payload_width*n)/8/1024
      if path_size+payload_size < hvc_min_size:
        hvc_min_params = hvc_arity
        hvc_min_size = path_size+payload_size
        hvc_params[hvc_arity] = {"arity" : hvc_arity, "width" : width, "payload width" : payload_width ,"beta" : beta_path, "q" : q_hvc, "SIS beta" : beta_hvc, "SIS width" : 2*width, "path size" : path_size, "payload size" : payload_size, "failure prob" : fail_prob}

  if verbose:
    to_tabulate = []
    for p in hvc_params.values():
      if p["arity"] == hvc_min_params:
        to_tabulate.append([Back.YELLOW+str(p["arity"]),str(p["width"]),str(p["payload width"]),str(p["beta"]),str(p["q"]),str(p["SIS beta"]),str(p["SIS width"]),("%.4f" % p["path size"]) + " KB",("%.4f" % p["payload size"]) + " KB",("%.4f" % (p["path size"]+p["payload size"])) + " KB","2^"+str(p["failure prob"])+Back.RESET])
      else:
        to_tabulate.append([str(p["arity"]),str(p["width"]),str(p["payload width"]),str(p["beta"]),str(p["q"]),str(p["SIS beta"]),str(p["SIS width"]),("%.4f" % p["path size"]) + " KB",("%.4f" % p["payload size"]) + " KB",("%.4f" % (p["path size"]+p["payload size"])) + " KB","2^"+str(p["failure prob"])])
    print(tabulate(to_tabulate,headers=["arity","width","payload width","beta_agg","q_HVC","beta_hvc","SIS width","path size","payload size","total size","fail prob"]))

#  f = open("summary.txt", "a")
#  f.write("security param: " + str(secpar) +  ", rho: " + str(rho) + ", tau: " + str(tau) + "\n")
#  f.write(str(hvc_params[hvc_min_params]) + "\n")
#  f.close()

#  f = open("hvc_sec" + str(secpar) + "_rho" + str(rho) + "_tau" + str(tau) + ".log", "w")
#  f.write(tabulate(to_tabulate,headers=["arity","q_HVC","beta_agg","beta_HVC","width","size","fail prob"]))
#  f.close()
  return hvc_params[hvc_min_params]

def find_param(n, secpar, rho, tau, fail_prob_target, verbose):
  print("finding param for sec =", secpar, " rho =", rho, " and tau =", tau)

  # Find Hamming weight for weights in random linear combination
  alpha_w = findHammingWeight(n,2^secpar)

  # Find parameters for the key homomorphic one-time signature scheme compatible with the given constraints.
  kots_param = find_kots_params(n, secpar, rho, alpha_w, fail_prob_target+1, verbose)
  # Find parameters for the homomorphic vector commitment compatible with the given constraints and the KOTS parameters.
  hvc_param = find_hvc_params(n, secpar, rho, tau, alpha_w, kots_param, fail_prob_target+1, verbose)
  return (alpha_w,kots_param,hvc_param)
  
def find_params(n,secpars,taus,rhos,fail_prob_target,verbosity):
  params = {}
  to_tabulate = []
  for secpar in secpars:
    params[secpar] = {}
    for tau in taus:
      params[secpar][tau] = {}
      for rho in rhos:
#      cProfile.run('find_param(n, secpar, rho, tau, fail_prob_target, verbose)',sort='cumulative')
        params[secpar][tau][rho] = find_param(n, secpar, rho, tau, fail_prob_target, bool(verbosity>1))
        to_tabulate.append([str(secpar),str(tau),str(rho),str(params[secpar][tau][rho][0]),str(params[secpar][tau][rho][1]["alpha_H"]),str(params[secpar][tau][rho][1]["delta"]),str(params[secpar][tau][rho][1]["phi"]),str(params[secpar][tau][rho][1]["gamma"]),str(params[secpar][tau][rho][1]["beta_sigma"]),str(params[secpar][tau][rho][1]["q"]),str(params[secpar][tau][rho][2]["arity"]),str(params[secpar][tau][rho][2]["beta"]),str(params[secpar][tau][rho][2]["q"]),("%.4f" % (params[secpar][tau][rho][1]["size"]+params[secpar][tau][rho][2]["path size"]+params[secpar][tau][rho][2]["payload size"]))+" KB"])
  if verbosity > 0:
    print(tabulate(to_tabulate,headers=["secpar","tau","rho","alpha_w","alpha_H","delta","phi","gamma","beta_sigma","q_kots","eta","beta_open","q_hvc","size"],tablefmt="simple_outline"))
  f = open("params_table.tex","w")
  f.write(tabulate(to_tabulate,headers=["$\\secpar$","$\\tau$","$\\rho$","$\\alpha_w$","$\\alpha_H$","$\\delta$","$\\varphi$","$\\gamma$","$\\beta_\\sigma$","$\\qkots$","$\\eta$","$\\betaopen$","$\\qhvc$","Size"],tablefmt="latex_booktabs").replace("\\textbackslash{}","\\").replace("\\$","$").replace("\\_","_"))
  f.close()
  return params





# security parameter
secpars = [112,128]
# polynomial degree
n = 512
# number of users
rhos = [1024, 8192, 131072]
# height of the tree
taus = [21, 24, 26]

verbosity = 1

fail_prob_target = 20


find_params(n,secpars,taus,rhos,fail_prob_target,verbosity)

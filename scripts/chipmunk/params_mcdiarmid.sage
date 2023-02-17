from tabulate import tabulate

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
  
  return ZZ(ceil((((3*secpar+delta)/n)+log(q,2))/log(phi+.5,2)))

def findNTTFriendlyPrime(n,beta):
  """ Finds the smallest prime q > beta, such that Z[x]/(x^n + 1) is NTT friendly.
  """
  
  q = next_prime(beta)
  # NTT friendliness requires that q ≡ 1 (mod 2*n)
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
  """ Uses handwaving to translates the security parameter into a root Hermite factor.
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
    # this, we guess gamma and check afterwards if the guess was correct.
    guessed_gamma = 1
    gamma_too_big = True
    while gamma_too_big:
      guessed_gamma +=1
      # The proof requires a union bound over gamma*n coefficients. To achieve 
      # overall failure probability below 2^(-fail_prob_target) we thus require 
      # a per-coefficient failure probability below 
      #   (2^{-fail_prob_target})/(gamma*n) 
      # = (2^{-(fail_prob_target + log_2(gamma) + log_2(n))}
      # (See Lemma TODO)
      epsilon = fail_prob_target + log(n*guessed_gamma,2)
      # Use McDiarmid bound to choose a norm bound that a single coefficient 
      # will exceed with probability at most 2^{-epsilon}
      beta_sigma = ZZ(ceil(beta_fresh*sqrt((epsilon + 1) * log(2) * 2 * rho * alpha_w)))
      # The norm bound of the corresponding SIS-instance (Theorem TODO)
      beta_kots = 2*beta_sigma+4*alpha_w*phi*alpha_H
      # Find a large enough NTT friendly prime q. The first bound is required 
      # for SIS to be non-trivial, the second one for Lemma TODO 19 of the squirrel paper
      q = findNTTFriendlyPrime(n,max(2*beta_kots,16*alpha_w*alpha_H*phi))
      # Choose gamma according to the conditions in Lemma TODO
      gamma = findGamma(secpar,delta,n,q,phi)
      # Check if we guessed gamma correctly
      gamma_too_big = bool(guessed_gamma < gamma)
    
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
      to_tabulate.append([str(p["alpha_H"]),str(p["delta"]), str(p["phi"]),str(p["gamma"]),str(p["beta_sigma"]),str(p["q"]),str(p["SIS beta"]),str(p["SIS width"]),("%.4f" % p["size"]) + " KB","2^"+("%.6f" % p["failure prob"])])
    print(tabulate(to_tabulate,headers=["alpha_H", "delta", "phi","gamma","beta_sigma","q_kots","beta_kots","SIS width","signature size","sig fail prob"],tablefmt="simple_outline"),"\n")
  return params[min_params]

def find_hvc_params(n, secpar, rho, tau, alpha_w, payload_width, q_payload, fail_prob_target, verbose):
  """ Finds parameters for the homomorphic vector commitment compatible with the inputs.
  
  Specifically, the parameters should result in a vector commitment with secpar bits security that supports vectors of length 2^tau of payloads consisting of payload_width R_{q_payload} elements and aggregation of up to rho openings.
  When aggregating using uniformly random ternary polynomials with Hamming weight alpha_w, the aggregated opening will verify correctly with probability at least 1-2^(-fail_prob_target).
  """
  
  # root hermite factor
  c = getRootHermiteFactor(secpar)

  hvc_params = {}
  hvc_min_params = 0
  hvc_min_size = oo

  hvc_arity = 1
  sis_is_hard = True
  while sis_is_hard:
    hvc_arity+=2
    # The payloads we're hashing into the leaves consist of payload_width many R_{q_payload} elements 
    # decompsed into ceil(log_arity(q_payload)) R_q elements each.
    decomposed_payload_width = payload_width*ceil(log(q_payload,hvc_arity))
    # We have a circular dependecy between width, beta_path, and q_hvc. To solve 
    # this, we guess width and check afterwards if the guess was correct.
    guessed_width = 1
    width_too_big = True
    beta_fresh = ZZ((hvc_arity-1)/2)
    while width_too_big:
      guessed_width+=1
      # The proof requires a union bound over all coefficients in all ring 
      # elements making up the payload and the path. The payload consists of
      # decomposed_payload_width many ring elements and the path consists of
      # 2*guessed_width*tau ring elements. Since each one has n coefficients, to
      # achieve an overall failure probability below 2^(-fail_prob_target) we 
      # require a per-coefficient failure probability below 
      #   (2^{-fail_prob_target})/(n*(decomposed_payload_width+2*guessed_width*tau)) 
      # = (2^{-(fail_prob_target + log_2(n*(decomposed_payload_width+2*guessed_width*tau)))}
      # (See Lemma TODO)
      epsilon = fail_prob_target + log(n*(decomposed_payload_width+2*guessed_width*tau),2)
      # Use McDiarmid bound to choose a norm bound that a single coefficient 
      # will exceed with probability at most 2^{-epsilon}
      beta_path = ceil(beta_fresh*sqrt((epsilon + 1) * log(2) * 2 * rho * alpha_w))
      # The norm bound of the corresponding SIS-instance (Theorem TODO)
      beta_hvc = 4 * beta_path
      # Find a large enough NTT friendly prime q_hvc. 
      # The bound is required for SIS to be non-trivial.
      q_hvc = findNTTFriendlyPrime(n,2*beta_hvc)
      # Compute actual width of a decomposed ring element.
      width = ZZ(ceil(log(q_hvc, hvc_arity)))
      # Check if we guessed width correctly.
      width_too_big = bool(guessed_width < width)
    
    # Check whether the resulting two SIS-instances are hard. The first one is 
    # used to hash along the path. the second one is used to hash the payload 
    # into the leaves.
    (path_sis_is_hard, path_error_msg) = rsisIsHard(beta_hvc, q_hvc, n, 2 * width, c)
    (payload_sis_is_hard, payload_error_msg) = rsisIsHard(beta_hvc, q_hvc, n, decomposed_payload_width, c)
    sis_is_hard = (path_sis_is_hard and payload_sis_is_hard)
    if sis_is_hard:
      # Determine tighter bound on failure probability by inverting McDiarmid
      fail_prob = -(beta_path^2 / (beta_fresh^2 * log(2) * 2 * rho * alpha_w) - 1) + log(n*(decomposed_payload_width+2*width*tau),2)
      # A path consists of 2*width*tau ring elements with n coefficients each.
      # Each coefficient is bounded by beta_path, requiring 
      # log_2(2*beta_path+1) bits to store.
      path_size = (ceil(log(2*beta_path+1,2))*n*2*width*tau)/8/1024
      # A payload consists of decomposed_payload_width ring elements with n 
      # coefficients each. Each coefficient is bounded by beta_path, requiring 
      # log_2(2*beta_path+1) bits to store.
      payload_size = (ceil(log(beta_path*2+1,2))*decomposed_payload_width*n)/8/1024
      if path_size+payload_size < hvc_min_size:
        hvc_min_params = hvc_arity
        hvc_min_size = path_size+payload_size
        hvc_params[hvc_arity] = {"arity" : hvc_arity, "width" : width, "payload width" : decomposed_payload_width ,"beta" : beta_path, "q" : q_hvc, "SIS beta" : beta_hvc, "SIS width" : 2*width, "path size" : path_size, "payload size" : payload_size, "failure prob" : fail_prob}

  if verbose:
    to_tabulate = []
    for p in hvc_params.values():
      to_tabulate.append([str(p["arity"]),str(p["width"]),str(p["payload width"]),str(p["beta"]),str(p["q"]),str(p["SIS beta"]),str(p["SIS width"]),("%.4f" % p["path size"]) + " KB",("%.4f" % p["payload size"]) + " KB",("%.4f" % (p["path size"]+p["payload size"])) + " KB","2^"+("%.6f" % p["failure prob"])])
    print(tabulate(to_tabulate,headers=["arity","width","payload width","beta_agg","q_HVC","beta_hvc","SIS width","path size","payload size","total size","fail prob"],tablefmt="simple_outline"),"\n")
  return hvc_params[hvc_min_params]

def find_param(n, secpar, rho, tau, fail_prob_target, verbose):
  print("Finding param for sec = " + str(secpar) + " rho = " + str(rho) + " tau = " + str(tau) + ", and failure probability 2^-" + str(fail_prob_target))

  # Find Hamming weight for weights in random linear combination
  alpha_w = findHammingWeight(n,2^secpar)

  # Find parameters for the key homomorphic one-time signature scheme compatible with the given constraints.
  kots_param = find_kots_params(n, secpar, rho, alpha_w, fail_prob_target+1, verbose)
  # Find parameters for the homomorphic vector commitment compatible with the given constraints and the KOTS parameters.
  hvc_param = find_hvc_params(n, secpar, rho, tau, alpha_w, 2, kots_param["q"], fail_prob_target+1, verbose)
  return (alpha_w,kots_param,hvc_param)
  
def find_params(n,secpars,taus,rhos,fail_prob_targets,verbosity):
  params = {}
  to_tabulate = []
  for secpar in secpars:
    params[secpar] = {}
    for tau in taus:
      params[secpar][tau] = {}
      for rho in rhos:
        params[secpar][tau][rho] = {}
        for fail_prob_target in fail_prob_targets:
#        cProfile.run('find_param(n, secpar, rho, tau, fail_prob_target, verbose)',sort='cumulative')
         params[secpar][tau][rho][fail_prob_target] = find_param(n, secpar, rho, tau, fail_prob_target, bool(verbosity>1))
         to_tabulate.append([
           str(secpar),
           str(tau),
           str(rho),
           str(fail_prob_target),
           str(params[secpar][tau][rho][fail_prob_target][0]),
           str(params[secpar][tau][rho][fail_prob_target][1]["alpha_H"]),
           str(params[secpar][tau][rho][fail_prob_target][1]["delta"]),
           str(params[secpar][tau][rho][fail_prob_target][1]["phi"]),
           str(params[secpar][tau][rho][fail_prob_target][1]["gamma"]),
           str(params[secpar][tau][rho][fail_prob_target][1]["beta_sigma"]),
           str(params[secpar][tau][rho][fail_prob_target][1]["q"]),
           str(params[secpar][tau][rho][fail_prob_target][2]["arity"]),
           str(params[secpar][tau][rho][fail_prob_target][2]["beta"]),
           str(params[secpar][tau][rho][fail_prob_target][2]["q"]),
           ("%.4f" % (params[secpar][tau][rho][fail_prob_target][1]["size"]+params[secpar][tau][rho][fail_prob_target][2]["path size"]+params[secpar][tau][rho][fail_prob_target][2]["payload size"]))+" KB"
         ])
  if verbosity > 0:
    print(tabulate(to_tabulate,headers=[
        "secpar",
        "tau",
        "rho",
        "epsilon",
        "alpha_w",
        "alpha_H",
        "delta",
        "phi",
        "gamma",
        "beta_sigma",
        "q_kots",
        "eta",
        "beta_open",
        "q_hvc",
        "size"
      ],tablefmt="simple_outline"))
  return params

# All constraints referenced within this script refer to the numbered constraints found in Table 3 of the Chipmunk paper.
from tabulate import tabulate

def cardinality_of_set_of_ternary_poly(n,alpha):
  """ Determines the size of the set of ternary polynomials of degree n and Hamming weight alpha.
  """
  
  return binomial(n,alpha)*2^alpha

def find_hamming_weight(n,l):
  """ Finds the minimal value alpha, such that the set of ternary polynomials of degree n and Hamming weight alpha has size at least l.
  """
  
  for alpha in range(n+1):
    if cardinality_of_set_of_ternary_poly(n,alpha) >= l:
      return alpha
  raise ValueError("There does not exist a Hamming weight satisfying the specified conditions.")

def get_gamma(secpar,delta,n,q,phi):
  """ Determines the minimal value gamma, such that comstraint 6 is satisfied.
  """
  
  return ZZ(ceil((((3*secpar+delta)/n)+log(q,2))/log(phi+.5,2)))
  
def get_alpha_w(secpar,n):
  """ Determines the minimal value alpha_w, such that constraint 11 is satisfied.
  """
  return find_hamming_weight(n,2^secpar)
  
def get_alpha_H_and_delta(secpar,n):
  """ Determines the minimal values alpha_H and delta, such that constraints 5 and 7 are satisfied.
  """
  alpha_H = find_hamming_weight(n, 2^(2*secpar))
  delta = ZZ(ceil(log(cardinality_of_set_of_ternary_poly(n,alpha_H),2))-2*secpar)
  return (alpha_H,delta)

def get_beta_sigma(n,alpha_H,alpha_w,rho,phi,gamma,epsilon):
  """ Determines the minimal value beta_sigma, such that constraint 4 is satisfied.
  """
  return ZZ(ceil(4*phi*alpha_H*sqrt(.5*alpha_w*rho*log(2*n*gamma/epsilon))))

def get_beta_kots(alpha_w,alpha_H,phi,beta_sigma):
  """ Calculates the norm bound of the SIS instance corresponding to the given KOTS parameters.
  """
  return ZZ(2*beta_sigma+4*alpha_w*alpha_H*phi)

def get_beta_agg(n,tau,xi,eta,kappa,kappaprime,alpha_w,rho,epsilon):
  """ Determines the minimal value beta_agg, such that constraint 1 is satisfied.
  """
  return ZZ(ceil(eta*sqrt(2*alpha_w*rho*(log(2*n/epsilon)+log(2*tau*kappa+xi*kappaprime)))))


def find_ntt_friendly_prime(n,beta):
  """ Finds the smallest prime q > beta, such that Z[x]/(x^n + 1) is NTT friendly.
  """
  
  q = next_prime(beta)
  # NTT friendliness requires that q ≡ 1 (mod 4*n)
  while q % (4*n) != 1:
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

def get_root_hermite_factor(secpar):
  """ Uses handwaving to translates the security parameter into a root Hermite factor.
  """
  
  if secpar == 128:
    return 1.004
  else:
    if secpar == 112:
      return 1.005
    else:
      raise ValueError("Input security parameter should be either 112 or 128")

def find_kots_params(n, secpar, rho, alpha_w, epsilon, verbose):
  """ Finds parameters for the key homomorphic one-time signature scheme compatible with the inputs.
  
  Specifically, the parameters should result in a scheme with secpar bits security that supports aggregation of up to rho signatures.
  When aggregating using uniformly random ternary polynomials with Hamming weight alpha_w, the aggregated signature will verify with probability at least 1 - epsilon.
  If verbose is set to true, the function will print a table of all found parameter sets.
  """
  
  # root hermite factor
  c = get_root_hermite_factor(secpar)

  # Get alpha_H and delta satisfying constraints 5 and 7
  (alpha_H,delta) = get_alpha_H_and_delta(secpar,n)
  
  params = {}
  min_params = 0
  min_size = oo

  sis_is_hard = True
  # Larger values for phi lead to smaller signatures.
  # But larger values for phi also make the related SIS problem easier.
  # We therefore try increasingly large values for phi until SIS becomes easy.
  phi = 1
  while sis_is_hard:
    # We have a circular dependecy between gamma, beta_sigma, and q.
    # To solve this, we guess gamma and check if the guess was correct.
    guessed_gamma = 1
    gamma_too_big = True
    while gamma_too_big:
      guessed_gamma +=1
      # Get beta_sigma satisfying constraint 4
      beta_sigma = get_beta_sigma(n,alpha_H,alpha_w,rho,phi,guessed_gamma,epsilon)
      # The norm bound of the corresponding SIS-instance.
      beta_kots = get_beta_kots(alpha_w,alpha_H,phi,beta_sigma)
      # Find a large enough NTT friendly prime q. The first bound is required 
      # for SIS to be non-trivial, the second one to satisfy constraint 8.
      q = find_ntt_friendly_prime(n,max(2*beta_kots,16*alpha_w*alpha_H*phi))
      # Get gamma satisfying constraint 6
      gamma = get_gamma(secpar,delta,n,q,phi)
      # Check if we guessed gamma correctly
      gamma_too_big = bool(guessed_gamma < gamma)
    
    #Check whether the resulting SIS-Instance is hard
    (sis_is_hard, sis_check_msg) = rsisIsHard(beta_kots, q, n, gamma, c)
    if sis_is_hard:
      
      # A Signature consists of gamma ring elements with n coefficients each.
      # Each coefficients is bounded by beta_sigma, requiring 
      # log_2(2*beta_sigma+1) bits to store.
      size = (gamma * n) * ceil(log(2 * beta_sigma + 1,2))
      if size < min_size:
        min_params = phi
        min_size = size
        params[phi] = {"alpha_H" : alpha_H, "delta" : delta, "phi" : phi, "gamma" : gamma, "beta_sigma" : beta_sigma, "q'" : q, "size" : size}
    phi += 1
    
  to_tabulate = []
  if verbose:
    for p in params.values():
      to_tabulate.append([str(p["alpha_H"]),str(p["delta"]), str(p["phi"]),str(p["gamma"]),str(p["beta_sigma"]),str(p["q'"]),str(ZZ(ceil(p["size"]/8/1024))) + " KB"])
    print(tabulate(to_tabulate,headers=["alpha_H", "delta", "phi","gamma","beta_sigma","q'","signature size"],tablefmt="simple_outline"),"\n")
  return params[min_params]

def find_hvc_params(n, secpar, rho, tau, alpha_w, xi, qprime, epsilon, verbose):
  """ Finds parameters for the homomorphic vector commitment compatible with the inputs.
  
  Specifically, the parameters should result in a vector commitment with secpar bits security that supports vectors of length 2^tau of payloads consisting of xi R_{qprime} elements and aggregation of up to rho openings.
  When aggregating using uniformly random ternary polynomials with Hamming weight alpha_w, the aggregated opening will verify correctly with probability at least 1-epsilon.
  If verbose is set to true, the function will print a table of all found parameter sets to  standard output.
  """
  
  # root hermite factor
  c = get_root_hermite_factor(secpar)

  hvc_params = {}
  hvc_min_params = 0
  hvc_min_size = oo
  
  # Larger values for eta lead to smaller decommitments.
  # But larger values for eta also make the related SIS problem easier.
  # We therefore try increasingly large values for eta until SIS becomes easy.
  eta = 2
  sis_is_hard = True
  while sis_is_hard:
    kappaprime = ceil(log(qprime,2*eta+1))
    # We have a circular dependecy between kappa, beta_agg, and q. To solve 
    # this, we guess kappa and check if the guess was correct.
    guessed_kappa = 0
    kappa_too_big = True
    while kappa_too_big:
      guessed_kappa+=1
      beta_agg = get_beta_agg(n,tau,xi,eta,guessed_kappa,kappaprime,alpha_w,rho,epsilon)
      # The norm bound of the corresponding SIS-instance.
      beta_hvc = 4 * beta_agg
      # Find a large enough NTT friendly prime q_hvc. 
      # The bound is required for SIS to be non-trivial.
      q = find_ntt_friendly_prime(n,2*beta_hvc)
      # Compute actual kappa.
      kappa = ceil(log(q,2*eta+1))
      # Check if we guessed kappa correctly.
      kappa_too_big = bool(guessed_kappa < kappa)
    
    # Check whether the resulting two SIS-instances are hard. The first one is 
    # used to hash along the path. the second one is used to hash the payload 
    # into the leaves.
    (path_sis_is_hard, path_error_msg) = rsisIsHard(beta_hvc, q, n, 2 * kappa, c)
    (payload_sis_is_hard, payload_error_msg) = rsisIsHard(beta_hvc, q, n, xi*kappaprime, c)
    sis_is_hard = (path_sis_is_hard and payload_sis_is_hard)
    if sis_is_hard:
      # A path consists of 2*tau*kappa ring elements.
      # A payload consists of xi*kappa' ring elements.
      # Each ring elements consists of n coefficients each.
      # Each coefficient is bounded by beta_agg, requiring 
      # log_2(2*beta_agg+1) bits to store.
      size = (2*tau*kappa+xi*kappaprime)*n*ceil(log(2*beta_agg+1,2))
      if size < hvc_min_size:
        hvc_min_params = eta
        hvc_min_size = size
        hvc_params[eta] = {"eta" : eta, "kappa" : kappa, "kappa'" : kappaprime ,"beta_agg" : beta_agg, "q" : q, "SIS beta" : beta_hvc, "SIS width" : 2*kappa, "size" : size, "epsilon" : epsilon}
    eta+=1

  if verbose:
    to_tabulate = []
    for p in hvc_params.values():
      to_tabulate.append([str(p["eta"]),str(p["kappa"]),str(p["kappa'"]),str(p["beta_agg"]),str(p["q"]),str(p["SIS beta"]),str(p["SIS width"]),str(ZZ(ceil(p["size"]/8/1024))) + " KB",p["epsilon"]])
    print(tabulate(to_tabulate,headers=["eta","kappa","kappa'","beta_agg","q","beta_hvc","SIS width","opening size","epsilon"],tablefmt="simple_outline"),"\n")
  return hvc_params[hvc_min_params]

def find_param(n, secpar, rho, tau, epsilon, verbose):
  """
    Finds parameters for the Chipmunk multi-signature scheme compatible with the inputs.
    
    Specifically, the parameters should result in a synchronized multi-signature scheme with secpar bits security that supports 2^tau time periods and aggregation of up to rho signatures, where any individual aggregation attempt will fail with probability at most epsilon.
  """
  print("Finding params for secpar = " + str(secpar) + " tau = " + str(tau) + " rho = " + str(rho) + ", and epsilon=" + str(epsilon))

  # Find alpha_w satisfying constraint 11
  alpha_w = get_alpha_w(secpar,n)
  # Compute chi according to contraint 10
  chi = ZZ(ceil(secpar/log(1/(2*epsilon))))

  # Find parameters for the key homomorphic one-time signature scheme compatible with the given constraints.
  kots_param = find_kots_params(n, secpar, rho, alpha_w, epsilon, verbose)
  # Find parameters for the homomorphic vector commitment compatible with the given constraints and the KOTS parameters. Here xi is set to 2 satisfying constraint 9.
  hvc_param = find_hvc_params(n, secpar, rho, tau, alpha_w, 2, kots_param["q'"], epsilon, verbose)
  return (alpha_w,chi,kots_param,hvc_param)
  
def find_params(n,secpars,taus,rhos,epsilons,verbosity):
  """
    Returns a dictionary containing valid parameter sets for all combinations of input constraints.
    
    If verbose is set to 1, the function will print a table of the found parameter sets to standard output.
    If verbose is set to 2, additional information about the process of finding the parameter sets will also be printed to standard output.
  """
  params = {}
  to_tabulate = []
  for secpar in secpars:
    params[secpar] = {}
    for tau in taus:
      params[secpar][tau] = {}
      for rho in rhos:
        params[secpar][tau][rho] = {}
        for epsilon in epsilons:
         params[secpar][tau][rho][epsilon] = find_param(n, secpar, rho, tau, epsilon, bool(verbosity>1))
         to_tabulate.append([
           str(secpar),
           str(tau),
           str(rho),
           str(epsilon),
           str(params[secpar][tau][rho][epsilon][0]),
           str(params[secpar][tau][rho][epsilon][1]),
           str(params[secpar][tau][rho][epsilon][2]["alpha_H"]),
           str(params[secpar][tau][rho][epsilon][2]["phi"]),
           str(params[secpar][tau][rho][epsilon][2]["gamma"]),
           str(params[secpar][tau][rho][epsilon][2]["beta_sigma"]),
           str(params[secpar][tau][rho][epsilon][2]["q'"]),
           str(params[secpar][tau][rho][epsilon][3]["eta"]),
           str(params[secpar][tau][rho][epsilon][3]["beta_agg"]),
           str(params[secpar][tau][rho][epsilon][3]["q"]),
           str(round((params[secpar][tau][rho][epsilon][2]["size"]+params[secpar][tau][rho][epsilon][3]["size"]+ceil(log(params[secpar][tau][rho][epsilon][1],2)))/8/1024))+" KB"
         ])
  if verbosity > 0:
    print(tabulate(to_tabulate,headers=[
        "secpar",
        "tau",
        "rho",
        "epsilon",
        "alpha_w",
        "chi",
        "alpha_H",
        "phi",
        "gamma",
        "beta_sigma",
        "q'",
        "eta",
        "beta_agg",
        "q",
        "size"
      ],tablefmt="simple_outline"))
  return params
  
# security parameter
secpars = [112]
# polynomial degree
n = 512
# number of users
rhos = [4096]
# height of the tree
taus = [21,24,26]
# targeted failure probability
epsilons = [2^(-10),2^(-15),2^(-16)]

verbosity = 1

params = find_params(n,secpars,taus,rhos,epsilons,verbosity)

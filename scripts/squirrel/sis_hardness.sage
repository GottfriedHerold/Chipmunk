q = 0x80201
rho = 1000
alpha = 23
n = 256
logq = 20
dim = n*(2*logq+1)

search_space = RR(log(factorial(256) / factorial(256-alpha) / factorial(alpha), 2) + alpha)
print(search_space)

gaussian_heuristic = sqrt(dim/2/pi/e) * q^(1/(2*logq + 1))
target = 4*alpha*rho*sqrt(2*n*logq)

c = RR((target/gaussian_heuristic)^(1/dim))
print(c)
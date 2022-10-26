n = 256

acc = 0 


for alpha in range (25):
    acc += factorial(n)/factorial(n-alpha)/factorial(alpha) * 2^alpha

    print(alpha, RR(log(acc,2)))


alpha = 12
prb = factorial(n)/factorial(alpha)/factorial(alpha)/factorial(n-2*alpha) 
print(alpha, RR(log(prb,2)))
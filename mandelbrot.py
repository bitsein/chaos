import matplotlib.pyplot as plt
import numpy as np

def f(x,c):
    return x*x + c

z0 = 0
n = 500
a = int(2*n)
X = []
Y = []

for i in range(-a, a, 1):
    for j in range(-a, a, 1):
        z = z0
        c = (i/n) +(j/n)*1j
        for k in range(1000):
            z = f(z,c)
            if np.abs(z) > 1.5:
                break
        else:
            X.append(i/n)
            Y.append(j/n)
        
plt.scatter(X, Y, s = 0.5)
plt.xlabel('Re')
plt.ylabel('Im')
plt.savefig('man.eps')

import matplotlib.pyplot as plt
import numpy as np

s = 2.5
n = 15000
t = 0.0001

X = []
Y = []

def h(x,a):
    return a*x*(1-x)

x0 = np.random.rand() #初期値
for i in range(n):
    a = s + t * i
    x = x0
    for j in range(0, 10000):
        x = h(x,a) 
    for k in range(0, 20):
        x = h(x,a)
        X.append(a)
        Y.append(x)

plt.scatter(X, Y, s = 0.01)
plt.xlabel('a')
plt.ylabel('limit')
plt.savefig('logistic.eps')
plt.savefig('logistic.png')

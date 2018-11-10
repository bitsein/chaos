#Chua回路の４次ルンゲクッタ法によるシミュレーション
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys

args = sys.argv
argc = len(args)
if argc != 2:
    print('ERROR:Input_G_value')
    quit()

C_1 = 1/9
C_2 = 1
L = 1/7
G = float(args[1])
m_0 = -0.5
m_1 = -0.8
Bp = 1


def g(h):
    return m_0*h+(1/2)*(m_1-m_0)*abs(h+Bp)+(1/2)*(m_0-m_1)*abs(h-Bp)


vc_1 = []
vc_2 = []
iL = []

x = 0.01
y = 0.01
z = 0.01

vc_1.append(x)
vc_2.append(y)
iL.append(z)

t = 0.01
count = 10000

for var in range(0, count):

    a_1 = (G*(y-x)-g(x))/C_1
    b_1 = (G*(x-y)+z)/C_2
    c_1 = -y/L
    a_2 = (G*((y+b_1*t/2)-(x+a_1*t/2))-g((x+a_1*t/2)))/C_1
    b_2 = (G*((x+a_1*t/2)-(y+b_1*t/2))+(z+c_1*t/2))/C_2
    c_2 = -(y+b_1*t/2)/L
    a_3 = (G*((y+b_2*t/2)-(x+a_2*t/2))-g((x+a_2*t/2)))/C_1
    b_3 = (G*((x+a_2*t/2)-(y+b_2*t/2))+(z+c_2*t/2))/C_2
    c_3 = -(y+b_2*t/2)/L
    a_4 = (G*((y+b_3*t)-(x+a_3*t))-g((x+a_3*t)))/C_1
    b_4 = (G*((x+a_3*t)-(y+b_3*t))+(z+c_3*t))/C_2
    c_4 = -(y+b_3*t)/L

    x += (a_1+2*a_2+2*a_3+a_4)*t/6
    y += (b_1+2*b_2+2*b_3+b_4)*t/6
    z += (c_1+2*c_2+2*c_3+c_4)*t/6

    vc_1.append(x)
    vc_2.append(y)
    iL.append(z)

fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(vc_1, vc_2, iL, s=3)
ax.set_xlabel("vc_1")
ax.set_ylabel("vc_2")
ax.set_zlabel("iL")
ax.set_title(G)
'''plt.show()'''
plt.savefig('RungeKutta.eps')

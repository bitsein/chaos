#２重振り子（陽的オイラー法）
import numpy as np
import sys
import matplotlib.pyplot as plt

args = sys.argv
argc = len(args)
if argc != 5:
    print('ERROR:Input_four_initial_value.')
    quit()

g = 9.8
l1 = 0.5
l2 = 0.5
m1 = 0.5
m2 = 0.3

theta1 = float(args[1])
theta2 = float(args[2])
theta1dot = float(args[3])
theta2dot = float(args[4])
omega = np.sqrt(g / l1)
t = 0.01
count = 9999
T = np.arange(0, 100, 0.01)


def h(t1, t2, t1dot, t2dot):
    return 0.5 * m1 * l1 * l1 * t1dot * t1dot + 0.5 * m2 * (l1 * l1 * t1dot * t1dot + l2 * l2 * t2dot * t2dot + 2 * l1 * l2 * t1dot * t2dot * np.cos(t1 - t2)) \
           - m1 * l1 * g * np.cos(t1) - m2 * g * (l1 * np.cos(t1) + l2 * np.cos(t2))


M = m2 / (m1 + m2)
L = l2 / l1


def a(t1, t2, t1dot, t2dot):
    return (omega * omega * L * (-np.sin(t1) + M * np.cos(t1 - t2) * np.sin(t2)) - M * L * (t1dot * t1dot * np.cos(t1 - t2) + L * t2dot * t2dot) * np.sin(t1 - t2))\
            / (L - M * L * np.cos(t1 - t2) * np.cos(t1 - t2))


def b(t1, t2, t1dot, t2dot):
    return (omega * omega * np.cos(t1 - t2) * np.sin(t1) - omega * omega * np.sin(t2) + (t1dot * t1dot + M * L * t2dot * t2dot * np.cos(t1 - t2)) * np.sin(t1 - t2))\
            / (L - M * L * np.cos(t1 - t2) * np.cos(t1 - t2))


THETA1 = []
THETA2 = []
THETA1DOT = []
THETA2DOT = []
H = []

THETA1.append(theta1)
THETA2.append(theta2)
THETA1DOT.append(theta1dot)
THETA2DOT.append(theta2dot)
H.append(h(theta1, theta2, theta1dot, theta2dot))

x = theta1
y = theta2
z = theta1dot

for var in range(0, count):

    theta1 += t * theta1dot
    theta2 += t * theta2dot
    theta1dot += t * a(x, y, theta1dot, theta2dot)
    theta2dot += t * b(x, y, z, theta2dot)
    x = theta1
    y = theta2
    z = theta1dot
    THETA1.append(theta1)
    THETA2.append(theta2)
    THETA1DOT.append(theta1dot)
    THETA2DOT.append(theta2dot)
    H.append(h(theta1, theta2, theta1dot, theta2dot))

plt.xlabel("time")
plt.ylabel("theta2dot")
plt.plot(T, THETA2DOT)
plt.show()

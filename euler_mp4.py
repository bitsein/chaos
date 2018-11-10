#陽的オイラー法の結果をmp4形式で出力
#fpsとかの設定がいまいちわからなかった
import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.animation import FuncAnimation

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
count = 2500
T = np.arange(0, t*count, t)


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
X1 = []
Y1 = []
X2 = []
Y2 = []

THETA1.append(theta1)
THETA2.append(theta2)
THETA1DOT.append(theta1dot)
THETA2DOT.append(theta2dot)
H.append(h(theta1, theta2, theta1dot, theta2dot))
X1.append(l1 * np.sin(theta1))
Y1.append(- l1 * np.cos(theta1))
X2.append(l1 * np.sin(theta1) + l2 * np.sin(theta2))
Y2.append(- l1 * np.cos(theta1) - l2 * np.cos(theta2))


x = theta1
y = theta2
z = theta1dot

for var in range(0, count-1):

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
    X1.append(l1 * np.sin(theta1))
    Y1.append(- l1 * np.cos(theta1))
    X2.append(l1 * np.sin(theta1) + l2 * np.sin(theta2))
    Y2.append(- l1 * np.cos(theta1) - l2 * np.cos(theta2))



fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(-1.5, 1.5), ylim=(-1.5, 1.5))
ax.grid()

line, = plt.plot([], [], 'ro-', animated = True)

time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def init():
    time_text.set_text('')
    return line, time_text

def update(i):
    next_x = [0, X1[i], X2[i]]
    next_y = [0, Y1[i], Y2[i]]
    line.set_data(next_x, next_y)

    time_text.set_text(time_template % (i*t))
    return line, time_text

FRAME_INTERVAL = t*10 # [msec] interval between frames
FPS = 10 / FRAME_INTERVAL
ani = FuncAnimation(fig, update, frames=len(T),
                    interval=FRAME_INTERVAL, init_func=init, blit=True)
ani.save('single_pendulum.mp4', fps=FPS, extra_args=['-vcodec', 'libx264'])
plt.show()

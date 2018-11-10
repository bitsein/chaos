#symp.pyの結果をmp4形式にしたもの
import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.animation import FuncAnimation

args = sys.argv
argc = len(args)
if argc != 5:
    print('ERROR:Input_four_initial_value.')
    quit()


def newton( xn, h, l1, l2, l, m1, m2, m, q1, q2, p1, p2):
	F1 = l2 * p1 - l1 * p2 * np.cos(xn);
	F2 = (m1 * l1 + m2 * l1) * p2 - m2 * l2 * p1 * np.cos(xn);
	G = m1 + m2 - m2 * np.cos(xn) * np.cos(xn);
	F1dot = l1 * p2 * np.sin(xn);
	F2dot = m2 * l2 * p2 * np.sin(xn);
	Gdot = 2 * m2 * np.cos(xn) * np.sin(xn);
	A = h * F1 / (l1 * l1 *l2 * G);
	B = h * F2 / (m2 * l1 * l2 * l2 * G);
	fx = xn - (q1 + A - q2 - B);
	fdotx = 1 - h * ((F1dot * G - F1 * Gdot) / (l1 * l1 * l2 * G * G) - (F2dot * G - F2 * Gdot) / (m2 * l1 * l2 * l2 * G * G));
	return xn - fx / fdotx;


h = 0.01;
th1 = float(args[1]);
th2 = float(args[2]);
dth1 = float(args[3]);
dth2 = float(args[4]);
datasize=2000; #time = datasize * h

l1 = 0.5;
l2 = 0.5;
m1 = 0.5;
m2 = 0.3;
g = 9.8;
tol=pow(10, -12);
	
m = m2/(m1+m2);
l = l2/l1;
w = np.sqrt(g/l1);

X1 = []
Y1 = []
X2 = []
Y2 = []
LEN = np.arange(0, datasize * h, h)

X1.append(l1 * np.sin(th1))
Y1.append(- l1 * np.cos(th1))
X2.append(l1 * np.sin(th1) + l2 * np.sin(th2))
Y2.append(- l1 * np.cos(th1) - l2 * np.cos(th2))

q1=th1;
q2=th2;
dth = th1-th2;
p1=(m1+m2)*l1*l1*dth1+m2*l1*l2*np.cos(dth)*dth2;
p2=m2*l1*l2*np.cos(dth)*dth1+m2*l2*l2*dth2;
T=m1*l1*l1*dth1*dth1/2.0+m2*(l1*l1*dth1*dth1+l2*l2*dth2*dth2+2.0*l1*l2*dth1*dth2*np.cos(dth))/2.0;
U=-m1*l1*g*np.cos(th1)-m2*g*(l1*np.cos(th1)+l2*np.cos(th2));
Ham=T+U;
	
for var in range(0, datasize):
	bth1=th1;
	bth2=th2;
	bdth1=dth1;
	bdth2=dth2;
	xn=q1-q2;
	xn1=newton(xn,h,l1,l2,l,m1,m2,m,q1,q2,p1,p2);
	
	while abs(xn1-xn)>tol:
		xn=xn1;
		xn1=newton(xn,h,l1,l2,l,m1,m2,m,q1,q2,p1,p2);
	x=xn1;
	F1 = l2*p1-l1*p2*np.cos(x);
	F2 = (m1*l1 + m2*l1)*p2-m2*l2*p1*np.cos(x);
	G = m1 + m2 - m2*np.cos(x)*np.cos(x);
	fq1=q1+h*F1/(l1*l1*l2*G);
	fq2=q2+h*F2/(m2*l1*l2*l2*G);
	fp1=p1+h*(-F1*F2*np.sin(x)/(l1*l1*l2*l2*G*G)-(m1+m2)*l1*g*np.sin(fq1));
	fp2=p2+h*(F1*F2*np.sin(x)/(l1*l1*l2*l2*G*G)-m2*l2*g*np.sin(fq2));
	q1=fq1;
	q2=fq2;
	p1=fp1;
	p2=fp2;
	th1=q1;
	th2=q2;
	dth = th1-th2;
	M11=(m1+m2)*l1*l1;
	M12=m2*l1*l2*np.cos(dth);
	M21=m2*l1*l2*np.cos(dth);
	M22=m2*l2*l2;
	detM=M11*M22-M21*M12;
	invM11=M22/detM;
	invM12=-M12/detM;
	invM21=-M21/detM;
	invM22=M11/detM;
	dth1=invM11*p1+invM12*p2;
	dth2=invM21*p1+invM22*p2;
	T=m1*l1*l1*dth1*dth1/2.0+m2*(l1*l1*dth1*dth1+l2*l2*dth2*dth2+2.0*l1*l2*dth1*dth2*np.cos(dth))/2.0;
	U=-m1*l1*g*np.cos(th1)-m2*g*(l1*np.cos(th1)+l2*np.cos(th2));
	Ham=T+U;
	while th1<-np.pi:
		th1=th1+2.0*np.pi;
	while th1>np.pi:
		th1=th1-2.0*np.pi;
	while th2<-np.pi:
		th2=th2+2.0*np.pi;
	while th2>np.pi:
		th2=th2-2.0*np.pi;
	X1.append(l1 * np.sin(th1))
	Y1.append(- l1 * np.cos(th1))
	X2.append(l1 * np.sin(th1) + l2 * np.sin(th2))
	Y2.append(- l1 * np.cos(th1) - l2 * np.cos(th2))


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

    time_text.set_text(time_template % (i*h))
    return line, time_text

FRAME_INTERVAL = h*10 # [msec] interval between frames
FPS = 10 / FRAME_INTERVAL
ani = FuncAnimation(fig, update, frames=len(LEN),
                    interval=FRAME_INTERVAL, init_func=init, blit=True)
ani.save('DP.mp4', fps=FPS, extra_args=['-vcodec', 'libx264'])
plt.show()

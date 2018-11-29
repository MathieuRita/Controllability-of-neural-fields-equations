import numpy as np 
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.integrate as integrate


#Params
alpha=2
X0=0.
sigma=1.
lda=1.
X=np.linspace(-np.pi/2,np.pi/2,500)

#1. Developpement en serie de Fourier de la gaussienne en fonction de sigma

def gauss(X0=0,sigma=1):
	return 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (X - X0)**2 / (2 * sigma**2) )

def a(n,sigma=1,x0=0,T=np.pi,borne1=-np.pi/2,borne2=np.pi/2):
	def fonction(x):
		return  np.cos(n*((2*np.pi)/T)*x)* (1/(sigma * np.sqrt(2 * np.pi))) * (np.exp( - (x - x0)**2 / (2 * sigma**2)))
	inte=integrate.quad(fonction,borne1,borne2)
	if n==0:
		return (1/np.pi)*inte[0]
	else :
		return (2/np.pi)*inte[0] 

four=0*X 
for i in range(10):
	four+=a(i)*np.cos(2*i*X)

plt.plot(X,gauss())
plt.plot(X,four)
plt.show()

#2. Difference finie sur les equations de degre 3

#params simu
a0=a(0)
a1=a(1)
nbx=500
nbt=20
pasx=np.pi/nbx
T=1.
past=T/nbt
Vin=10.*gauss(0)+5.*gauss(np.pi/4,sigma=0.05)-8*gauss(np.pi/4.5,sigma=0.05)
#Vin=np.cos(X)
a=0.5
tau=past/a
print("tau"+str(tau))

def S(x):
	return 1/(1+np.exp(-x))

# Simulation
V=Vin*np.ones([nbt,nbx])

for t in range(1,nbt):
	for j in range(0,nbx):
		xj=-np.pi/2+(j*pasx)
		inte=0
		for i in range(nbx):
			xi=-np.pi/2+(i*pasx)
			inte+=(a0+a1*np.cos(alpha*(xi-xj))*S(lda*V[t-1,i]))
		V[t,j]=V[t-1,j]+(past/tau)*(-V[t-1,j]+(1/np.pi)*pasx*inte)
	print(t)

#Plots l'evolution des profils spatiaux
for time in range(0,nbt):
	plt.plot(X,V[time,:],label="t="+str(time))
plt.legend()
plt.show()

#Animation
fig, ax = plt.subplots()

line, = ax.plot(X, V[0,:])

def init():  # only required for blitting to give a clean slate.
    line.set_ydata([np.nan] * len(X))
    return line,


def animate(i):
    line.set_ydata(V[i,:])  # update the data.
    return line,


ani = animation.FuncAnimation(
    fig, animate, init_func=init, interval=1000, blit=True, save_count=None)

plt.show()

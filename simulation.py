import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.animation as animation

#Params
J0= -7.3
J1= 11
alpha=2
X0=0.
sigma=0.01
lda=1.

X=np.linspace(-np.pi/2,np.pi/2,500)
J=J0+J1*np.cos(alpha*X)

# plt.plot(X,J)
# plt.show()

a0=0.28
a1=0.063*2
a2=0
a3=0
#a2=-0.0072*2
#a3=0.0032*2

def gauss(X0,sigma=0.01):
	return 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (X - X0)**2 / (2 * sigma**2) )

four=a0+a1*np.cos(2*X)

 #V0(y) entre -pi/2 et pi/2

#params simu

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

V=Vin*np.ones([nbt,nbx])

#Simulations 1 des equations normales
for t in range(1,nbt):
	for j in range(0,nbx):
		xj=-np.pi/2+(j*pasx)
		inte=0
		for i in range(nbx):
			xi=-np.pi/2+(i*pasx)
			inte+=(a0+a1*np.cos(alpha*(xi-xj))*S(lda*V[t-1,i]))
		V[t,j]=V[t-1,j]+(past/tau)*(-V[t-1,j]+(1/np.pi)*pasx*inte)
	print(t)

for time in range(0,nbt):
	plt.plot(X,V[time,:],label="t="+str(time))
plt.legend()
plt.show()

plt.plot(np.linspace(0.,10,nbt),V[:,100])
plt.plot(np.linspace(0.,10,nbt),V[:,0])
plt.plot(np.linspace(0.,10,nbt),V[:,200])
plt.show()

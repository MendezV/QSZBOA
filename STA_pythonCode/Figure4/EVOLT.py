
###################################
#PACKAGES
###################################

import numpy as np
from scipy import linalg as la
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import sys

###################################

###################################





###############################################
##CONSTANTS
###############################################

#constants
cosa1=201
points_x=cosa1
points_t=4*int(int(points_x**2 /50)/4.0) #para quesea multiplo de 4
DD=-20
steps=int(cosa1*5)
hbar=1.0
m=1.0
W=10.0 #size of the well
d=0.5
cent=-0.0


# create grid
tf=1000
W=W/2.0
xf=W
x0=-W
dt=tf/points_t
dx=(xf-x0)/(points_x)
x=np.linspace(x0,xf,points_x) #debe tener un numero impar de elemntos para que sirva NInteg
t=np.linspace(0,tf,points_t)


#############################################

#############################################



############################################
#POTENTIALS
############################################


def window(xvec,d,cent):
    return 0.5*((np.sign(xvec+d/2.0-cent))-(np.sign(xvec-d/2.0-cent)))

def Vpike(A,xvec,d,cent):
    return A*((1-2*(xvec-cent)/d)*window(xvec,d/2,cent+d/4)+(1+2*(xvec-cent)/d)*window(xvec,d/2,cent-d/4))

def Vsq(A,xvec,d,cent):
    return A*window(xvec,d,cent)

def VJar(A1,A2,A3, xvec,cent):
    return A1*(xvec**4 - 0.5*(A2*(xvec)**2*(np.sign( (xvec-cent) )+1)+A3*(xvec)**2*(np.sign( -(xvec-cent) )+1)))

def Vharm(mw2,xvec,cent):
    return 0.5*mw2*(xvec-cent)**2

def Vsqdouble(A1,A3,d1,d3,xvec,cent1,cent2):
    return A1*window(xvec,d1,cent1)+A3*window(xvec,d3,cent2)





#############################################

#############################################






#############################################
##FUNCTIONS
#############################################

#########INTEGRATION########################
from scipy.interpolate import UnivariateSpline

def Ninteg(x,func,x0,xf,dx):
    spl = UnivariateSpline(x,func, k=3, s=0)
    derivativas=spl.derivative()
    intfunc=np.zeros(np.size(xf))
    for i in range(np.size(xf)):
        intfunc[i]=spl.integral(x0, xf[i])
 
    return intfunc

def Ninteg2(x,func,x0,xf,dx):
    spl = UnivariateSpline(x,func, k=3, s=0)
    intfunc=spl.integral(x0, xf)
    return intfunc

def eta(t):
    return ((t**3.0)/(tf**3.0))*(1 + 3.0*(1 - (t/tf)) + 6.0*(1 - (t/tf))**2.0)
#############################################


#########DIFFERENTIATION########################

def Nderivat(func,x):
    spl = UnivariateSpline(x,func, k=3, s=0)
    derivativas=spl.derivative()
    return derivativas(x)

def Nderivat2(func,x):
    Nderivat(func,x)
    return  Nderivat( Nderivat(func,x),x)

#############################################



#############################################

#############################################






#############################################
#INITIAL CONDITION
#############################################



psigrid=[]
norm=[]
onesies=np.diag(np.ones(points_x))
T=-0.5*(-2.0*np.diag(np.ones(points_x))+np.diag(np.ones(points_x-1),1)+np.diag(np.ones(points_x-1),-1))/(dx**2)+np.diag(Vsqdouble(DD,DD,2.5,2.5,x,-5+1.25,5-1.25))   #first hamiltonian

#T=-0.5*(-2.0*np.diag(np.ones(points_x))+np.diag(np.ones(points_x-1),1)+np.diag(np.ones(points_x-1),-1))/(dx**2)+np.diag( Vharm(1,x,0))   #first hamiltonian




T2=-0.5*(-2.0*np.diag(np.ones(points_x))+np.diag(np.ones(points_x-1),1)+np.diag(np.ones(points_x-1),-1))/(dx**2)+np.diag(Vsqdouble(0.05*DD,0.05*DD,2.5,2.5,x,-5+1.25,5-1.25))   #second hamiltonian


values2=la.eigh(T)
values=la.eigh(T2)

'''
psig=values2[1].T[0] #ground state for first potental
psie=values2[1].T[1] #excited state for first potental
Eg=values2[0][0]   #ground state energy for first potental
Ee=values2[0][1]   #excited state energy for first potental


psigprime=values[1].T[0] #ground state for first potental
psieprime=values[1].T[1] #excited state for first potental
Egprime=values[0][0] #ground state energy for first potental
Eeprime=values[0][1] #excited state energy for first potental
'''

##for 2-7

psig=(1/np.sqrt(2))*(values2[1].T[0]-values2[1].T[1]) #excited state for first potental #ground state for first potental
psie=values2[1].T[6] #excited state for first potental
Eg=values2[0][1]   #ground state energy for first potental
Ee=values2[0][6]   #excited state energy for first potental


psigprime=(1/np.sqrt(2))*(values2[1].T[0]-values2[1].T[1]) #excited state for first potental #ground state for first potental
psieprime=values[1].T[6] #excited state for first potental
Egprime=values[0][1] #ground state energy for first potental
Eeprime=values[0][6] #excited state energy for first potental



#interpolated states "prime"
psi1=psig/np.sqrt( Ninteg2(x,psig**2,x0,xf,dx))
psi2=psie/np.sqrt( Ninteg2(x,psie**2,x0,xf,dx))
E1=Eg
E2=Ee

print(Ninteg(x,psi1**2,x0,[xf],dx)[0])
print(Ninteg(x,psi2**2,x0,[xf],dx)[0])
'''
plt.plot(x,psi1)
plt.plot(x,psi2)
plt.show()
'''
#############################################

#############################################






###############################################################
filename=sys.argv[1]
V=np.loadtxt(filename,delimiter=',')
psigrid=[]
norm=[]
psigrid.append(psi1)
onesies=np.diag(np.ones(points_x))
T=-0.5*(-2.0*np.diag(np.ones(points_x))+np.diag(np.ones(points_x-1),1)
        +np.diag(np.ones(points_x-1),-1))/(dx**2)

bounds=1
plt.imshow(V, interpolation='nearest',aspect='auto')
plt.title(r'$\rho(x,t)$')
plt.ylabel(r'$t/t_{f}$')
plt.xlabel(r'$x$')
plt.xticks(np.arange(0,(np.shape(V)[1]-2*bounds+1),(np.shape(V)[1]-2*bounds)/(2*int(xf))),np.arange(-int(xf),2*int(xf)+1))
plt.yticks(np.arange(0,np.shape(V)[0],np.shape(V)[0]/4.0),np.linspace(0,1,5))
plt.colorbar()
plt.show()


for i in range(points_t):
    psipres=psigrid[i]
    #T=-0.5*Nderivat2(psipres,x)
    Vc=np.diag(V[i,:])
    Hpsi=T+Vc
    matt=onesies-1j*Hpsi*dt
    mattinv=la.inv( onesies +1j*Hpsi*dt)
    psinew=(mattinv@(matt@psipres))
    psigrid.append(psinew)

'''
for i in range(2*itertevol+points_t,2*itertevol+points_t+5*itertevol):
    psipres=psigrid[i]
    #T=-0.5*Nderivat2(psipres,x)
    Vc=np.diag(V[-1,:])
    Hpsi=T+Vc
    matt=onesies-1j*Hpsi*dt
    mattinv=la.inv( onesies +1j*Hpsi*dt)
    psinew=(mattinv@(matt@psipres))
    psigrid.append(psinew)

'''

import math
def polarThe(z):
    a= z.real
    b= z.imag
    theta = math.atan2(b,a)
    return theta 

def polarR(z):
    a= z.real
    b= z.imag
    r = math.hypot(a,b)
    return r

psirho=np.array([[polarR(psigrid[i][j]) for i in range(points_t)] for j in range(points_x)]).T
phiphase=np.array([[polarThe(psigrid[i][j]) for i in range(points_t)] for j in range(points_x)]).T

#psirho= np.sqrt(abs(   np.conj(np.array(psigrid))  *np.array(psigrid)   ))
#phiphase=np.real( 1j*np.log(   np.array(psigrid) / psirho   )  )




plt.plot(psirho[-1]**2)
plt.plot(psi2**2)
plt.title("final state norm")
plt.show()

print(Ninteg(x,psirho[0]**2,x0,[xf],dx)[0])
print(Ninteg(x,psirho[-1]**2,x0,[xf],dx)[0])

plt.plot(psirho[0]**2)
plt.title("Initial state norm")
plt.plot(psi1**2)
plt.show()


plt.plot(np.real(psigrid[-1]))
plt.plot(np.imag(psigrid[-1]))
plt.title("final state real and imag")
plt.show()

plt.plot(psigrid[-1])
plt.plot(psi2)
plt.title("final state and target")
plt.show()

plt.plot(phiphase[0])
plt.plot(phiphase[-1])
plt.title("phase")
plt.show()


'''

exp_energySTA=np.zeros(np.size(t)+2*itertevol)

for i in range(0,2*itertevol+points_t):

	Tcal=-psirho[i,:]*np.cos(phiphase[i,:])*Nderivat2(psirho[i,:]*np.sin(phiphase[i,:]),x)-psirho[i,:]*np.cos(phiphase[i,:])*Nderivat2(psirho[i,:]*np.sin(phiphase[i,:]),x)
	T_exp=Ninteg(x,Tcal,x0,[xf],dx)[0]
	Vcal=psirho[i,:]*psirho[i,:]*V[-1,:]
	V_exp=Ninteg(x,Vcal,x0,[xf],dx)[0]
	exp_energySTA[i]=T_exp+V_exp


for i in range(0,2*itertevol+points_t):
	Tcal=np.real(sum((np.conj(np.array(psigrid))[i,:])*(T@ np.array(psigrid)[i,:])*dx))
	Vcal=np.real(sum(psirho[i,:]*psirho[i,:]*Vdot[i,:])*dx)
	exp_energySTA[i]=Tcal+Vcal



plt.title('expectation value of the energy')
plt.ylabel(r'$\langle E(t) \rangle / \hbar \omega$')
plt.xlabel(r'$t/t_{f}$')
plt.plot(exp_energySTA,c='r')
plt.show() 

'''


bounds=1
plt.imshow(psirho, interpolation='nearest',aspect='auto')
plt.title(r'$\rho(x,t)$')
plt.ylabel(r'$t/t_{f}$')
plt.xlabel(r'$x$')
plt.xticks(np.arange(0,(np.shape(V)[1]-2*bounds+1),(np.shape(V)[1]-2*bounds)/(2*int(xf))),np.arange(-int(xf),2*int(xf)+1))
plt.yticks(np.arange(0,np.shape(V)[0],np.shape(V)[0]/4.0),np.linspace(0,1,5))
plt.colorbar()
plt.show()


plt.imshow(phiphase, interpolation='nearest',aspect='auto')
plt.title(r'$\phi(x,t)$')
plt.ylabel(r'$t/t_{f}$')
plt.xlabel(r'$x$')
plt.xticks(np.arange(0,(np.shape(V)[1]-2*bounds+1),(np.shape(V)[1]-2*bounds)/(2*int(xf))),np.arange(-int(xf),2*int(xf)+1))
plt.yticks(np.arange(0,np.shape(V)[0],np.shape(V)[0]/4.0),np.linspace(0,1,5))
plt.colorbar()
plt.show()


np.savetxt('rho'+filename, psirho, delimiter=',')
np.savetxt('phi'+filename, phiphase, delimiter=',')

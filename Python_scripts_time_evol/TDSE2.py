
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
cosa1=501
points_x=cosa1
#tau=int(int(points_x**2 /50)/4.0)
tau=500
points_t=4*tau #para que sea multiplo de 4
DD=-20
steps=int(cosa1*5)
hbar=1.0
m=1.0
W=10.0 #size of the well
d=0.5
cent=-0.0
enerlev=int(sys.argv[1])

# create grid
tf=float(sys.argv[2])
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
'''

def Vpike(A,xvec,d,cent):
    return A*((1-2*(xvec-cent)/d)*window(xvec,d/2,cent+d/4)+(1+2*(xvec-cent)/d)*window(xvec,d/2,cent-d/4))

def Vsq(A,xvec,d,cent):
    return A*window(xvec,d,cent)

def VJar(A1,A2,A3, xvec,cent):
    return A1*(xvec**4 - 0.5*(A2*(xvec)**2*(np.sign( (xvec-cent) )+1)+A3*(xvec)**2*(np.sign( -(xvec-cent) )+1)))

def Vharm(mw2,xvec,cent):
    return 0.5*mw2*(xvec-cent)**2
'''

def Vsqdouble(A1,A3,d1,d3,xvec,cent1,cent2):
    return A1*window(xvec,d1,cent1)+A3*window(xvec,d3,cent2)





#############################################

#############################################




#############################################
#INITIAL CONDITION
#############################################

psigrid=[]
norm=[]
onesies=np.diag(np.ones(points_x))
H=-0.5*(-2.0*np.diag(np.ones(points_x))+np.diag(np.ones(points_x-1),1)+np.diag(np.ones(points_x-1),-1))/(dx**2)
values2=la.eigh(H)

psi1=values2[1].T[enerlev]/np.sqrt(dx)
#psi2=values2[1].T[1]/np.sqrt(dx)
E1=values2[0][enerlev]
#E2=values2[0][2]


#######SETTING UP POTENTIAL

##evolution of the energy levels step 1
TrayV4=np.zeros([points_x,points_t])

for j in range(tau):
	
	V0=Vsqdouble(DD*(j/float(tau-1)),0,2.5,2.5,x,-5+1.25,5-1.25)
	TrayV4[:,j]=V0

for j in range(tau):
	
	V0=Vsqdouble(DD,DD*(j/float(tau-1)),2.5,2.5,x,-5+1.25,5-1.25)
	TrayV4[:,j+tau]=V0


for j in range(tau):
	
	V0=Vsqdouble(DD*(1-j/float(tau-1)),DD,2.5,2.5,x,-5+1.25,5-1.25)
	TrayV4[:,j+2*tau]=V0


for j in range(tau):
	
	V0=Vsqdouble(0,DD*(1-j/float(tau-1)),2.5,2.5,x,-5+1.25,5-1.25)
	TrayV4[:,j+3*tau]=V0

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




###############################################################
V=TrayV4.T
psigrid=[]
norm=[]
psigrid.append(psi1)
onesies=np.diag(np.ones(points_x))
T=-0.5*(-2.0*np.diag(np.ones(points_x))+np.diag(np.ones(points_x-1),1)
        +np.diag(np.ones(points_x-1),-1))/(dx**2)


#plt.plot(psi2**2)
#plt.show()

print("...Initial Norm:",np.sum(np.conj(psigrid[0])*(psigrid[0]))*dx)
print("...Initial eigen-Energy:",E1)
print("...Initial Energy:",np.sum(np.conj(psigrid[0])*(T@psigrid[0]))*dx)



for i in range(points_t):
    psipres=psigrid[i]
    Vc=np.diag(V[i,:])*0
    Hpsi=T+Vc
    matt=onesies-1j*Hpsi*dt
    mattinv=la.inv( onesies +1j*Hpsi*dt)
    psinew=(mattinv@(matt@psipres))
    psigrid.append(psinew)



for i in range(points_t,2*points_t):
    psipres=psigrid[i]
    Vc=np.diag(V[-1,:])*0
    Hpsi=T+Vc
    matt=onesies-1j*Hpsi*dt
    mattinv=la.inv( onesies +1j*Hpsi*dt)
    psinew=(mattinv@(matt@psipres))
    psigrid.append(psinew)



#psirho= np.sqrt(abs(   np.conj(np.array(psigrid))  *np.array(psigrid)   ))
#phiphase=abs( 1j*np.log(   np.array(psigrid) / psirho   )  )
#phiphase=np.arctan2(np.imag(np.array(psigrid)),np.real(np.array(psigrid)))

psirho= np.imag(np.array(psigrid))
#phiphase=abs( 1j*np.log(   np.array(psigrid) / psirho   )  )
phiphase= np.real(np.array(psigrid))

plt.plot(psirho[-1]**2)
#plt.plot(psi2**2)
plt.show()

plt.plot(psirho[0]**2)
plt.plot(psi1**2)
plt.show()


print("...Final Norm:",np.sum(np.conj(psigrid[-1])*(psigrid[-1]))*dx)
print("...Final Energy:",np.sum(np.conj(psigrid[-1])*(T@psigrid[-1]))*dx)

bounds=1
plt.imshow(psirho, interpolation='nearest',aspect='auto')
plt.title(r'$\rho(x,t)$')
plt.ylabel(r'$t/t_{f}$')
plt.xlabel(r'$x$')
plt.xticks(np.arange(0,(np.shape(V)[1]-2*bounds+1),(np.shape(V)[1]-2*bounds)/(2*int(xf))),np.arange(-int(xf),2*int(xf)+1))
plt.yticks(np.arange(0,2*np.shape(V)[0],2*np.shape(V)[0]/4.0),np.linspace(0,1,5))
plt.colorbar()
plt.savefig('rho.png')
plt.show()


plt.imshow(phiphase, interpolation='nearest',aspect='auto')
plt.title(r'$\phi(x,t)$')
plt.ylabel(r'$t/t_{f}$')
plt.xlabel(r'$x$')
plt.xticks(np.arange(0,(np.shape(V)[1]-2*bounds+1),(np.shape(V)[1]-2*bounds)/(2*int(xf))),np.arange(-int(xf),2*int(xf)+1))
plt.yticks(np.arange(0,2*np.shape(V)[0],2*np.shape(V)[0]/4.0),np.linspace(0,1,5))
plt.colorbar()
plt.savefig('phi.png')
plt.show()


#np.savetxt('rho'+filename, psirho, delimiter=',')
#np.savetxt('phi'+filename, phiphase, delimiter=',')

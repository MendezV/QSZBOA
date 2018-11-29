
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
#optimal600

cosa1=501
points_x=cosa1
points_t=4*int(int(points_x**2 /30)/4.0) #para quesea multiplo de 4
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


#########MATRIXMULTIPLICATION########################
def Band_mult(A,x):
	n=np.size(x)
	res=np.zeros(n,dtype=complex)
	res[0]=A[1][0]*x[0]+A[0][0]*x[1]
	for i in range(1,n-1):
		res[i]=A[2][i]*x[i-1]+A[1][i]*x[i]+A[0][i]*x[i+1]

	res[n-1]=A[2][n-1]*x[n-2]+A[1][n-1]*x[n-1]
	return res
#############################################




########TRIDIAGONAL SOLVER#############
def TDMASolve(a, b, c, d):
    nmax = np.size(d)#n in the numbers is rows
    # Modify the first-row coefficients
    c[0] =c[0]/b[0] #Division by zero risk.
    d[0] = d[0]/b[0]
    for i in range(1, nmax):
        ptemp = b[i] - (a[i] * c[i-1])
        c[i] /= ptemp
        d[i] = (d[i] - a[i] * d[i-1])/ptemp
    #Back Substitution
    xres = np.zeros(nmax,dtype=complex)
    xres[-1] = d[-1]
    for i in range(-2,-nmax-1,-1):
	    xres[i] = d[i] - c[i] * xres[i+1]
    return xres

#############################################

#############################################



#############################################
#INITIAL CONDITION
#############################################



psi1=(1/np.sqrt(W))*(np.sin(np.pi*1*(x-W)/(2*W)))
psi2=(1/np.sqrt(W))*(np.sin(np.pi*4*(x-W)/(2*W)))
E1=0.5*(np.pi*1/(2*W))**2
E2=0.5*(np.pi*4/(2*W))**2





######SETTUP FOR TIME EVOLUTION#######
onesies=np.diag(np.ones(points_x))
T=-0.5*(-2.0*np.diag(np.ones(points_x))+np.diag(np.ones(points_x-1),1)
		+np.diag(np.ones(points_x-1),-1))/(dx**2)


A0p=-0.5*np.ones(points_x)/(dx**2)
A0p[0]=0
A1p=np.ones(points_x)/(dx**2)
A2p=-0.5*np.ones(points_x)/(dx**2)
A2p[points_x-1]=0
#######################################


print("...Initial Norm:",np.sum(np.conj(psi1)*(psi1))*dx)
print("...Initial eigen-Energy:",E1)

##plot initial states
'''
plt.plot(x,psi1)
plt.plot(x,psi2)
plt.show()
'''


###############################################################
filename=sys.argv[1]
V=np.loadtxt(filename,delimiter=',')



##plot potential
'''
bounds=1
plt.imshow(V, interpolation='nearest',aspect='auto')
plt.title(r'$\rho(x,t)$')
plt.ylabel(r'$t/t_{f}$')
plt.xlabel(r'$x$')
plt.xticks(np.arange(0,(np.shape(V)[1]-2*bounds+1),(np.shape(V)[1]-2*bounds)/(2*int(xf))),np.arange(-int(xf),2*int(xf)+1))
plt.yticks(np.arange(0,np.shape(V)[0],np.shape(V)[0]/4.0),np.linspace(0,1,5))
plt.colorbar()
plt.show()
'''

psipres=psi1
for i in range(points_t):
	Vc=V[i,:]
	Hpsi=A1p+Vc
	
	A0matt=-1j*A0p*dt
	A1matt=1-1j*Hpsi*dt
	A2matt=-1j*A2p*dt
	mattmult=Band_mult([A2matt,A1matt,A0matt],psipres)
	
	A0mattinv=1j*A0p*dt
	A1mattinv=1+1j*Hpsi*dt
	A2mattinv=1j*A2p*dt
	psinew=TDMASolve(A0mattinv, A1mattinv, A2mattinv, mattmult)
	
	psipres=psinew


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

psirho1=np.array([polarR(psi2[j]) for j in range(points_x)]).T
psirho2=np.array([polarR(psinew[j]) for j in range(points_x)]).T
phiphase1=np.array([polarThe(psi2[j]) for j in range(points_x)]).T
phiphase2=np.array([polarThe(psinew[j]) for j in range(points_x)]).T
print(Ninteg(x,psirho1**2,x0,[xf],dx)[0])
print(Ninteg(x,psirho2**2,x0,[xf],dx)[0])
#psirho= np.sqrt(abs(   np.conj(np.array(psigrid))  *np.array(psigrid)   ))
#phiphase=np.real( 1j*np.log(   np.array(psigrid) / psirho   )  )


##plots

plt.plot(psirho1**2)
plt.title("Initial and final states norm")
plt.plot(psirho2**2)
plt.show()


##plots

plt.plot(np.real(psinew))
plt.plot(np.imag(psinew))
plt.title("final state real and imag")
plt.show()

##plots

plt.plot(phiphase1**2)
plt.title("Initial and final states phase")
plt.plot(phiphase2**2)
plt.show()


##fidelity
print((dx*np.trapz(np.conj(psinew)*psi2))*np.conj(dx*np.trapz(np.conj(psinew)*psi2)),"...fidelity6")


###################################
#PACKAGES
###################################

import numpy as np
from scipy import linalg as la
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import matplotlib
from scipy.linalg import solve_banded
import sys
import seaborn as sns
sns.set_style("white")
###################################

###################################
font = {'size'   : 20 }




###############################################
##CONSTANTS
###############################################

#constants
cosa1=801
points_x=cosa1
tau=int(int(points_x**2)/20.0)
#tau=5000
points_t=4*tau #para que sea multiplo de 4
DD=-0.8

steps=2000
hbar=1.0
m=1.0
W=10.0
d=4
cent=-0.0
# divide by two so a well from -W to W is of input width
W=W/2.0
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

print(dt,dx)
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
frames=tau
##evolution of the energy levels step 1
TrayV4=np.zeros([points_x,points_t])

for j in range(tau):
    
    V0=Vsqdouble(DD*(j/float(frames-1)),0,W-d/2,W-d/2,xvec,-W+(W-d/2)/2,W-(W-d/2)/2)
    TrayV4[:,j]=V0

for j in range(tau):
    
    V0=Vsqdouble(DD,DD*(j/float(frames-1)),W-d/2,W-d/2,xvec,-W+(W-d/2)/2,W-(W-d/2)/2)
    TrayV4[:,j+tau]=V0


for j in range(tau):
    
    V0=Vsqdouble(DD*(1-j/float(frames-1)),DD,W-d/2,W-d/2,xvec,-W+(W-d/2)/2,W-(W-d/2)/2)
    TrayV4[:,j+2*tau]=V0


for j in range(tau):
    
    V0=Vsqdouble(0,DD*(1-j/float(frames-1)),W-d/2,W-d/2,xvec,-W+(W-d/2)/2,W-(W-d/2)/2)
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

A0p=-0.5*np.ones(points_x)/(dx**2)
A0p[0]=0

A1p=np.ones(points_x)/(dx**2)


A2p=-0.5*np.ones(points_x)/(dx**2)
A2p[points_x-1]=0



for i in range(points_t):
	psipres=psigrid[i]
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

	psigrid.append(psinew)



for i in range(points_t,int(3*points_t/2)):
	psipres=psigrid[i]
	Hpsi=A1p
	
	A0matt=-1j*A0p*dt
	A1matt=1-1j*Hpsi*dt
	A2matt=-1j*A2p*dt
	mattmult=Band_mult([A2matt,A1matt,A0matt],psipres)
	
	A0mattinv=1j*A0p*dt
	A1mattinv=1+1j*Hpsi*dt
	A2mattinv=1j*A2p*dt
	psinew=TDMASolve(A0mattinv, A1mattinv, A2mattinv, mattmult)
	
	psigrid.append(psinew)


psirho= np.sqrt(abs(   np.conj(np.array(psigrid))  *np.array(psigrid)   ))
#phiphase=abs( 1j*np.log(   np.array(psigrid) / psirho   )  )
phiphase=np.arctan2(np.imag(np.array(psigrid)),np.real(np.array(psigrid)))

#psiim= np.imag(np.array(psigrid))
#psire= np.real(np.array(psigrid))
enerprime=0
if enerlev>5:
	enerprime=enerlev-5
else:
	enerprime=enerlev+5

ns=np.arange(1,11,1)

'''
plt.plot(psirho[points_t])
plt.plot(phiphase[points_t])
plt.show()


fig = plt.figure()
ax1 = fig.add_subplot()
plt.imshow(psirho, interpolation='nearest',aspect='auto')
#plt.title(r'$\rho(x,t)$')
plt.ylabel(r'$t/4\tau$',size=25)
plt.xlabel(r'$x/L$',size=25)
plt.xticks(np.arange(0,(np.shape(V)[1]+1),(np.shape(V)[1])/(int(xf))),np.arange(-int(xf),2*int(xf)+1,2)/10,fontsize=18)
plt.yticks(np.arange(0,2*np.shape(V)[0],2*np.shape(V)[0]/4.0),np.linspace(0,2,5),fontsize=18)
cb=plt.colorbar()
plt.text(600,3700, r'$|\psi|^{2}$', rotation=-90,fontsize=25)
plt.tight_layout()
plt.show()
'''

weights=[]
for i in range(1,11):
	cn=np.sum(np.array(psigrid[points_t])*np.sqrt(2.0/10.0)*np.sin(i*np.pi*(x+5)/10.0) )*dx
	cn2=abs(np.conj(cn)*cn)
	weights.append(cn2)
print(np.sum(np.array(weights)))
plt.bar(ns,weights,color='b')
plt.xlabel('n',size=25)
plt.ylabel(r"$P_{n \leftarrow 7}(4 \tau )$",size=25)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tight_layout()
plt.savefig("Figure_3b2.png",dpi=500)
plt.show()


print("...Final Norm:",np.sum(np.conj(psigrid[-1])*(psigrid[-1]))*dx)
print("...Final Energy:",np.sum(np.conj(psigrid[-1])*(T@psigrid[-1]))*dx)



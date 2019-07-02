import numpy as np
from scipy import linalg as la
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import sys



#constants
tf=0.24*np.pi
tfFin=tf+tf
xif=1.0
xii=2.0
MaxReach=10
points_x=801
leng=max([xif,xii])
x0,xf=-MaxReach/np.sqrt(leng),MaxReach/np.sqrt(leng)
points_t=int(points_x**2 /50)
print(points_t)
#points_t=500
x=np.linspace(x0,xf,points_x) #debe tener un numero impar de elemntos para que sirva NInteg
t=np.linspace(0,tf,points_t)
itertevol=int(points_t/2)
t2=np.linspace(-tf,tfFin,points_t+2*itertevol)
dt2=tfFin/points_t
dx=(xf-x0)/(points_x)
dt=tf/points_t



########functions

##boundary condditions for the states

def psi(x):
    return (np.pi**(-1.0/4.0))*np.exp(-(x**2)/2.0)

def psiI(x):
    return (xii**(1.0/4.0))*psi(np.sqrt(xii)*x)

def psiF(x):
    return (xif**(1.0/4.0))*psi(np.sqrt(xif)*x)
    #return x*psiI(x)/sqrt(2)

##interpolating function
def eta(t):
    return ((t**3.0)/(tf**3.0))*(1 + 3.0*(1 - (t/tf)) + 6.0*(1 - (t/tf))**2.0)


def Ninteg(x,func,x0,xf,dx):
    spl = UnivariateSpline(x,func, k=3, s=0)
    intfunc=np.zeros(np.size(xf))
    for i in range(np.size(xf)):
        intfunc[i]=spl.integral(x0, xf[i])
 
    return intfunc


def Nderivat(func,x):
    spl = UnivariateSpline(x,func, k=3, s=0)
    derivativas=spl.derivative()
    return derivativas(x)

def Nderivat2(func,x):
    Nderivat(func,x)
    return  Nderivat( Nderivat(func,x),x)


##Time dependent Density
def rho(x,t):
    rho0=(1 - eta(t))*psiI(x) + eta(t)*psiF(x) 
    Z=Ninteg(x,rho0**2,x0,[xf],dx)[0]
    return rho0/np.sqrt(Z)




##energy terms change if the wave functions are not gaussian 
def phi(t):
        return (t/tf)*(1 - (t/tf))*((0.5*xif + 0.5*xii)*t - 0.5*xii*tf)

    
def VI(x):
    return 0.5*(xii**2)*x**2


def VF(x):
    return 0.5*(xif**2)*x**2


def kt(t):
    return (1 - eta(t))*xii + eta(t)*xif

def VT(x,t):
    return 0.5*kt(t)*x**2


###############################################################
filename=sys.argv[1]
V=np.loadtxt(filename,delimiter=',')
psigrid=[]
norm=[]
psigrid.append(psiI(x))
onesies=np.diag(np.ones(points_x))
T=-0.5*(-2.0*np.diag(np.ones(points_x))+np.diag(np.ones(points_x-1),1)
        +np.diag(np.ones(points_x-1),-1))/(dx**2)

print(np.shape(V),points_t,points_x )
#itertevol=1000
for i in range(0,points_t):
    psipres=psigrid[i]
    #T=-0.5*Nderivat2(psipres,x)
    Vc=np.diag(V[i,:])
    Hpsi=T+Vc
    matt=onesies-1j*Hpsi*dt
    mattinv=la.inv( onesies +1j*Hpsi*dt)
    psinew=(mattinv@(matt@psipres))
    psigrid.append(psinew)

for i in range(points_t,points_t+2*itertevol):
    psipres=psigrid[i]
    #T=-0.5*Nderivat2(psipres,x)
    Vc=np.diag(V[-1,:])
    Hpsi=T+Vc
    matt=onesies-1j*Hpsi*dt
    mattinv=la.inv( onesies +1j*Hpsi*dt)
    psinew=(mattinv@(matt@psipres))
    psigrid.append(psinew)


psirho= np.sqrt(abs(   np.conj(np.array(psigrid))  *np.array(psigrid)   ))
phiphase=abs( 1j*np.log(   np.array(psigrid) / psirho   )  )




plt.plot(psirho[-1]**2)
plt.plot(psiF(x)**2)
plt.show()


print(Ninteg(x,psirho[itertevol+points_t]**2,-xf,[xf],dx)[0])

plt.plot(psirho[0]**2)
plt.plot(psiI(x)**2)
plt.show()

####V derivative 
transV=V.T
Temp7=[]
for element7 in transV:
    Temp7.append(Nderivat(element7,t2))
Vdot=np.array(Temp7).T

exp_energySTA=np.zeros(np.size(t)+2*itertevol)

'''
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
plt.xticks(np.arange(0,(np.shape(psirho)[1]-2*bounds+1),(np.shape(psirho)[1]-2*bounds)/(2*int(xf))),np.arange(-int(xf),2*int(xf)+1))
plt.yticks(np.arange(0,np.shape(psirho)[0],np.shape(psirho)[0]/4.0),np.linspace(0,1,5))
plt.colorbar()
plt.show()


plt.imshow(phiphase, interpolation='nearest',aspect='auto')
plt.title(r'$\phi(x,t)$')
plt.ylabel(r'$t/t_{f}$')
plt.xlabel(r'$x$')
plt.xticks(np.arange(0,(np.shape(phiphase)[1]-2*bounds+1),(np.shape(phiphase)[1]-2*bounds)/(2*int(xf))),np.arange(-int(xf),2*int(xf)+1))
plt.yticks(np.arange(0,np.shape(phiphase)[0],np.shape(phiphase)[0]/4.0),np.linspace(0,1,5))
plt.colorbar()
plt.show()


np.savetxt('rho'+filename, psirho, delimiter=',')
np.savetxt('phi'+filename, phiphase, delimiter=',')


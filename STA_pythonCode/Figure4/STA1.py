import numpy as np
from scipy import linalg as la
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt




#constants
tf=0.24*np.pi
tfFin=tf+0.5*tf
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

##boundary conditions for the states
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




######### STA ######################


### state density
rhogrid=np.array([rho(x,tau) for tau in t] )
rhosq=rhogrid**2
transprhogrid=rhosq.T


Temp=[]
for element in transprhogrid:
    Temp.append(Nderivat(element,t))

transtemp= np.array(Temp).T   

Temp2=[]
for element2 in transtemp:
    Temp2.append(Ninteg(x,element2,0,x,dx))


up=np.array(Temp2)
u=up/rhosq
transu=u.T

## first term
Temp3=[]
for element3 in transu:
    Temp3.append(Nderivat(element3,t))

transtemp3=np.array(Temp3).T
Temp4=[]
for element4 in transtemp3:
    Temp4.append(Ninteg(x,element4,0,x,dx))
first=np.array(Temp4)

## second term
Temp5=[]
for element5 in rhogrid:
    Temp5.append(Nderivat2(element5,x))

second=0.5*np.array(Temp5)/rhogrid
## third term
third=-0.5*u**2
##fourth term
fourth=-np.array([Nderivat(phi(t),t) for r in x]).T


###potential
V=second+fourth+third+first


plt.plot(x,V[-1,:])
plt.plot(x,V[0,:])
plt.show()



preV=np.array([VI(x) for i in range(itertevol)])
postV=np.array([VF(x) for i in range(itertevol)])
fullV=np.concatenate((np.concatenate((preV, V), axis=0), postV), axis=0)


bounds=1
plt.imshow(fullV, interpolation='nearest', aspect='auto')
plt.title('V (x,t)')
plt.ylabel(r'$t/t_{f}$')
plt.xlabel(r'$x$') 
plt.xticks(np.arange(0,(np.shape(V)[1]-2*bounds+1),(np.shape(V)[1]-2*bounds)/(2*int(xf))),np.arange(-int(xf),2*int(xf)+1))
plt.yticks(np.arange(0,np.shape(V)[0],np.shape(V)[0]/4.0),np.linspace(0,1,5))
plt.colorbar()

plt.show()




from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X,Y= np.meshgrid(x, t2, sparse=False, indexing='ij')
Z=fullV
ax.plot_surface(X, Y, Z.T)
plt.show()


np.savetxt('VSTA.dat', V, delimiter=',')








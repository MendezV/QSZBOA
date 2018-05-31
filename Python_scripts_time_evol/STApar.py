import numpy as np
from scipy import linalg as la
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
from datetime import datetime
from multiprocessing import Pool, cpu_count



#constants
tf=0.24*np.pi
xif=2.0
xii=1.0
MaxReach=4.5
points_x=2001
x0,xf=-MaxReach/np.sqrt(xif),MaxReach/np.sqrt(xif)
#points_t=points_x**2 /20
points_t=5000
x=np.linspace(x0,xf,points_x) #debe tener un numero impar de elemntos para que sirva NInteg
t=np.linspace(0,tf,points_t)
dx=(xf-x0)/(points_x)
dt=tf/points_t



########functions

##boundary condditions for the states
def psi(x):
    return (np.pi**(-1.0/4.0))*np.exp(-(x**2)/2.0)

#### initial wave function
def psiI(x):
    return (xii**(1.0/4.0))*psi(np.sqrt(xii)*x)

#### final wave function
def psiF(x):
    return (xif**(1.0/4.0))*psi(np.sqrt(xif)*x)
    #return x*psiI(x)/sqrt(2)

##interpolating function
def eta(t):
    return ((t**3.0)/(tf**3.0))*(1 + 3.0*(1 - (t/tf)) + 6.0*(1 - (t/tf))**2.0)

#### definite integral from 0 to each point in the domain
def NintegXarr(func):
    spl = UnivariateSpline(x,func, k=3, s=0)
    intfunc=np.zeros(np.size(xf))
    for i in range(np.size(xf)):
        intfunc[i]=spl.integral(0, x[i])
 
    return intfunc

#### integral over all x in the domain
def NintegXfull(func):
    spl = UnivariateSpline(x,func, k=3, s=0)
    return spl.integral(x0, xf)


#### integral for all time when STA is acting (size f is size t)
def NintegTfull(func):
    spl = UnivariateSpline(t,func, k=3, s=0)
    return spl.integral(0.0, tf)


##### time derivative
def NderivatT(func):
    spl = UnivariateSpline(t,func, k=3, s=0)
    derivativas=spl.derivative()
    return derivativas(t)

#### x derivative for all points in the domain
def NderivatX(func):
    spl = UnivariateSpline(x,func, k=3, s=0)
    derivativas=spl.derivative()
    return derivativas(x)


##### 
def NderivatXX(func):
    return  NderivatX( NderivatX(func))



##Time dependent Density
def rho(x,t):
    rho0=(1 - eta(t))*psiI(x) + eta(t)*psiF(x) 
    Z=NintegXfull(rho0**2)
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

### arr must be a list of arrays and func must only have one argument
def ParalFun(func, arr, num_processes=None):
    ''' Apply a function separately to each column in a dataframe, in parallel.'''
    
    # If num_processes is not specified, default to minimum(#columns, #machine-cores)
    if num_processes==None:
        num_processes = min(np.shape(arr)[1], cpu_count())
    
    # 'with' context manager takes care of pool.close() and pool.join() for us
    with Pool(num_processes) as pool:
       
        # pool.map returns results as a list
        results_list = pool.map(func, arr)
        
        # return list of processed columns, concatenated together as a new dataframe
        return results_list

###############################################################




######### STA ######################


### state density
rhogrid=np.array([rho(x,tau) for tau in t] )
rhosq=rhogrid**2
transprhogrid2=list(rhosq.T)


start = datetime.now()
paraltranstemp = list(np.array(ParalFun(NderivatT, transprhogrid2)).T)
print( datetime.now() - start )

start = datetime.now()
u = np.array(ParalFun( NintegXarr, paraltranstemp))/rhosq
print( datetime.now() - start )



start = datetime.now()
Temp=[]
for element in transprhogrid2:
    Temp.append(NderivatT(element))

transtemp= np.array(Temp).T 
print( datetime.now() - start )

start = datetime.now()
Temp21=[]
for element2 in transtemp:
    Temp21.append(NintegXarr(element2))


up=np.array(Temp21)
print( datetime.now() - start )




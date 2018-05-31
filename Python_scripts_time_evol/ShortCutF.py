
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
cosa1=1001
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
tf=400
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










##########################################
#GENERATE POTENTIAL IN THE USUAL WAY
##########################################
TrayV4=np.zeros([points_x,points_t])



frames=int(points_t/4)
for j in range(frames):
    
    V=Vsqdouble(DD*(j/float(frames-1)),0,2.5,2.5,x,-5+1.25,5-1.25)
    TrayV4[:,j]=V
    
for j in range(frames):
    
    V=Vsqdouble(DD,DD*(j/float(frames-1)),2.5,2.5,x,-5+1.25,5-1.25)
    TrayV4[:,j+frames]=V
    


for j in range(frames):
 
    V=Vsqdouble(DD*(1-j/float(frames-1)),DD,2.5,2.5,x,-5+1.25,5-1.25)
    TrayV4[:,j+2*frames]=V
   

for j in range(frames):

    V=Vsqdouble(0,DD*(1-j/float(frames-1)),2.5,2.5,x,-5+1.25,5-1.25)
    TrayV4[:,j+3*frames]=V



#############################################

#############################################








################################################################
#SHORTCUT TO ADIABATICITY 
################################################################




#########INTEGRATION########################
from scipy.interpolate import UnivariateSpline

def Ninteg(x,func,x0,xf,dx):
    spl = UnivariateSpline(x,func, k=3, s=0)
    derivativas=spl.derivative()
    intfunc=np.zeros(np.size(xf))
    for i in range(np.size(xf)):
        intfunc[i]=spl.integral(0, xf[i])
 
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






###########FUNCTION###########
def Shortcut(x,psiI,psiF, E_I, E_F, eta):
    PsiI= UnivariateSpline(x,psiI, k=5, s=0)
    PsiF= UnivariateSpline(x,psiF, k=5, s=0)
    
    def rho(x,t):
        rho0=(1 - eta(t))*PsiI(x) + eta(t)*PsiF(x) 
        Z=np.trapz(rho0**2)*dx
        return rho0/np.sqrt(Z)
##energy terms change if the wave functions are not gaussian 
    def phi(t, E_I, E_F):
        return (t/tf)*(1 - (t/tf))*(( E_F + E_I)*t - E_I*tf)
    rhogrid=np.array([rho(x,tau) for tau in t] )

    Temp=[]
    for element in rhogrid.T:
        Temp.append(Nderivat(element**2,t))
    

    Temp21=[]
    for element2 in np.array(Temp).T:
        Temp21.append(Ninteg(x,element2,0,x,dx))
        
    up=np.array(Temp21)
    rhosq=rhogrid**2
    #rhosq[where(rhosq<10**-9)]=0.0
    #up[where(abs(up)<10**-9)]=0.0
    u=up/rhosq
    #u[isnan(u)]=0.0
    #u[isinf(u)]=0.0
    plt.title("hydrodinamic velocity vs x")
    plt.ylabel("u(x,0)")
    plt.xlabel("x")
    plt.plot(x,u[0,:])
    plt.show()

    plt.imshow(u, interpolation='nearest', aspect='auto')
    plt.colorbar()
    plt.show()
    
     ## first term

    Temp3=[]
    for element3 in u.T:
        Temp3.append(Nderivat(element3,t))
    Temp4=[]
    for element4 in np.array(Temp3).T:
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
    #print(np.shape(first),np.shape(second),np.shape(third), np.shape(t), np.shape(tf),np.shape(E_I),np.shape(E_F))
    fourth=-np.array([Nderivat(phi(t, E_I, E_F),t) for r in x]).T
    print(np.shape(first),np.shape(second),np.shape(third),np.shape(fourth))
    Vres=second+fourth
    #Vres=second+fourth+first+third
    
    c=14*10**5
    Vres[np.where( (Vres)>c )]=c
    Vres[np.where( (Vres)<-c )]=-c
    slices=9
    for i in range(0,slices):
    #title("Interpolating function for the state density")
        plt.ylabel(r"$V$")
        plt.xlabel("x")
        plt.ylim([-100,100])
        #plt.xlim([-0.7,0.7])
        plt.plot(x,Vres[int(i*(np.size(t)-1)/(slices-1)),:],label="interpolating",c='r')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.savefig(str(i)+'.png')
        plt.show()



    return Vres
##############################




####################################################################

####################################################################



psigrid=[]
norm=[]
onesies=np.diag(np.ones(points_x))
T=-0.5*(-2.0*np.diag(np.ones(points_x))+np.diag(np.ones(points_x-1),1)
        +np.diag(np.ones(points_x-1),-1))/(dx**2)+np.diag(Vsqdouble(DD,DD,2.5,2.5,x,-5+1.25,5-1.25))

values2=la.eigh(T)

'''
for i in range(10):
	psi2=10*values2[1].T[i]/np.sqrt(dx)
	E2=values2[0][i]
	plt.plot(x,psi2+E2*10)
	plt.plot(x,10*E2*psi2/psi2)
plt.show()
'''

psi1=values2[1].T[0]/np.sqrt(dx)
psi2=values2[1].T[2]/np.sqrt(dx)
E1=values2[0][0]
E2=values2[0][2]
plt.plot(x,psi1)
plt.plot(x,psi2)
plt.show()

def ShortcutR(x,psiI,psiF, E_I, E_F, eta):
    PsiI= UnivariateSpline(x,psiI, k=5, s=0)
    PsiF= UnivariateSpline(x,psiF, k=5, s=0)
    
    def rho(x,t):
        rho0=(1 - eta(t))*PsiI(x) + eta(t)*PsiF(x) 
        Z=np.trapz(rho0**2)*dx
        return rho0/np.sqrt(Z)
##energy terms change if the wave functions are not gaussian 
    def phi(t, E_I, E_F):
        return (t/tf)*(1 - (t/tf))*(( E_F + E_I)*t - E_I*tf)
    # Two subplots, unpack the axes array immediately
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    ax1.plot(x,np.sqrt(PsiI(x)**2))
    ax1.set_title('Initial Density')
    ax1.set_ylabel(r'$\rho (x,0)$')
    ax1.set_xlabel(r'$x$')
    ax2.plot(x,np.sqrt(PsiF(x)**2),c='g')
    ax2.set_title('Target Density')
    ax2.set_ylabel(r'$\rho (x,t_{f})$')
    ax2.set_xlabel(r'$x$')

    plt.show()
    rhogrid=np.array([rho(x,tau) for tau in t] )

    

    return  rhogrid

R=ShortcutR(x,psi1,psi2, E1, E2, eta)

bounds=1
plt.imshow(R[:,bounds:-bounds], interpolation='nearest', aspect='auto')
plt.title('V (x,t)')
plt.ylabel(r'$t/t_{f}$')
plt.xlabel(r'$x$') 
plt.xticks(np.arange(0,(np.shape(R)[1]-2*bounds+1),(np.shape(R)[1]-2*bounds)/4),np.arange(-2,4))
plt.yticks(np.arange(0,np.shape(R)[0],np.shape(R)[0]/4.0),np.linspace(0,1,5))
plt.colorbar()

plt.show() 

VRR=Shortcut(x,psi1,psi2, E1, E2, eta)

plt.imshow(VRR[:,bounds:-bounds], interpolation='nearest', aspect='auto')
plt.title('V (x,t)')
plt.ylabel(r'$t/t_{f}$')
plt.xlabel(r'$x$') 
plt.xticks(np.arange(0,(np.shape(R)[1]-2*bounds+1),(np.shape(R)[1]-2*bounds)/4),np.arange(-2,4))
plt.yticks(np.arange(0,np.shape(R)[0],np.shape(R)[0]/4.0),np.linspace(0,1,5))
plt.colorbar()

plt.show() 


np.savetxt('VST.dat', VRR, delimiter=',')




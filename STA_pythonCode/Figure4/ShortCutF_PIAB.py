
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
cosa1=601
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
    spl = UnivariateSpline(x,func, k=5, s=0)
    intfunc=np.zeros(np.size(xf))
    for i in range(np.size(xf)):
        intfunc[i]=spl.integral(x0, xf[i])
 
    return intfunc

def Ninteg2(x,func,x0,xf,dx):
    spl = UnivariateSpline(x,func, k=5, s=0)
    intfunc=spl.integral(x0, xf)
    return intfunc


def eta(t):
    return ((t**3.0)/(tf**3.0))*(1 + 3.0*(1 - (t/tf)) + 6.0*(1 - (t/tf))**2.0)
#############################################










#########DIFFERENTIATION########################

def Nderivat(func,x):
    spl = UnivariateSpline(x,func, k=5, s=0)
    derivativas=spl.derivative()
    return derivativas(x)

def Nderivat2(func,x):
    Nderivat(func,x)
    return  Nderivat( Nderivat(func,x),x)

#############################################








#####################################

###########FUNCTION################



def Shortcut(x,psiI,psiF, E_I, E_F, eta):
    PsiI= UnivariateSpline(x,psiI, k=5, s=0)
    PsiF= UnivariateSpline(x,psiF, k=5, s=0)
    
    def rho(x,t):
        rho0=(1 - eta(t))*PsiI(x) + eta(t)*PsiF(x) 
        Z=Ninteg(x,rho0**2,x0,[xf],dx)[0]
        return rho0/np.sqrt(Z)


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



    Temp=[]
    for element in rhogrid:
        Temp.append(Ninteg(x,element**2,0,x,dx))
    

    Temp21=[]
    for element2 in np.array(Temp).T:
        Temp21.append(Nderivat(element2,t))


    up=np.array(Temp21).T


    rhosq=rhogrid**2
    u=up/rhosq
    u[np.where( np.isnan(u) )]=0
    '''
    plt.title("hydrodinamic velocity vs x")
    plt.ylabel("u(x,0)")
    plt.xlabel("x")
    plt.plot(x,u[0,:])
    plt.show()

    plt.imshow(u, interpolation='nearest', aspect='auto')
    plt.colorbar()
    plt.show()
    '''
    

    ## first term
    Temp3=[]
    for element3 in u.T:
        Temp3.append(Nderivat(element3,t))
    Temp4=[]
    for element4 in np.array(Temp3).T:
        Temp4.append(Ninteg(x,element4,0,x,dx))
    first=np.array(Temp4)
    first[np.where( np.isnan(first) )]=0

    ## second term
    Temp5=[]
    for element5 in rhogrid:
        Temp5.append(Nderivat2(element5,x))
    second=0.5*np.array(Temp5)/rhogrid
    second[np.where( np.isnan(second) )]=0



    ## third term
    third=-0.5*u**2
    third[np.where( np.isnan(third) )]=0

    ##fourth term
    fourth=-np.array([Nderivat(phi(t, E_I, E_F),t) for r in x]).T
    fourth[np.where( np.isnan(fourth) )]=0
    print(np.shape(first),np.shape(second),np.shape(third),np.shape(fourth))
    


    ##final result
    #Vres=second+fourth
    Vres=second+fourth+third+first
    Vres[np.where( np.isnan(Vres) )]=0

    ##cutoff
    c=10
    Vres[np.where( np.isnan(Vres) )]=0
    Vres[np.where( (Vres)>c )]=c
    Vres[np.where( (Vres)<-c )]=-c
    '''
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
    '''
    for element5 in u:
        Temp5.append(Ninteg(x,element5,0,x,dx))
    phires=np.array(Temp4)

    
    return Vres, phires
###############################################################

####################################################################







####################################################################
################# boundary conditions #####################


psigrid=[]
norm=[]
onesies=np.diag(np.ones(points_x))
T=-0.5*(-2.0*np.diag(np.ones(points_x))+np.diag(np.ones(points_x-1),1)+np.diag(np.ones(points_x-1),-1))/(dx**2)+np.diag(Vsqdouble(DD,DD,2.5,2.5,x,-5+1.25,5-1.25))   #first hamiltonian

#T=-0.5*(-2.0*np.diag(np.ones(points_x))+np.diag(np.ones(points_x-1),1)+np.diag(np.ones(points_x-1),-1))/(dx**2)+np.diag( Vharm(1,x,0))   #first hamiltonian




T2=-0.5*(-2.0*np.diag(np.ones(points_x))+np.diag(np.ones(points_x-1),1)+np.diag(np.ones(points_x-1),-1))/(dx**2)+np.diag(Vsqdouble(0.05*DD,0.05*DD,2.5,2.5,x,-5+1.25,5-1.25))   #second hamiltonian


values2=la.eigh(T)
values=la.eigh(T2)

'''
#correct for general potential
psig=values2[1].T[1] #ground state for first potental
psie=values2[1].T[6] #excited state for first potental
Eg=values2[0][1]   #ground state energy for first potental
Ee=values2[0][6]   #excited state energy for first potental


psigprime=values[1].T[1] #ground state for first potental
psieprime=values[1].T[6] #excited state for first potental
Egprime=values[0][1] #ground state energy for first potental
Eeprime=values[0][6] #excited state energy for first potental
'''

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
'''


'''
##for 1-2
#correction for exponentially small splittings
psig=(1/np.sqrt(2))*(values2[1].T[0]+values2[1].T[1]) #ground state for first potental
psie=(1/np.sqrt(2))*(values2[1].T[0]-values2[1].T[1]) #excited state for first potental
Eg=values2[0][0]   #ground state energy for first potental
Ee=values2[0][1]   #excited state energy for first potental


psigprime=(1/np.sqrt(2))*(values[1].T[0]+values[1].T[1]) #ground state for first potental
psieprime=(1/np.sqrt(2))*(values[1].T[0]-values[1].T[1]) #excited state for first potental
Egprime=values[0][0] #ground state energy for first potental
Eeprime=values[0][1] #excited state energy for first potental
'''


'''
#interpolated states "prime"
psi1=psig/np.sqrt( Ninteg2(x,psig**2,x0,xf,dx) )
psi2=psie/np.sqrt( Ninteg2(x,psie**2,x0,xf,dx) )
E1=Eg
E2=Ee

#print(Ninteg(x,psi1**2,x0,[xf],dx)[0])
#print(Ninteg(x,psi2**2,x0,[xf],dx)[0])
'''
##for 2-7

psi1=(1/np.sqrt(W))*(np.sin(np.pi*1*(x-W)/(2*W)))
psi2=(1/np.sqrt(W))*(np.sin(np.pi*5*(x-W)/(2*W)))
E1=0.5*(np.pi*1/(2*W))**2
E2=0.5*(np.pi*5/(2*W))**2
'''
plt.plot(x,psi1)
plt.plot(x,psi2)
plt.show()
'''

#################################################################

##################################################################





def ShortcutR(x,psiI,psiF, E_I, E_F, eta):
    PsiI= UnivariateSpline(x,psiI, k=5, s=0)
    PsiF= UnivariateSpline(x,psiF, k=5, s=0)
    
    def rho(x,t):
        rho0=(1 - eta(t))*PsiI(x) + eta(t)*PsiF(x) 
        Z=Ninteg(x,rho0**2,x0,[xf],dx)[0]
        return rho0/np.sqrt(Z)
##energy terms change if the wave functions are not gaussian 
    def phi(t, E_I, E_F):
        return (t/tf)*(1 - (t/tf))*(( E_F + E_I)*t - E_I*tf)
    # Two subplots, unpack the axes array immediately
    '''
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    ax1.plot(x,PsiI(x))
    ax1.set_title('Initial Density',size=30)
    ax1.set_ylabel(r'$\rho (x,0)$',size=30)
    ax1.set_xlabel(r'$x (a.u)$',size=30)
    ax2.plot(x,PsiF(x),c='g')
    ax2.set_title('Target Density',size=30)
    ax2.set_ylabel(r'$\rho (x,t_{f})$',size=30)
    ax2.set_xlabel(r'$x$',size=30)
    plt.tight_layout()

    plt.show()
	'''
    rhogrid=np.array([rho(x,tau) for tau in t] )
    '''
    f, ax1 = plt.subplots()
    ax1.plot(x/10,PsiI(x),c='b',label='t=0')
    ax1.plot(x/10,PsiF(x),c='g',label='$t=t_f$')
    ax1.set_ylabel(r'$\rho (x,t)$',size=30)
    ax1.set_xlabel(r'$x/L$',size=30)
    plt.legend(fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.tight_layout()
    plt.show()
	'''
    return  rhogrid

R=ShortcutR(x,psi1,psi2, E1, E2, eta)
'''
bounds=1
plt.imshow(R[:,bounds:-bounds], interpolation='nearest', aspect='auto')
plt.title(r'$\rho (x,t)$',size=30)
plt.ylabel(r'$t/t_{f}$',size=30)
plt.xlabel(r'$x/L$',size=30) 
plt.xticks(np.arange(0,(np.shape(R)[1]-2*bounds+1),(np.shape(R)[1]-2*bounds)/4),(0.5/2)*np.arange(-2,4),fontsize=25)
plt.yticks(np.arange(0,np.shape(R)[0],np.shape(R)[0]/4.0),np.linspace(0,1,5),fontsize=25)
plt.colorbar()
plt.tight_layout()
plt.show() 

'''

VRR, PhiR=Shortcut(x,psi1,psi2, E1, E2, eta)
'''
bounds=1
plt.imshow(VRR[:,bounds:-bounds], interpolation='nearest', aspect='auto')
plt.title('V(x,t)')
plt.ylabel(r'$t/t_{f}$')
plt.xlabel(r'$x$') 
plt.xticks(np.arange(0,(np.shape(R)[1]-2*bounds+1),(np.shape(R)[1]-2*bounds)/4),np.arange(-2,4))
plt.yticks(np.arange(0,np.shape(R)[0],np.shape(R)[0]/4.0),np.linspace(0,1,5))
plt.colorbar()

plt.show()
'''
'''
plt.imshow(PhiR[:,bounds:-bounds], interpolation='nearest', aspect='auto')
plt.title('phi (x,t)')
plt.ylabel(r'$t/t_{f}$')
plt.xlabel(r'$x$')
plt.xticks(np.arange(0,(np.shape(R)[1]-2*bounds+1),(np.shape(R)[1]-2*bounds)/4),np.arange(-2,4))
plt.yticks(np.arange(0,np.shape(R)[0],np.shape(R)[0]/4.0),np.linspace(0,1,5))
plt.colorbar()

plt.show() 

plt.plot(PhiR[-1,bounds:-bounds])
plt.plot(PhiR[0,bounds:-bounds])
plt.ylim([-1,7])
plt.show()
'''

np.savetxt('VST27_100.dat', VRR, delimiter=',')



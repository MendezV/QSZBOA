

###################################
#PACKAGES
###################################

import numpy as np
from scipy import linalg as la
import sys

###################################

###################################

###############################################
##CONSTANTS
###############################################

#constants
cosa1=1001
points_x=cosa1
tau=25*int(int(points_x**2)/2.0)
points_t=4*tau #para que sea multiplo de 4
DD=-20
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
#t=np.linspace(0,tf,points_t)

print(tf,dt,dx)


#############################################

#############################################





############################################
#POTENTIALS
############################################


def window(xvec,d,cent):
    return 0.5*((np.sign(xvec+d/2.0-cent))-(np.sign(xvec-d/2.0-cent)))

def Vsqdouble(A1,A3,d1,d3,xvec,cent1,cent2):
    return A1*window(xvec,d1,cent1)+A3*window(xvec,d3,cent2)


#############################################

#############################################


#############################################
#INITIAL CONDITION
#############################################

#psi1=values2[1].T[enerlev]/np.sqrt(dx)
psi1=np.sqrt(2.0/(2*W))*np.sin(np.pi*(enerlev+1)*(x+W)/(2*W))
E1=0.5*np.pi*np.pi*(enerlev+1)*(enerlev+1)/((2*W)**2)


#############################################

#############################################



#############################################
##FUNCTIONS
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
###MAIN
###############################################################




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
print("...Initial Energy:",np.sum(np.conj(psi1)*(np.matmul(T.T, psi1)))*dx)
print(".....Initial state 1+"+sys.argv[1])

psipres=psi1
for i in range(tau):
        Vc=Vsqdouble(DD*(i/float(tau-1)),0,2.5,2.5,x,-5+1.25,5-1.25)
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



for i in range(tau):
        Vc=Vsqdouble(DD,DD*(i/float(tau-1)),2.5,2.5,x,-5+1.25,5-1.25)
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

for i in range(tau):
        Vc=Vsqdouble(DD*(1-i/float(tau-1)),DD,2.5,2.5,x,-5+1.25,5-1.25)
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

for i in range(tau):
        Vc=Vsqdouble(0,DD*(1-i/float(tau-1)),2.5,2.5,x,-5+1.25,5-1.25)
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

E0=np.abs(np.sum(np.conj(np.array(psipres))*(np.matmul(T,np.array(psipres)) ) )*dx)
print("...Final Norm:",np.sum(np.conj(psinew)*(psinew))*dx)
print("...Final Energy:",E0)
###################################################################################

###################################################################################




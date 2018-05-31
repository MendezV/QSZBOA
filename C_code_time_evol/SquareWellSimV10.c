// motion of a wavepacket incident on a potential
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>


const double PI = 3.14159265358979323846264338327950;

void parameters(double *x0, double *k0, double *width, double *V0, double *a,
				double *xmin, double *xmax, double *n, double *dx, double *dx2,
				double *dt, double *k1, double *nt);
void initial_packet(double complex Psi[], double x0, double k0,
					double width, double xmin, double n, double dx, double dt);
void initial_packet2(double complex Psi[], double x0, double k1,
					 double a, double xmin, double n, double dx, double dt);
void initial_packet3(double complex Psi[], double x0, double k0,
					 double width, double xmin, double n, double dx, double dt);
void evolve(double complex Psi[], double complex Psiold[], double *t,
			double V0, double a, double dx, double dx2, double dt,
			double xmin, double n, double nt, int j);

double V(double x, double V0, double a, double tau, double Tt);
double Len2( double a, double tau, double Tt);
double Len( double a, double tau, double Tt);

void evolve2(double complex Psi[], double complex Psiold[], double *t,
			 double V0, double a, double dx, double dx2, double dt,
			 double xmin, double n, double Tt, int j);
double Omegsq( double a, double tau, double Tt);

double Lenprime( double a, double tau, double Tt);
double complex determinant(double complex *Uhalf, double k);

int main(void)
{
	double magnitud=0.0;
	int i,j,k;
	double Tt=0.0001;
	FILE *in;
	int var;
	int test;
	char filename[100]="packet201.dat";
	double complex *Psi;
	double complex *Psiold;
	
	
	
	double x0, k0, width, V0, a, xmin, xmax, n, dx, dx2,nt, dt,k1;
	parameters(&x0, &k0, &width, &V0, &a, &xmin, &xmax, &n, &dx, &dx2, &dt,&k1,&nt);
	
	Psiold=malloc(n*sizeof(double complex));
	Psi=malloc(n*sizeof(double complex));
	
	/*initial_packet(Re, Im, x0, k0, width, xmin, n, dx, dt);*/
	initial_packet2(Psi, a/2.0, k1, a, xmin, n, dx, dt);
	/*initial_packet3(Re, Im, x0, k0, width, xmin, n, dx, dt);*/
	
	
	in = fopen(filename,"w");
	if(!in){
		printf("problems opening the file %s\n", filename);
		exit(1);
	}
	
	
	double t = 0.0;
	for(j=0;j<nt;j++){
		
		for (int k=0; k < n; k++) {
			
			Psiold[k]=Psi[k];
		}
		
		evolve(Psi , Psiold, &t, V0, a, dx, dx2, dt, xmin, n, Tt,j);
		if (j%100==0) {
			for(i=0;i<n;i++){
				magnitud=conj(Psi[i])*Psi[i];
				fprintf(in,"%f ",magnitud);
				
			}
			fprintf(in,"\n");
		}
		
	}
	
	for(j=0;j<0;j++){
		evolve2(Psi, Psiold, &t, V0, a, dx, dx2, dt, xmin, n, Tt,j);
		if (j%100==0) {
			for(i=0;i<n;i++){
				magnitud=conj(Psi[i])*Psi[i];
				fprintf(in,"%f ",magnitud);
				
			}
			fprintf(in,"\n ");
		}
		
	}
	
	
	
	
	
	fclose(in);
	
	
	
	return EXIT_SUCCESS;
}

void parameters(double *x0, double *k0, double *width, double *V0, double *a,
				double *xmin, double *xmax, double *n, double *dx, double *dx2,
				double *dt, double *k1, double *nt)
{
	*x0 = 0;
	*width = 1;
	*k0 = 4;
	*xmax = 2;
	*xmin = 0;
	*V0 = 0;
	*a = 1.0; /* initial size of the well*/
	*k1= PI/(*a);
	*dx = 0.01;
	*dx2 = (*dx)*(*dx);
	*n = (*xmax - *xmin)/(*dx);
	*nt = 100000;
	*dt=0.0001/(*nt);
}

void initial_packet(double complex Psi[], double x0, double k0,
					double width, double xmin, double n, double dx, double dt)
// initial Gaussian wavepacket
{
	
	double delta2 = 0.01;
	double A = pow(2*PI*delta2, -0.25);
	/*double b = 0.5*k0*dt;*/
	for (int i=1; i <= n; i++) {
		double x = xmin + (i-1)*dx;
		double arg = 0.25*pow(x - x0, 2)/delta2;
		double e = exp(-arg);
		Psi[i] = A*cos(1000*k0*(x - x0))*e+I*A*sin(1000*k0*(x - x0+2))*e;
		
	}
}

void initial_packet2(double complex Psi[], double x0, double k1,
					 double a, double xmin, double n, double dx, double dt)
// initial sine wavepacket
{
	
	double A = pow(2.0/a, 0.5);
	/*double b = 0.5*k0*dt;*/
	for (int i=1; i <= n; i++) {
		double x = xmin + (i-1)*dx;
		Psi[i] = 0.0;

		
	}
	for (int i=(int)(-xmin)/dx ; i <= (int)(a-xmin)/dx; i++) {
		double x = (i-1)*dx;
		/*Re[i] = A*sin(k1*(x - a))*cos((100)*(x - x0));
		 Im[i] = A*sin(k1*(x - a))*sin((100)*(x - x0));*/
		Psi[i] = A*sin(k1*(x));
		
		
	}
}

void initial_packet3(double complex Psi[], double x0, double k0,
					 double width, double xmin, double n, double dx, double dt)
// initial Gaussian wavepacket
{
	
	double delta2 = width*width;
	double A = pow(2*PI*delta2, -0.25);
	/*double b = 0.5*k0*dt;*/
	for (int i=1; i <= n; i++) {
		double x = xmin + (i-1)*dx;
		Psi[i] = 0.0;
		
	}
	for (int i=(int)(-xmin)/dx ; i <= (int)(1.0-xmin)/dx; i++) {
		double x = (i-1)*dx;
		double arg = 0.25*pow(x , 2)/delta2;
		double e = exp(-arg);
		Psi[i] = A*cos(k0*(x - x0))*e+I*A*sin(k0*(x - x0))*e;;

	}
}

/* np sirve con funciones de onda discontinuas, ni deberia*/
void evolve(double complex Psi[], double complex Psiold[], double *t,
			double V0, double a, double dx, double dx2, double dt,
			double xmin, double n, double Tt, int j)
{
	
	
	double tau=*t;
	double *Uhalf;
	double complex *coefs;
	int n_x=(int)(Len(a, Tt, Tt)/dx);
	int n_y=(int)(Len(a, Tt, Tt)/dx);
	int pos,pos2,pos3,pos4,pos5,pos6;
	int i,k;
	int init=(int)((-xmin)/dx );
	int fini=(int)((Len(a, Tt, Tt)-xmin)/dx);
	double complex sum;
	double x;
	
	
	printf("%f %f %f \n", tau/Tt , Tt*Tt*Omegsq( a, tau, Tt), Len(a, tau, Tt));
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/* setting up the time evolution matrix that will be invrted for the  crank-nicolson method*/

	Uhalf=malloc(n*n*sizeof(double complex));
	coefs=malloc(6*n*sizeof(double complex));
	
	/*first everything is zero*/
	for (int i=0; i < n; i++) {
		for (int k=0; k < n; k++) {
			
			pos = k + (n * i) ;  /*position in the array*/
			
			Uhalf[pos]=0.0;
		}
		
		
	}
	
	/*submatrix of elements that correspond to the subspace that will be evolved half a step*/
	
	pos = init + (n_y * init);
	Uhalf[pos]=1.0+I*dt*(I*1.0/dx2  -V(x,V0,a,tau,Tt) )/2.0;
	pos = init + 1 + (n_y * init);
	Uhalf[pos]=-I*dt*(I*0.5/dx2)/2.0;
	pos = init + (n_y * (init+1));
	Uhalf[pos]=-I*dt*(I*0.5/dx2)/2.0;

	
	pos = fini-1 + (n_y * (fini-1));
	Uhalf[pos]=1.0+I*dt*(I*1.0/dx2  -V(x,V0,a,tau,Tt))/2.0;
	pos = fini-2 + (n_y * (fini-1));
	Uhalf[pos]=-I*(I*0.5*dt/dx2)/2.0;
	pos = fini-1 + (n_y * (fini-2));
	Uhalf[pos]=-I*(I*0.5*dt/dx2)/2.0;
	
	for (int i=init+1; i < fini-1; i++) {
		x = (i-1)*dx;
		
		for (int k=init+1; k < fini-1; k++) {
			
			pos = k + (n_y * i);
			
			if(i==k){
				Uhalf[pos]=1.0+I*dt*(I*1.0/dx2)/2.0;
			}
			if(i==(k+1)){
				Uhalf[pos]=-I*dt*(I*0.5/dx2)/2.0;
			}
			if(i==(k-1)){
				Uhalf[pos]=-I*dt*(I*0.5/dx2)/2.0;
			}
			
		}
		
		
	}
	
	for (int i=init; i < fini; i++) {
		sum=0.0;
		for (int k=init; k < fini; k++) {
			pos = i + (n * k) ;
			sum+=Uhalf[pos]*Psiold[k];
			
		}
		pos = 3 + (6 * i) ;
		coefs[pos]=sum;
		/*printf("%f \n",creal(Psi[i]));*/
	}
	
	
	/*setting up coefs of the linear system to perform tridiagonal inversion with the Thomas algorithm*/
	x = (init-1)*dx;
	pos = 0 + (6 * init);
	coefs[pos]=0;
	pos = 1 + (6 * init);
	coefs[pos]=1.0-I*dt*(I*1.0/dx2  -V(x,V0,a,tau,Tt))/2.0;
	pos = 2 + (6 * init);
	coefs[pos]=I*dt*(I*0.5/dx2)/2.0;
	
	
	for (int i=init+1; i < fini-1; i++) {
		x = (i-1)*dx;
		
		pos = 0 + (6 * i);
		coefs[pos]=I*dt*(I*0.5/dx2)/2.0;
		pos = 1 + (6 * i);
		coefs[pos]=1.0-I*dt*(I*1.0/dx2 -V(x,V0,a,tau,Tt))/2.0;
		pos = 2 + (6 * i);
		coefs[pos]= I*dt*(I*0.5/dx2)/2.0;
		
		
	}
	x = ((fini-1)-1)*dx;
	pos = 0 + (6 * (fini-1));
	coefs[pos]= I*dt*(I*0.5/dx2)/2.0;
	pos = 1 + (6 * (fini-1));
	coefs[pos]=1.0-I*dt*(I*1.0/dx2 -V(x,V0,a,tau,Tt))/2.0;
	pos = 2 + (6 * (fini-1));
	coefs[pos]=0;
	
	
	
	/*Explicit thomas algorithm for finding the coeficients that solve the linear system*/
	
	pos = 1 + (6 * init);
	pos2 = 2 + (6 * init);
	pos3= 4 + (6 * init);
	coefs[pos3]=coefs[pos2]/coefs[pos];
	
	
	pos = 1 + (6 * init);
	pos2 = 3 + (6 * init);
	pos3 = 5 + (6 * init);
	coefs[pos3]= coefs[pos2]/coefs[pos];
	
	
	for (int i=init+1; i < fini-1; i++) {
		
		pos = 0 + (6 * i);
		pos2 = 1 + (6 * i);
		pos3= 2 + (6 * i);
		pos4 = 4 + (6 * (i-1));
		pos5 = 4 + (6 * i);
		coefs[pos5]= coefs[pos3]/(coefs[pos2]-coefs[pos]*coefs[pos4]);
		
		
		pos = 0 + (6 * i);
		pos2 = 1 + (6 * i);
		pos3= 3 + (6 * i);
		pos4 = 4 + (6 * (i-1));
		pos5 = 5 + (6 * (i-1));
		pos6= 5 + (6 * i);
		coefs[pos6]=(coefs[pos3]-coefs[pos5]*coefs[pos])/(coefs[pos2]-coefs[pos]*coefs[pos4]);
		/*printf("%f %f \n",creal((coefs[pos3]-coefs[pos6]*coefs[pos])),creal((coefs[pos2]-coefs[pos]*coefs[pos4])));*/
		
	}
	
	
	pos3= 4+ (6 * (fini-1));
	coefs[pos3]=0;
	
	
	pos = 0 + (6 * (fini-1));
	pos2 = 1 + (6 * (fini-1));
	pos3= 3 + (6 * (fini-1));
	pos4 = 4 + (6 * (fini-2));
	pos5 = 5 + (6 * (fini-2));
	pos6= 5 + (6 * (fini-1));
	coefs[pos6]=(coefs[pos3]-coefs[pos5]*coefs[pos])/(coefs[pos2]-coefs[pos]*coefs[pos4]);
	
	/*for(int i=init; i < fini; i++) {
		printf("%f %f %f %f %f %f\n",creal(coefs[0 + (6 * i)]),creal(coefs[1 + (6 * i)]),creal(coefs[2 + (6 * i)]),creal(coefs[3 + (6 * i)]),creal(coefs[4 + (6 * i)]),creal(coefs[5 + (6 * i)]));
	}*/
	
	/* solving the linear system  with the thomas coefficients*/
	
	for (int i=fini-1; i >= init; i--) {
		pos= 4 + (6 * i);
		pos2= 5 + (6 * i);
	    Psi[i]=coefs[pos2]-coefs[pos]*Psi[i+1];
		/*printf("%f \n",creal(Psi[i]));*/
		
	}

	*t += dt;
	
	
}




void evolve2(double complex Psi[], double complex Psiold[], double *t,
			 double V0, double a, double dx, double dx2, double dt,
			 double xmin, double n, double Tt, int j)
{
	
}


double V(double x, double V0, double a, double tau, double Tt)
{
	
	
	return x*x*Omegsq( a, tau, Tt)/2.0;
	/*return 0;*/
	
}

double Omegsq( double a, double tau, double Tt)
{
	double gamma=1;
	double s=tau/Tt;
	double gammadprime=(gamma-1)*(60.0*s*(1.0+s*(2.0*s-3.0)))/(Tt*Tt);
	
	return -gammadprime/Len(1.0, tau, Tt);
	
	
}

double Len( double a, double tau, double Tt){
	
	double gamma=1;
	double s=tau/Tt;
	return a*(1.0+(gamma-1.0)*s*s*s*(10.0+3.0*s*(2.0*s-5.0)));
}

double Len2( double a, double tau, double Tt){
	
	double gamma=1.0;
	double s=tau/Tt;
	return a*(gamma)*s+a*(1-s);
}

double Lenprime( double a, double tau, double Tt){
	
	double gamma=1;
	double s=tau/Tt;
	return 30*a*(gamma-1.0)*(s*s+s*s*s*s-2*s*s*s)/Tt;
}




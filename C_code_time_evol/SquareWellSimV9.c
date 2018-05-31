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
	double Tt=0.1;
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
	*xmax = 1;
	*xmin = 0;
	*V0 = 0;
	*a = 1.0; /* initial size of the well*/
	*k1= PI/(*a);
	*dx = 0.01;
	*dx2 = (*dx)*(*dx);
	*n = (*xmax - *xmin)/(*dx);
	*nt = 100000;
	*dt=0.01/(*nt);
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
	double *Uhalfprime;
	int n_x=(int)(Len(a, Tt, Tt)/dx);
	int n_y=(int)(Len(a, Tt, Tt)/dx);
	int pos;
	int i,k;
	int init=(int)((-xmin)/dx );
	int fini=(int)((Len(a, Tt, Tt)-xmin)/dx);
	double complex sum;
	
	
	
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/* setting up the time evolution matrix that will be invrted for the  crank-nicolson method*/

	Uhalf=malloc(n*n*sizeof(double));
	Uhalfprime=malloc(n*n*sizeof(double));
	
	/*first everything is zero*/
	for (int i=0; i < n; i++) {
		for (int k=0; k < n; k++) {
			
			pos = k + (n * i) ;  /*position in the array*/
			
			Uhalf[pos]=0.0;
		}
		
		
	}
	
	/*submatrix of elements that correspond to the subspace that will be evolved*/
	pos = init + (n_y * init);
	Uhalf[pos]=1.0+I*dt*(I*1.0/dx2);
	pos = init + 1 + (n_y * init);
	Uhalf[pos]=-I*dt*(I*0.5/dx2);
	pos = init + (n_y * (init+1));
	Uhalf[pos]=-I*dt*(I*0.5/dx2);

	
	pos = fini-1 + (n_y * (fini-1));
	Uhalf[pos]=1.0+I*dt*(I*1.0/dx2);
	pos = fini-2 + (n_y * (fini-1));
	Uhalf[pos]=-I*(I*0.5*dt/dx2);
	pos = fini-1 + (n_y * (fini-2));
	Uhalf[pos]=-I*(I*0.5*dt/dx2);
	
	for (int i=init+1; i < fini-1; i++) {
		double x = (i-1)*dx;
		
		for (int k=init+1; k < fini-1; k++) {
			
			pos = k + (n_y * i);  /*position in the array*/
			
			if(i==k){
				Uhalf[pos]=1.0+I*dt*(I*1.0/dx2);
			}
			if(i==(k+1)){
				Uhalf[pos]=-I*dt*(I*0.5/dx2);
			}
			if(i==(k-1)){
				Uhalf[pos]=-I*dt*(I*0.5/dx2);
			}
			
		}
		
		
	}
	
	
	
	 
	 
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
	
	
	
	
	
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/* inverting the unitary operator with half a time step*/
	
	for (int i=0; i < n; i++) {
		sum=0.0;
		for (int k=0; k < n; k++) {
			pos = i + (n * k) ;
			sum+=Uhalf[pos]*Psiold[k];
			
		}
		Psi[i]=sum;
		
	}
	
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
	
	
	
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/* performing time evolution step*/
	
	
	
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
	/*
	sum=0;
	for (int k=0; k < n; k++) {
		sum=Psi[k]-Psiold[k];
		printf("%f\n",creal(sum*conj(sum)));
	}
	*/

	
	
	
	
	
	
	
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
	/*return 0;*/
	
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








/*
	double complex z1 = 1.0 + 3.0 * I;
	double complex z2 = 1.0 - 4.0 * I;
	
	printf("Working with complex numbers:\n\v");
	
	printf("Starting values: Z1 = %.2f + %.2fi\tZ2 = %.2f %+.2fi\n", creal(z1), cimag(z1), creal(z2), cimag(z2));
	double complex sum = z1 + z2;
	printf("The sum: Z1 + Z2 = %.2f %+.2fi\n", creal(sum), cimag(sum));
	
	double complex difference = z1 - z2;
	printf("The difference: Z1 - Z2 = %.2f %+.2fi\n", creal(difference), cimag(difference));
	
	double complex product = z1 * z2;
	printf("The product: Z1 x Z2 = %.2f %+.2fi\n", creal(product), cimag(product));
	
	double complex quotient = z1 / z2;
	printf("The quotient: Z1 / Z2 = %.2f %+.2fi\n", creal(quotient), cimag(quotient));
	
	double complex conjugate = conj(z1);
	printf("The conjugate of Z1 = %.2f %+.2fi\n", creal(conjugate), cimag(conjugate));
	*/



/*
	
 for (int i=0; i < n; i++) {
 for (int k=0; k < n; k++) {
 
 pos = k + (n * i) ;
 
 printf("%f \n",creal(Uhalf[pos]));
 
 }
 printf("\n");
 
 }*/

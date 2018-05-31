// motion of a wavepacket incident on a potential
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


const double PI = 3.14159265358979323846264338327950;

void parameters(double *x0, double *k0, double *width, double *V0, double *a,
				double *xmin, double *xmax, double *n, double *dx, double *dx2,
				double *dt, double *k1, double *nt);
void initial_packet(double Re[], double Im[], double x0, double k0,
					double width, double xmin, double n, double dx, double dt);
void initial_packet2(double Re[], double Im[], double x0, double k1,
					double a, double xmin, double n, double dx, double dt);
void initial_packet3(double Re[], double Im[], double x0, double k0,
					double width, double xmin, double n, double dx, double dt);
void evolve(double Re[], double Im[], double Imold[], double *t,
			double V0, double a, double dx, double dx2, double dt,
			double xmin, double n, double nt, int j);

double V(double x, double V0, double a, double tau, double Tt);
double Len2( double a, double tau, double Tt);
double Len( double a, double tau, double Tt);

void evolve2(double Re[], double Im[], double Imold[], double *t,
			 double V0, double a, double dx, double dx2, double dt,
			 double xmin, double n, double Tt, int j);
double Omegsq( double a, double tau, double Tt);

int main(void)
{
	double magnitud=0.0;
	int i,j;
	double *Re;
	double *Im;
	double *Norm;
	double Tt=1.0;
	FILE *in;
	int var;
	int test;
	char filename[100]="packet2.dat";
	
	
	
	
	
	double x0, k0, width, V0, a, xmin, xmax, n, dx, dx2,nt, dt,k1;
	parameters(&x0, &k0, &width, &V0, &a, &xmin, &xmax, &n, &dx, &dx2, &dt,&k1,&nt);
	
	Re=malloc(n*sizeof(double));
	Im=malloc(n*sizeof(double));
	Norm=malloc(n*sizeof(double));
	
	
	/*initial_packet(Re, Im, x0, k0, width, xmin, n, dx, dt);*/
	initial_packet2(Re, Im, a/2.0, k1, a, xmin, n, dx, dt);
	
	/*initial_packet3(Re, Im, x0, k0, width, xmin, n, dx, dt);*/
	
	
	double t = 0.0;
	for(j=0;j<nt;j++){
		evolve(Re, Im,  Norm, &t, V0, a, dx, dx2, dt, xmin, n, Tt,j);
	
	}
	
	for(j=0;j<0*nt;j++){
		evolve2(Re, Im, Norm, &t, V0, a, dx, dx2, dt, xmin, n, Tt,j);
		
	}
	
	
	
	/*opens the file, writes, closes the file*/
	/*printf("Writing to file: %s\n", filename);*/
	in = fopen(filename,"w");
	if(!in){
		printf("problems opening the file %s\n", filename);
		exit(1);
	}
	
	for(i=0;i<n;i++){
		magnitud=Re[i]*Re[i]+Im[i]*Im[i];
		fprintf(in," %f %f %f %f\n",xmin+i*dx,magnitud,Re[i],Im[i]);
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
	*xmax = 5;
	*xmin = -*xmax;
	*V0 = 0;
	*a = 1.0; /* initial size of the well*/
	*k1= PI/(*a);
	*dx = 0.01;
	*dx2 = (*dx)*(*dx);
	*n = (*xmax - *xmin)/(*dx);
	*nt = 100000;
	*dt=1.0/(*nt);
}

void initial_packet(double Re[], double Im[], double x0, double k0,
					double width, double xmin, double n, double dx, double dt)
// initial Gaussian wavepacket
{
	
	double delta2 = width*width;
	double A = pow(2*PI*delta2, -0.25);
	/*double b = 0.5*k0*dt;*/
	for (int i=1; i <= n; i++) {
		double x = xmin + (i-1)*dx;
		double arg = 0.25*pow(x - x0, 2)/delta2;
		double e = exp(-arg);
		Re[i] = A*cos(k0*(x - x0))*e;
		Im[i] = A*sin(k0*(x - x0))*e;
		/*Im[i] = A*sin(k0*(x - x0 - 0.5*b))*e;*/
	}
}

void initial_packet2(double Re[], double Im[], double x0, double k1,
					double a, double xmin, double n, double dx, double dt)
// initial sine wavepacket
{

	double A = pow(2.0/a, 0.5);
	/*double b = 0.5*k0*dt;*/
	for (int i=1; i <= n; i++) {
		double x = xmin + (i-1)*dx;
		Re[i] = 0.0;
		Im[i] = 0.0;
		
	}
	for (int i=(int)(-xmin)/dx ; i <= (int)(a-xmin)/dx; i++) {
		double x = (i-1)*dx;
		/*Re[i] = A*sin(k1*(x - a))*cos((100)*(x - x0));
		Im[i] = A*sin(k1*(x - a))*sin((100)*(x - x0));*/
		Re[i] = A*sin(k1*(x));
		Im[i] = 0.0;
		/*Re[i] = A*sin(k1*(x - a));
		Im[i] =0.0;*/

	}
}

void initial_packet3(double Re[], double Im[], double x0, double k0,
					double width, double xmin, double n, double dx, double dt)
// initial Gaussian wavepacket
{
	
	double delta2 = width*width;
	double A = pow(2*PI*delta2, -0.25);
	/*double b = 0.5*k0*dt;*/
	for (int i=1; i <= n; i++) {
		double x = xmin + (i-1)*dx;
		Re[i] = 0.0;
		Im[i] = 0.0;
		
	}
	for (int i=(int)(-xmin)/dx ; i <= (int)(1.0-xmin)/dx; i++) {
		double x = (i-1)*dx;
		double arg = 0.25*pow(x , 2)/delta2;
		double e = exp(-arg);
		Re[i] = A*cos(k0*(x - x0))*e;
		Im[i] = A*sin(k0*(x - x0))*e;
		/*Im[i] = A*sin(k0*(x - x0 - 0.5*b))*e;*/
	}
}

/* np sirve con funciones de onda discontinuas, ni deberia*/
void evolve(double Re[], double Im[], double Imold[], double *t,
			double V0, double a, double dx, double dx2, double dt,
			double xmin, double n, double Tt, int j)
{
	double tau=*t;
	
	for (int i=(int)((-xmin)/dx +1); i <= (int)((Len(a, tau, Tt)-xmin)/dx -1 ); i++) {
		double x = (i-1)*dx;
		double HRe = V(x,V0,a,tau,Tt)*Re[i] - 0.5*(Re[i+1] - 2*Re[i] + Re[i-1])
		/ (dx2);
		Im[i] += HRe*dt;
	}
	
	
	
	for (int i=(int)((-xmin)/dx +1) ; i <= (int)((Len(a, tau, Tt)-xmin)/dx -1 ); i++) {
		double x = (i-1)*dx;
		double HIm = V(x,V0,a,tau,Tt)*Im[i] - 0.5*(Im[i+1] - 2*Im[i] + Im[i-1])
		/(dx2);
		Re[i] -= HIm*dt;
	}
	
	for (int i=(int)((Len(a, tau, Tt)-xmin)/dx +1 );i<n; i++) {
		Re[i] =0;
		Im[i] =0;
	}
	

	/*integration*/
	double sum =0;
	int i;
	for(i = 0;i < n;i++){
		sum += (Re[i]*Re[i]+Im[i]*Im[i])*dx;
	}
	/*normalization */
	for(i = 0;i < n;i++){
		Re[i]/= sqrt(sum);
		Im[i]/= sqrt(sum);
	}
/*	printf("%f %f %f \n",tau, Omegsq( a, tau, Tt),Len(a, tau, Tt));*/
	
	*t += dt;        // time of real part
}




void evolve2(double Re[], double Im[], double Imold[], double *t,
			double V0, double a, double dx, double dx2, double dt,
			double xmin, double n, double Tt, int j)
{
	double tau=*t;
	for (int i=(int)(-xmin)/dx +1; i <= (int)((Len(a, Tt, Tt)-xmin)/dx-1); i++) {
		double x = (i-1)*dx;
		double HRe =  - 0.5*(Re[i+1] - 2*Re[i] + Re[i-1])
		/ dx2;
		Im[i] += HRe*dt;
	}
	
	for (int i=(int)(-xmin)/dx +1; i <= (int)((Len(a, Tt, Tt)-xmin)/dx-1); i++) {
		double x = (i-1)*dx;
		double HIm = - 0.5*(Im[i+1] - 2*Im[i] + Im[i-1])
		/ dx2;
		Re[i] -= HIm*dt;
	}
	

	*t += dt;        // time of real part
	/*integration*/
	double sum =0;
	int i;
	for(i = 0;i < n;i++){
		sum += (Re[i]*Re[i]+Im[i]*Im[i])*dx;
	}
	/*normalization */
	for(i = 0;i < n;i++){
		Re[i]/= sqrt(sum);
		Im[i]/= sqrt(sum);
	}
}


double V(double x, double V0, double a, double tau, double Tt)
{

	
	return 0.15*x*x*Omegsq( a, tau, Tt)/2.0;
	/*return 0;*/
	
}

double Omegsq( double a, double tau, double Tt)
{
	double gamma=1.0;
	double s=tau/Tt;
	double gammadprime=(gamma-1)*(60.0*s*(1.0+s*(2.0*s-3.0)))/(Tt*Tt);

	return -gammadprime/Len(1.0, tau, Tt);
	/*return 0;*/
	
}

double Len( double a, double tau, double Tt){

	double gamma=4.0;
	double s=tau/Tt;
	return a*(1.0+(gamma-1.0)*s*s*s*(10.0+3.0*s*(2.0*s-5.0)));
	}

double Len2( double a, double tau, double Tt){
	
	double gamma=5.0;
	double s=tau/Tt;
	return a*(gamma)*s+a*(1-s);
}




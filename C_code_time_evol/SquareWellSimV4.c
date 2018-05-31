// motion of a wavepacket incident on a potential
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


const double PI = 3.14159265358979323846264338327950;

void parameters(double *x0, double *k0, double *width, double *V0, double *a,
				double *xmin, double *xmax, double *n, double *dx, double *dx2,
				double *dt, double *k1);
void initial_packet(double Re[], double Im[], double x0, double k0,
					double width, double xmin, double n, double dx, double dt);
void initial_packet2(double Re[], double Im[], double x0, double k1,
					double a, double xmin, double n, double dx, double dt);
void evolve(double Re[], double Im[], double Imold[], double *t,
			double V0, double a, double dx, double dx2, double dt,
			double xmin, double n, double Tt, int j);
double V(double x, double V0, double a, double j, double Tt);
double V(double x, double V0, double a, double j, double Tt);

int main(void)
{
	double magnitud=0.0;
	int i,j;
	double *Re;
	double *Im;
	double *Imold;
	double Tt=10000.0;
	int T=10000;
	/*double Tt=100000.0;
	int T=100000;*/
	
	double x0, k0, width, V0, a, xmin, xmax, n, dx, dx2, dt,k1;
	parameters(&x0, &k0, &width, &V0, &a, &xmin, &xmax, &n, &dx, &dx2, &dt,&k1);
	
	Re=malloc(n*sizeof(double));
	Im=malloc(n*sizeof(double));
	Imold=malloc(n*sizeof(double));
	
	
	/*initial_packet(Re, Im, x0, k0, width, xmin, n, dx, dt);*/
	initial_packet2(Re, Im, a/2.0, k1, a, xmin, n, dx, dt);
	
	double t = 0.0;
	for(j=0;j<T;j++){
		evolve(Re, Im, Imold, &t, V0, a, dx, dx2, dt, xmin, n, Tt,j);
	
	}
	for(i=0;i<n;i++){
		magnitud=Re[i]*Re[i]+Im[i]*Im[i];
		printf(" %f %f \n",xmin+i*dx,magnitud);
	}
	

	return EXIT_SUCCESS;
}

void parameters(double *x0, double *k0, double *width, double *V0, double *a,
				double *xmin, double *xmax, double *n, double *dx, double *dx2,
				double *dt, double *k1)
{
	*x0 = 0;
	*width = 20;
	*k0 = 4;
	*xmax = 50;
	*xmin = -*xmax;
	*V0 = 0;
	*a = 5.0; /* initial size of the well*/
	*k1= PI/(*a);
	*dx = 0.01;
	*dx2 = (*dx)*(*dx);
	*n = (*xmax - *xmin)/(*dx);
	*dt = 0.0001;
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
		/*Re[i] = A*sin(k1*(x - a))*cos((10)*(x - x0));
		Im[i] = A*sin(k1*(x - a))*sin((10)*(x - x0));*/
		Re[i] = A*sin(k1*(x - a));
		Im[i] =0.0;

	}
}


/* np sirve con funciones de onda discontinuas, ni deberia*/
void evolve(double Re[], double Im[], double Imold[], double *t,
			double V0, double a, double dx, double dx2, double dt,
			double xmin, double n, double Tt, int j)
{
	
	for (int i=1; i < n; i++) {
		double x = xmin + (i-1)*dx;
		double HIm = V(x,V0,a,(double)j,Tt)*Im[i] - 0.5*(Im[i+1] - 2*Im[i] + Im[i-1])
		/ dx2;
		// real part defined at multiples of dt
		Re[i] += HIm*dt;
	}
	
	for (int i=1; i < n; i++) {
		double x = xmin + (i-1)*dx;
		// dt/2 earlier than real part
		Imold[i] = Im[i];
		double HRe = V(x,V0,a,(double)j,Tt)*Re[i] - 0.5*(Re[i+1] - 2*Re[i] + Re[i-1])
		/ dx2;
		// dt/2 later than real part
		Im[i] -= HRe*dt;
	}
	*t += dt;        // time of real part
}

double V(double x, double V0, double a, double j, double Tt)
{
	double gamma=5;
	// step potential
	/*if ( x < a*(1+(gamma-1)*(j/Tt)*(j/Tt)*(j/Tt)*(10+3*(j/Tt)*(2*(j/Tt)-5)))) {*/
	if ( fabs(x) < a) {
	 
	/*if ( fabs(x) < a+(gamma-1)*a*j/Tt) {*/
		return 0;/*-(1E+7)*x*x*(gamma-1)*(60*(j/Tt)*(1+(j/Tt)*(2*(j/Tt)-3)))/(2*Tt*Tt*(1+(gamma-1)*(j/Tt)*(j/Tt)*(j/Tt)*(10+3*(j/Tt)*(2*(j/Tt)-5)))); */
	}
	else{
		return 0;
	}
}

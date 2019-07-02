#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>    /* Standard Library of Complex Numbers */


#define MIN(x,y) x<y?x:y

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


double ran1(long *idum);
double gasdev(long *idum);
void write_to_file(double *arr, int size);
double *init_array(int n_steps, double val);

double *linspace( double x0, double xf, int n_steps);

float *load_matrix(char *filename, int *n, int *m);
double *init_matrix( int n_row, int n_cols, double val);
void print_array(double *x, int n_steps);
void print_array_int(int *x, int n_steps);
void print_matrix(double *A,  int n_rowA, int n_colsA);
double gammln(double xx);
double factrl(int n);
double factrl_doub(double n);
void matrix_multiply(double *matrix_C ,double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
void matrix_multiply_doub_int(double *matrix_C ,double *A,int *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
void matrix_multiply_int_doub(double *matrix_C ,int *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
double sum_array(double *x, int N);
double std_array(double *x, int N);
void std_array_along(double *matrix_C ,double *A, int axis, int n_rowA, int n_colsA);
void sum_array_along(double *matrix_C ,double *A, int axis, int n_rowA, int n_colsA);
void array_multiply(double *matrix_C ,double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
void array_multiply_doub_int(double *matrix_C ,double *A, int *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
void array_divide(double *matrix_C ,double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
void array_sum(double *matrix_C ,double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
void array_log(double *matrix_C ,double *A,int n_rowA, int n_colsA);
void array_exp(double *matrix_C ,double *A,int n_rowA, int n_colsA);
void array_factrl(double *matrix_C ,int *A,int n_rowA, int n_colsA);

void array_sin(double *matrix_C ,double *A,int n_rowA, int n_colsA);
void array_cos(double *matrix_C ,double *A,int n_rowA, int n_colsA);
void array_Window(double *matrix_C, double *x, double d, double cent ,int *A,int n_rowA, int n_colsA);
void array_sign(double *matrix_C, double *x, int *A,int n_rowA, int n_colsA);

void array_factrl_doub(double *matrix_C ,double *A,int n_rowA, int n_colsA);
void scalar_sum(double *matrix_C ,double a, double *A,int n_rowA, int n_colsA);
void scalar_multiply(double *matrix_C ,double a, double *A,int n_rowA, int n_colsA);
void array_normal(double *matrix_C ,long *idum ,int n_rowA, int n_colsA);
void array_assign(double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
void array_assign_floa(double *A,float *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);
void diag(double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB);


void V(double *Vxt, double *temp1, double *x, double DD, double t, double tau );
/*inner declaration of d1 d2 cent1 cent2*/
void Psi_init(double *psi,  double enerlev, double L, double x0);
void Band_Mult(double complex *U_diag, double complex *C_diag, double complex *L_diag,double complex *vector);
void TDMASolve(double complex *L_diag, double complex *C_diag, double complex *U_diag, double complex *vector);





int main(int argc, char **argv){
	
	
	int i,j;
	int x_points, t_points, tau_points;
	double tau,t,T;
	double DD=-20;
	int enerlev;
	double dx,dt;
	double *x;
	double *Vxt;
	double *temp1;
	double *psi;
	
	double L;
	double x0,xf;
	
	x_points=atof(argv[1]);
	T=atof(argv[2]);
	enerlev=atof(argv[3]);
	
	
	/*variables for use of random number geneator*/
	long *idum;
	long ed=-50;
	idum = &ed;
	/*******************/
	
	
	double complex *z1;
	z1=malloc(20*sizeof(double complex));
	
	
	
	return 0;
}






double ran1(long *idum)
/*
 “Minimal” random number generator of Park and Miller with Bays-Durham shuffle and added
 safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter idum between
 successive deviates in a sequence. RNMX should approximate the largest floating value that is less than 1.*/
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;
	if (*idum <= 0 || !iy) { /*Initialize.*/
		if (-(*idum) < 1) *idum=1; /*Be sure to prevent idum = 0.*/
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) { /*Load the shuffle table (after 8 warm-ups).*/
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ; /*Start here when not initializing.*/
	*idum=IA*(*idum-k*IQ)-IR*k; /*Compute idum=(IA*idum) % IM without overflows*/
	if(*idum < 0) *idum += IM; /* by Schrage’s method.*/
	j=iy/NDIV; /*Will be in the range 0..NTAB-1.*/
	iy=iv[j]; /*Output previously stored value and refill the shuffle table.*/
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX; /*Because users don’t expect endpoint values.*/
	else return temp;
}

double gasdev(long *idum)
/*Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum)
 as the source of uniform deviates.*/
{
	double ran1(long *idum);
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;
	if (*idum < 0) iset=0; /*Reinitialize.*/
	if (iset == 0) { /*We don’t have an extra deviate handy, so*/
		do {
			v1=2.0*ran1(idum)-1.0; /*pick two uniform numbers in the square extending from -1 to +1 in each direction,*/
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2; /*see if they are in the unit circle, and if they are not, try again.*/
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		/*Now make the Box-Muller transformation to get two normal deviates. Return one and
		 save the other for next time.*/
		gset=v1*fac;
		iset=1; /*Set flag.*/
		return v2*fac;
	} else { /*We have an extra deviate handy,*/
		iset=0; /*so unset the flag,*/
		return gset; /*and return it.*/
	}
}

double * init_array(int n_steps,double val){
	int i;
	double *x;
	if(!(x=malloc(sizeof(double) * n_steps))){
		fprintf(stderr,"problem with malloc\n");
		exit(1);
	}
	for(i=0;i<n_steps;i++){
		x[i] = val;
	}
	return x;
}

float *load_matrix(char *filename, int *n, int *m){
	float *matrix;
	FILE *in;
	int n_row, n_cols;
	int i;
	int j;
	
	if(!(in=fopen(filename, "r"))){
		printf("Problem opening file %s\n", filename);
		exit(1);
	}
	
	fscanf(in, "%d %d\n", &n_row, &n_cols);
	/*printf("%d %d\n", n_row, n_cols);*/
	
	matrix = malloc(n_row * n_cols * sizeof(float));
	
	for(i=0;i<n_row;i++){
		for(j=0;j<n_cols;j++){
			fscanf(in, "%f", &matrix[i*n_cols + j]);
		}
	}
	*n = n_row;
	*m = n_cols;
	return matrix;
}


double * init_matrix( int n_row, int n_cols, double val){
	
	double *matrix;
	
	int i;
	int j;
	
	
	matrix = malloc(n_row * n_cols * sizeof(double));
	
	for(i=0;i<n_row;i++){
		for(j=0;j<n_cols;j++){
			matrix[i*n_cols + j]=val;
		}
	}
	
	return matrix;
	
}


void print_array(double *x, int n_steps){
	int i;
	for(i=0;i<n_steps-1;i++){
		fprintf(stdout, "%f,", x[i]);
	}
	fprintf(stdout, "%f\n", x[n_steps-1]);
}

void print_array_int(int *x, int n_steps){
	int i;
	for(i=0;i<n_steps-1;i++){
		fprintf(stdout, "%d,", x[i]);
	}
	fprintf(stdout, "%d\n", x[n_steps-1]);
}


void print_matrix(double *A,  int n_rowA, int n_colsA){
	int i,j;
	for(i=0;i<n_rowA;i++){
		for(j=0;j<n_colsA-1;j++){
			fprintf(stdout, "%f,", A[i*n_colsA + j]);
		}
		fprintf(stdout, "%f\n", A[i*n_colsA + (n_colsA -1)]);
	}
}



double gammln(double xx)
/*Returns the value ln[gamma(xx)] for xx > 0.*/
{
	/*Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
	 accuracy is good enough.*/
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

double factrl(int n)
/*Returns the value n! as a floating-point number.*/
{
	double gammln(double xx);
	void nrerror(char error_text[]);
	static int ntop=4;
	static double a[33]={1.0,1.0,2.0,6.0,24.0}; /*Fill in table only as required.*/
	int j;
	if (n < 0) printf("Negative factorial in routine factrl");
	if (n > 32) return exp(gammln(n+1.0));
	/*Larger value than size of table is required. Actually, this big a value is going to overflow
	 on many computers, but no harm in trying.*/
	while (ntop<n) { /*Fill in table up to desired value.*/
		j=ntop++;
		a[ntop]=a[j]*ntop;
	}
	return a[n];
}

double factrl_doub(double n)
/*Returns the value n! as a floating-point number.*/
{
	int m;
	m= (int)n;
	
	double gammln(double xx);
	void nrerror(char error_text[]);
	static int ntop=4;
	static double a[33]={1.0,1.0,2.0,6.0,24.0}; /*Fill in table only as required.*/
	int j;
	if (m < 0) printf("Negative factorial in routine factrl");
	if (m > 32) return exp(gammln(m+1.0));
	/*Larger value than size of table is required. Actually, this big a value is going to overflow
	 on many computers, but no harm in trying.*/
	while (ntop<m) { /*Fill in table up to desired value.*/
		j=ntop++;
		a[ntop]=a[j]*ntop;
	}
	return a[m];
}


void matrix_multiply(double *matrix_C ,double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB)
{
	int i, j,k;
	
	
	for(i=0;i<n_rowA;i++){
		for(j=0;j<n_colsB;j++){
			matrix_C[i*n_colsB + j]=0.0;
			for(k=0;k<n_colsA;k++){
				matrix_C[i*n_colsB + j] = matrix_C[i*n_colsB + j]+A[i*n_colsA + k]*B[k*n_colsB + j];
			}
		}
		
	}
	
}

void matrix_multiply_doub_int(double *matrix_C ,double *A,int *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB)
{
	int i, j,k;
	
	for(i=0;i<n_rowA;i++){
		for(j=0;j<n_colsB;j++){
			matrix_C[i*n_colsB + j]=0.0;
			for(k=0;k<n_colsA;k++){
				matrix_C[i*n_colsB + j] = matrix_C[i*n_colsB + j]+A[i*n_colsA + k]*B[k*n_colsB + j];
			}
			
		}
		
	}
	
	
}

void matrix_multiply_int_doub(double *matrix_C ,int *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB)
{
	int i, j,k;
	
	for(i=0;i<n_rowA;i++){
		for(j=0;j<n_colsB;j++){
			matrix_C[i*n_colsB + j]=0.0;
			for(k=0;k<n_colsA;k++){
				matrix_C[i*n_colsB + j] = matrix_C[i*n_colsB + j]+A[i*n_colsA + k]*B[k*n_colsB + j];
			}
		}
	}
}

double sum_array(double *x, int N){
	int i;
	double Sum=0.0;
	for (i=0; i<N; i++) {
		Sum+=x[i];
	}
	return Sum;
}

double std_array(double *x, int N)
{
	double sum = 0.0, mean, tmp = 0.0;
	int i;
	
	for(i=0; i<N; ++i)
	{
		sum += x[i];
	}
	mean = sum/N;
	
	for(i=0; i<N; ++i)
	tmp += pow(x[i] - mean, 2);
	
	return sqrt(tmp/N);
}


void std_array_along(double *matrix_C ,double *A, int axis, int n_rowA, int n_colsA)
{
	double sum = 0.0, mean, tmp = 0.0;
	int i,j;
	
	if(axis==0){
		/*matrix_C = malloc(n_colsA* sizeof(double));*/
		for(j=0;j<n_colsA;j++){
			
			matrix_C[j]=0.0;
			mean=0.0;
			
			for(i=0;i<n_rowA;i++){
				mean=matrix_C[j]+A[i*n_colsA + j];
			}
			mean=mean/n_rowA;
			
			for(i=0; i<n_rowA; ++i){
				tmp += pow(A[i*n_colsA + j] - mean, 2);
			}
			matrix_C[j]=sqrt(tmp/n_rowA);
		}
		
	}
	else if(axis==1){
		/*matrix_C = malloc(n_rowA * sizeof(double));*/
		for(i=0;i<n_rowA;i++){
			
			matrix_C[i]=0.0;
			mean=0.0;
			
			for(j=0;j<n_colsA;j++){
				matrix_C[i]=matrix_C[i]+A[i*n_colsA + j];
			}
			mean=mean/n_colsA;
			
			for(j=0; j<n_colsA; ++j){
				tmp += pow(A[i*n_colsA + j] - mean, 2);
			}
			
			matrix_C[i]=sqrt(tmp/n_colsA);
		}
		
	}
	else{
		printf("axis not valid just 1 or 0 allowed \n");
		exit(1);
		
	}
}

void  sum_array_along(double *matrix_C ,double *A, int axis, int n_rowA, int n_colsA){
	
	int i,j;
	
	if(axis==0){
		/*matrix_C = malloc(n_colsA* sizeof(double));*/
		for(j=0;j<n_colsA;j++){
			matrix_C[j]=0.0;
			for(i=0;i<n_rowA;i++){
				matrix_C[j]=matrix_C[j]+A[i*n_colsA + j];
			}
		}
		
	}
	else if(axis==1){
		/*matrix_C = malloc(n_rowA * sizeof(double));*/
		for(i=0;i<n_rowA;i++){
			matrix_C[i]=0.0;
			for(j=0;j<n_colsA;j++){
				matrix_C[i]=matrix_C[i]+A[i*n_colsA + j];
			}
		}
		
	}
	else{
		printf("axis not valid just 1 or 0 allowed \n");
		exit(1);
		
	}
}

void array_multiply(double *matrix_C ,double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB)
{
	int i,size;
	
	
	if(n_rowA!=n_rowB ||n_colsA!=n_colsB){
		printf("dimensions do not coincide \n");
		exit(1);
	}
	else{
		
		size=n_rowA*n_colsA;
		for( i=0;i<size;i++){
			matrix_C[i]=A[i]*B[i];
		}
		
	}
	
}

void array_divide(double *matrix_C ,double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB)
{
	int i,size;
	
	
	if(n_rowA!=n_rowB ||n_colsA!=n_colsB){
		printf("dimensions do not coincide \n");
		exit(1);
	}
	else{
		
		size=n_rowA*n_colsA;
		for( i=0;i<size;i++){
			matrix_C[i]=A[i]/B[i];
		}
		
	}
	
	
}


void array_sum(double *matrix_C ,double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB){
	
	int i,size;
	
	if(n_rowA!=n_rowB || n_colsA!=n_colsB){
		printf("dimensions do not coincide \n");
		exit(1);
	}
	else{
		
		size=n_rowA*n_colsA;
		for( i=0;i<size;i++){
			matrix_C[i]=A[i]+B[i];
		}
		
	}
	
	
}
void array_log(double *matrix_C ,double *A, int n_rowA, int n_colsA){
	
	int i,size;
	
	size=n_rowA*n_colsA;
	
	for( i=0;i<size;i++){
		matrix_C[i]=log(A[i]);
	}
	
}

void array_exp(double *matrix_C,double *A,int n_rowA, int n_colsA)
{
	int i,size;
	
	size=n_rowA*n_colsA;
	
	for( i=0;i<size;i++){
		matrix_C[i]=exp(A[i]);
	}
	
}

void array_factrl(double *matrix_C ,int *A, int n_rowA, int n_colsA){
	
	int i,size;
	
	size=n_rowA*n_colsA;
	
	for( i=0;i<size;i++ ){
		matrix_C[i]=factrl(A[i]);
	}
	
}


void array_factrl_doub(double *matrix_C ,double *A, int n_rowA, int n_colsA){
	
	int i,size;
	
	size=n_rowA*n_colsA;
	
	for( i=0;i<size;i++ ){
		matrix_C[i]=factrl_doub(A[i]);
	}
	
}

void scalar_sum(double *matrix_C ,double a, double *A,int n_rowA, int n_colsA){
	
	int i,size;
	
	size=n_rowA*n_colsA;
	
	for( i=0;i<size;i++ ){
		matrix_C[i]=A[i]+a;
	}
	
}

void scalar_multiply(double *matrix_C ,double a, double *A,int n_rowA, int n_colsA){
	
	int i,size;
	
	size=n_rowA*n_colsA;
	
	for( i=0;i<size;i++ ){
		matrix_C[i]=a*A[i];
	}
	
}


void array_normal(double *matrix_C ,long *idum,int n_rowA, int n_colsA){
	
	int i,size;
	double gasdev(long *idum);
	
	size=n_rowA*n_colsA;
	
	for( i=0;i<size;i++ ){
		matrix_C[i]=gasdev(idum);
		
	}
	
}

void array_assign(double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB)
{
	int i,size;
	
	if(n_rowA!=n_rowB || n_colsA!=n_colsB){
		printf("dimensions do not coincide \n");
		exit(1);
	}
	else{
		
		size=n_rowA*n_colsA;
		for( i=0;i<size;i++){
			A[i]=B[i];
		}
		
	}
}

void array_assign_floa(double *A,float *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB)
{
	int i,size;
	
	if(n_rowA!=n_rowB || n_colsA!=n_colsB){
		printf("dimensions do not coincide \n");
		exit(1);
	}
	else{
		
		size=n_rowA*n_colsA;
		for( i=0;i<size;i++){
			A[i]=B[i];
		}
		
	}
}

void diag(double *A,double *B,int n_rowA,int n_colsA,int n_rowB,int n_colsB)
{
	int i,size;
	
	
	
	size=n_rowA*n_colsA;
	for( i=0;i<size;i++){
		A[i]=B[i*n_colsB + i];
	}
	
	
}





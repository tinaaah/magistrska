#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>

#define ST_PONOVITEV 100
#define ST_DELCEV 85

#define fromfile 0


// printanje vektorja
void print_vector(gsl_vector *vector){
    printf("%lf\t%lf\t%lf\n", gsl_vector_get(vector, 0), gsl_vector_get(vector, 1), gsl_vector_get(vector, 2));
}

// nakljucno porazdeljene tocke na sferi
void randdist(gsl_vector **data, int N){
    for(int i=0; i<N; i++){
	double u = (double)rand()/RAND_MAX;
	double v = (double)rand()/RAND_MAX;

	double phi = 2*M_PI*u;
	double theta = acos(2*v-1);

	data[i] = gsl_vector_alloc(3);
	gsl_vector_set(data[i], 0, cos(phi)*sin(theta));
	gsl_vector_set(data[i], 1, sin(phi)*sin(theta));
	gsl_vector_set(data[i], 2, cos(theta));
    }
}

// Prebrani fajli
void loaddata(gsl_vector **podatki, FILE* file, int N){
    // Preberem file in ga shranim kot vektor v gsl_vector
    fgetc(file);
    for(int i=0; i<N; i++){
	podatki[i] = gsl_vector_alloc(3);
	double tmp1, tmp2, tmp3;
	fscanf(file,"%lf %lf %lf", &tmp1, &tmp2, &tmp3);
	gsl_vector_set(podatki[i], 0, tmp1);
	gsl_vector_set(podatki[i], 1, tmp2);
	gsl_vector_set(podatki[i], 2, tmp3);

	//normiran vektor
	double norma = gsl_blas_dnrm2(podatki[i]);
	gsl_blas_dscal(1/norma, podatki[i]);
    }
    fclose(file);
    // Prebran file brao
}

// strukturni faktor
double structure_factor(gsl_vector **data, int N, int l){
    double sum = 0;
    for(int i=0; i<N; i++){
	for(int j=0; j<N; j++){
	    double product;
	    gsl_blas_ddot(data[i], data[j], &product);
	    if(product > 1) product = 1;
	    product = gsl_sf_legendre_Pl(l, product);
	    sum += product;
	}
    }
    return sum/(double)N;
}

// stevilska varianca preko lefendrovih polinomov
double varianca_legendre(gsl_vector **data, double kot, int N){
    double produkt, vsota;
    vsota = 0;
    for(int l=1; l<=100; l++){
	produkt = gsl_sf_legendre_Pl(l+1,cos(kot)) - gsl_sf_legendre_Pl(l-1, cos(kot));
	produkt = structure_factor(data,N,l) * pow(produkt,2)/(2*l+1);
	vsota += produkt;
    }
    return vsota*N/4;
}

// stevilska varianca preko povprecenja, n je stevilo ponovitev, N je stevilo tock
void varianca_stetje(gsl_vector **data, double *m, int n, int N, double kot){
    double produkt;
    gsl_vector *sredisce;
    for(int i=0; i<n; i++){
	randdist(&sredisce, 1);
	for(int j=0; j<N; j++){
	    gsl_blas_ddot(sredisce, data[j], &produkt);
	    if(produkt < cos(kot)) {
		m[i]++;
	    }
	}
    }
}

//printf("%e\n", pow(gsl_vector_get(sredisce,0),2) + pow(gsl_vector_get(sredisce,1),2) + pow(gsl_vector_get(sredisce,2),2));

int main(){
    srand(time(NULL));

    gsl_vector **podatki = malloc(sizeof(gsl_vector*) * ST_DELCEV); // ST_DELCEV = N !!!!!!!

    // Ce imam podatke odkomentiram to:
#if fromfile==1

    FILE *file = fopen("fibonacci/85.txt", "r");
    //FILE *file = fopen("thomson/85.xyz", "r");
    int N;
    fscanf(file, "%d", &N);					    // ST_DELCEV = N !!!!!!

    loaddata(podatki, file, N);

    // Izracuna strukturni faktor za vrednosti l=(1,100)

    /*
    for(int i=1; i<=100; i++){
	double s = structure_factor(podatki, N, i);
	printf("%d\t%e\n", i, s);
    }
	*/


    // Izracuna varianco preko legendrovih polinomov za theta od 0 do pi/2
    double theta, napaka;
    for(int i=0; i<=100; i++){
	theta = i*M_PI/2/100;
	napaka = varianca_legendre(podatki, theta, ST_DELCEV);
	printf("%e\t%e\n", theta/M_PI, napaka);
    }

    /*
    // Izracuna varianco s povprecenjem n ponovitev za theta od 0 do pi/2 
    int n = 10000;
    double napaka;
    int x,x2;
    for(int i=0; i<=100; i++){
	double *so_znotraj = calloc(n,sizeof(double));
	x = 0;
	x2 = 0;
	varianca_stetje(podatki, so_znotraj, n, 85, i*M_PI/100/2); // ST_DELCEV = N !!!!
	for(int j=0; j<n; j++){
	    x += (int)pow(so_znotraj[j],2);
	    x2 += (int)so_znotraj[j];
	}
	//printf("%d\t%d\n", x, x2);
	napaka = ((double)x)/n - pow(((double)x2)/n, 2);
	printf("%e\t%e\n", ((double)i)/100/2, napaka);
    }
	*/
#endif
   
    // Ce imam random odkomentiram to in povprecim
#if fromfile==0

    // Izracuna strukturni faktor za vrednosti l=(1,100)
    /*
    double *stopnja = calloc(100,sizeof(double));
    for(int j=0; j<ST_PONOVITEV; j++){
	randdist(podatki, ST_DELCEV);
	for(int l=1; l<=100; l++){
	    double s = structure_factor(podatki, ST_DELCEV, l);
	    stopnja[l-1] += s;
	}
    }
    for(int l=1; l<=100; l++){
	printf("%d\t%e\n", l, stopnja[l-1]/ST_PONOVITEV);
    }


    // Izracuna varianco preko legendrovih polinomov za theta od 0 do pi/2
    double *napaka = calloc(100,sizeof(double));
    double theta;
    for(int j=0; j<ST_PONOVITEV; j++){
	randdist(podatki, ST_DELCEV);
	for(int i=0; i<=100; i++){
	    theta = i*M_PI/2/100;
	    napaka[i] += varianca_legendre(podatki, theta, ST_DELCEV);
	}
    }
    for(int i=0; i<=100; i++){
	printf("%e\t%e\n", (double)i/100/2, napaka[i]/ST_PONOVITEV);
    }

	*/
    // Izracuna varianco s povprecenjem n ponovitev za theta od 0 do pi/2 

    double *napaka = calloc(100,sizeof(double));
    int x,x2;
    int n = 10000;
    for(int j=0; j<ST_PONOVITEV; j++){
	randdist(podatki, ST_DELCEV);
	for(int i=0; i<=100; i++){
	    double *so_znotraj = calloc(n,sizeof(double));
	    x = 0;
	    x2 = 0;
	    varianca_stetje(podatki, so_znotraj, n, ST_DELCEV, i*M_PI/100/2); // ST_DELCEV = N !!!!
	    for(int j=0; j<n; j++){
		x += (int)pow(so_znotraj[j],2);
		x2 += (int)so_znotraj[j];
	    }
	    napaka[i] += ((double)x)/n - pow(((double)x2)/n, 2);
	}
    }
    for(int i=0; i<=100; i++){
	printf("%e\t%e\n", ((double)i)/100/2, napaka[i]/ST_PONOVITEV);
    }

#endif

}

#define _USE_MATH_DEFINES

#include <stdio.h>
#include <cmath>
#include <complex>

using namespace std;

const complex double im(0.0,1.0);

complex double *matrix_mul(complex double *A, complex double *B, dim3 dimsA, dim3 dimsB, size_t size){
	complex double *C = (complex double *)malloc(size);
	if(dimsA.x != dimsB.y){
		printf("Invalid matrix multiplication\n");
		return -1;
	}
	for(int i = 0; i = dimsA.y; ++i){
		for(int j = 0; j = dimsB.x; ++j){
			for(int k = 0; k = dimsA.x; ++k){
				C[i*dimsA.y+j] += A[i*dimsA.y+k] * B[k*dimsB.x+i];
			}
		}
	}
	return C;

}
double *array_geometry(int n_array, double d){
	//TODO
	//FILE* f_array = fopen(geometry_path, "r");
	double *ar_geometry = (double *)malloc(n_array);
	for(int i=0; i<n_array: ++i){
		ar_geometry[i] = i*d;
	}
	return ar_geometry;
};

complex double *matrix_transpose(complex double *A, int N, int M, size_t size){
	complex double *C = (complex double *)malloc(size);
	for(int i=0; i<N; ++i){
		for(int j=0; j<M; ++j){
			C[i*N+j] = conj(A[j*M+i]);
		}
	}
	return C;
}


complex double *steering_vector(int n_array, int n_sample, double *array_geometry){
	//TODO
	complex double *a = (complex double *)malloc(n_array*n_sample);
	const double theta = 2*M_PI/n_sample;
	for(int i=0; i<n_sample; ++i){ //n_sample by n_array matrix
		for(int j=0; j<n_array; ++j){
			a[i*n_array+j] = exp(-2 * im* M_PI * array_geometry[j] * sin(theta*i));
		}
	}
	return a;
}
complex double *cov_matrix(complex double *src_signal, int n_array, int n_sample){
	//complex double
}

complex double *eigen_wv();
double *music();
double *source_generator(int sample_rate, int n_source);


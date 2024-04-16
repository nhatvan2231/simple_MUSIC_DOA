#include <stdio.h>

double *matrix_mul(double *A, double *B, dim3 dimsA, dim3 dimsB, size_t size){
	double *C = (double *)malloc(size);
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
double *array_geometry(){
	//TODO
	FILE* f_array = fopen(geometry_path, "r");

};

double *matrix_transpose(double *A, int N, int M, size_t size){
	double *C = (double *)malloc(size);
	for(int i=0; i<N; ++i){
		for(int j=0; j<M; ++j){
			C[i*N+j] = A[j*M+i];
		}
	}
	return C;
}


double *steering_vector(int n_array, int n_sample, double *array_geometry){
	//TODO
}
double *eigen_wv();
double *music();
double *source_generator(int sample_rate, int n_source);


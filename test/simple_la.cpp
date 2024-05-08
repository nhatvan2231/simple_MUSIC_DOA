#include "simple_la.h"

using namespace std;

// copy matrix/vector
std::complex<double>* matrix_copy(std::complex<double> *src, int M, int N){
	std::complex<double> *result = (std::complex<double> *)malloc(sizeof(std::complex<double>) * M * N);
	for(int i=0; i<M*N; ++i){
		result[i] = src[i];
	}
	return result;
}

// calculate matrix multiplication
std::complex<double>* matrix_mul(std::complex<double> *srcA,
							std::complex<double> *srcB, 
							int M_A, int N_A, 
							int M_B, int N_B
							){
	std::complex<double> *result = (std::complex<double> *)malloc(sizeof(std::complex<double>) * M_A * N_B);
	if(N_A != M_B){
		printf("Invalid matrix multiplication\n");
		return result;
	}
	for(int i = 0; i < M_A; ++i){
		for(int j = 0; j < N_B; ++j){
			for(int k = 0; k < N_A; ++k){
				result[i*N_B+j] += srcA[i*N_A+k] * srcB[k*N_B+i];
			}
		}
	}
	return result;
}

// calculate the conjugate transpose of a matrix or vector
std::complex<double>* conj_transpose(std::complex<double> *src, int M, int N){
	std::complex<double> *result = (std::complex<double> *)malloc(sizeof(std::complex<double>) * N * M);
	for(int i=0; i<N; ++i){
		for(int j=0; j<M; ++j){
			result[i*M+j] = conj(src[j*N+i]);
		}
	}
	return result;
}

// calculate the covariance matrix without normalize
std::complex<double>* cov_matrix(std::complex<double> *src, int M, int N){
	std::complex<double> *srcT = conj_transpose(src, N, M);
	std::complex<double> *result = matrix_mul(src, srcT, M, N, N, M);
	free(srcT);
	return result;
}

// get the i column of the matrix
std::complex<double>* get_col(std::complex<double>* src, int M, int N, int col_i){
	std::complex<double>* col = (std::complex<double> *)malloc(sizeof(std::complex<double>) * M);
	for(int i=0; i<M; ++i){
		col[i] = src[i*N + col_i];
	}
	return col;
}

// get the identity matrix with size x size
double* identity_mat(int size){
	double* mat_I = (double*)malloc(sizeof(double) * size * size);
	for(int i=0; i<size; ++i){
		for(int j=0; j<size; ++j){
			mat_I[i*size +j] = 0;
			if(i==j) mat_I[i*size +j] = 1;
		}
	}
	return mat_I;
}

// calculate the norm of vector
double vec_norm(std::complex<double> *src, int size){
	double result;
	for(int i=0; i<size; ++i){
		result += norm(src[i]);
	}
	return sqrt(result);
}

// generate random number
double random_num(double min=1.0, double max=10.0){
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(min, max);
	return dis(gen);
}

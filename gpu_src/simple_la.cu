#include "simple_la.cuh"


// copy matrix/vector 
__global__ void matrix_copy(thrust::complex<double> *src, thrust::complex<double>* dst, int M){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idy = blockIdx.y * blockDim.y + threadIdx.y;
	dst[idx*M + idy] = src[idx*M + idy];
}

//// calculate matrix multiplication
//int matrix_mul(std::complex<double> *srcA,
//							std::complex<double> *srcB, 
//							int M_A, int N_A, 
//							int M_B, int N_B,
//							std::complex<double>* &dst){
//	if(N_A != M_B){
//		printf("Invalid matrix multiplication\n");
//		return 1;
//	}
//	for(int i = 0; i < M_A; ++i){
//		for(int j = 0; j < N_B; ++j){
//			for(int k = 0; k < N_A; ++k){
//				dst[i*N_B+j] += srcA[i*N_A+k] * srcB[k*N_B+j];
//			}
//		}
//	}
//	return 0;
//}
//
//// calculate the conjugate transpose of a matrix or vector
//int conj_transpose(std::complex<double> *src, int M, int N, std::complex<double>* &dst){
//	for(int i=0; i<N; ++i){
//		for(int j=0; j<M; ++j){
//			dst[i*M+j] = conj(src[j*N+i]);
//		}
//	}
//	return 0;
//}
//
//// calculate the covariance matrix without normalize
//int cov_matrix(std::complex<double> *src, int M, int N, std::complex<double>* &dst){
//	std::complex<double> *srcT = (std::complex<double>*)malloc(sizeof(std::complex<double>)*N*M);
//	conj_transpose(src, M, N, srcT);
//	matrix_mul(src, srcT, M, N, N, M, dst);
//	free(srcT);
//	return 0;
//}
//
//// get the i column of the matrix
//int get_col(std::complex<double>* src, int M, int N, int col_i, std::complex<double>* &dst){
//	for(int i=0; i<M; ++i){
//		dst[i] = src[i*N + col_i];
//	}
//	return 0;
//}
//
//// get the identity matrix with size x size
//int identity_mat(int size, std::complex<double>* &dst){
//	for(int i=0; i<size; ++i){
//		for(int j=0; j<size; ++j){
//			dst[i*size +j] = std::complex<double>(0,0);
//			if(i==j) dst[i*size +j] = std::complex<double>(1,0);
//		}
//	}
//	return 0;
//}
//
//// calculate the norm of vector
//double vec_norm(std::complex<double> *src, int size){
//	double result=0;
//	for(int i=0; i<size; ++i){
//		result += norm(src[i]);
//	}
//	return sqrt(result);
//}
//
//// generate random number
//double random_num(double min=1.0, double max=10.0){
//	std::random_device rd;
//	std::mt19937 gen(rd());
//	std::uniform_real_distribution<double> dis(min, max);
//	return dis(gen);
//}

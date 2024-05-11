#ifndef _SIMPLE_LA_CUH_
#define _SIMPLE_LA_CUH_
#include <stdio.h>
#include <complex>
#include <iostream>
#include <random>
#include <cuda_runtime.h>
#include <thrust/complex.h>

__global__ void matrix_copy(thrust::complex<double> *src, thrust::complex<double>* dst, int M);
//__global__ int matrix_mul(thrust::complex<double> *srcA,
//				thrust::complex<double> *srcB, 
//				int M_A, int N_A, 
//				int M_B, int N_B, 
//				thrust::complex<double>* &dst);
//__global__ int conj_transpose(thrust::complex<double> *src, int M, int N, thrust::complex<double>* &dst);
//__global__ int cov_matrix(thrust::complex<double> *src, int M, int N, thrust::complex<double>* &dst);
//__global__ int get_col(thrust::complex<double>* src, int M, int N, int col_i, thrust::complex<double>* &dst);
//__global__ int identity_mat(int size, thrust::complex<double>* &src);
//__global__ double vec_norm(thrust::complex<double> *src, int size);
//__global__ double random_num(double min, double max);

#endif

#include <stdio.h>
#include <complex>
#include <iostream>
#include <random>

std::complex<double>* matrix_copy(std::complex<double> *src, int M, int N);
std::complex<double>* matrix_mul(std::complex<double> *srcA,
							std::complex<double> *srcB, 
							int M_A, int N_A, 
							int M_B, int N_B);
std::complex<double>* conj_transpose(std::complex<double> *src, int M, int N);
std::complex<double>* cov_matrix(std::complex<double> *src, int M, int N);
std::complex<double>* get_col(std::complex<double>* src, int M, int N, int col_i);
double* identity_mat(int size);
double vec_norm(std::complex<double> *src, int size);
double random_num(double min, double max);

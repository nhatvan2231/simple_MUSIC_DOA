#ifndef _SIMPLE_LA_H_
#define _SIMPLE_LA_H_
#include <stdio.h>
#include <complex>
#include <iostream>
#include <random>

int matrix_copy(std::complex<double> *src, int M, int N, std::complex<double>* &dst);
int matrix_mul(std::complex<double> *srcA,
							std::complex<double> *srcB, 
							int M_A, int N_A, 
							int M_B, int N_B, 
							std::complex<double>* &dst);
int conj_transpose(std::complex<double> *src, int M, int N, std::complex<double>* &dst);
int cov_matrix(std::complex<double> *src, int M, int N, std::complex<double>* &dst);
int get_col(std::complex<double>* src, int M, int N, int col_i, std::complex<double>* &dst);
int identity_mat(int size, std::complex<double>* &src);
double vec_norm(std::complex<double> *src, int size);
double random_num(double min, double max);

#endif

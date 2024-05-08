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
double vec_norm(std::complex<double> *src, int size);
double random_num(double min=1.0, double max=10.0);

int eigen_power(std::complex<double> *src, int M, int N, double ep, std::complex<double>* &ei_value, std::complex<double>* &ei_vector);
int power_ei_vector(std::complex<double>* A, std::complex<double>* &v, int size, double ep);
int power_ei_value(std::complex<double>* A, std::complex<double>* v, int size, std::complex<double> &ei_value);

int householder_reflections(std::complex<double>* A, int size, std::complex<double>* &R);
int compute_P(std::complex<double>* colA, int size, int col, std::complex<double>* &mat_P);
int compute_n(std::complex<double>* colA, int size, std::complex<double>* &vec_n);
int compute_u(std::complex<double>* colA, int size, std::complex<double>* &vec_u);














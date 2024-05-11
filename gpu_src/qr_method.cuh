//reference source :https://kwokanthony.medium.com/detailed-explanation-with-example-on-qr-decomposition-by-householder-transformation-5e964d7f7656
//
#ifndef _QR_METHOD_H_
#define _QR_METHOD_H_
#include "simple_la.h"

int isRowEchelon(std::complex<double>* A, int size, double eps=1e-9);
int eigen_vv(std::complex<double>* A, int size, std::complex<double>* &ei_vectors, std::complex<double>* &ei_values);
int householder_reflections(std::complex<double>* A, int size, std::complex<double>* &Q, std::complex<double>* &R);
int compute_H(std::complex<double>* colA, int size, int col, std::complex<double>* &mat_H);
int compute_n(std::complex<double>* colA, int size, std::complex<double>* &vec_n);
int compute_u(std::complex<double>* colA, int size, std::complex<double>* &vec_u);
#endif

#include "simple_la.h"

int householder_reflections(std::complex<double>* A, int size, std::complex<double>* &R);
int compute_P(std::complex<double>* colA, int size, int col, std::complex<double>* &mat_P);
int compute_n(std::complex<double>* colA, int size, std::complex<double>* &vec_n);
int compute_u(std::complex<double>* colA, int size, std::complex<double>* &vec_u);
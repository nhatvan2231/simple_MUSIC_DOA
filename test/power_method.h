#include "simple_la.h"

int eigen_power(std::complex<double> *src, int M, int N, double ep, std::complex<double>* &ei_value, std::complex<double>* &ei_vector);
int power_ei_vector(std::complex<double>* A, std::complex<double>* &v, int size, double ep);
int power_ei_value(std::complex<double>* A, std::complex<double>* v, int size, std::complex<double> &ei_value);

#include "power_method.h"

/*****************************************
 *********Power Iteration method**********
*****************************************/
int eigen_power(std::complex<double> *src, int M, int N, double ep, std::complex<double>* &ei_value, std::complex<double>* &ei_vector){
	if (M != N) return -1;
	std::complex<double> *v = (std::complex<double> *)malloc(sizeof(std::complex<double>) * N);
	for(int i=0; i<N; ++i){
		v[i] = std::complex<double>(random_num(), random_num());
	}
	power_ei_vector(src, v, N, ep);
	power_ei_value(src, v, N, ei_vector[0]);

	ei_value = v;
	return 0;
}

int power_ei_vector(std::complex<double>* A, std::complex<double>* &v, int size, double ep){
	// copy the vector 
	std::complex<double> *v_old = (std::complex<double> *)malloc(sizeof(std::complex<double>) * size);
	for(int i=0; i<size; ++i) v_old[i] = v[i];

	std::complex<double>* Av = matrix_mul(A, v, size, size, size, 1);
	double Av_norm = vec_norm(Av, size);

	// estimate eigenvector
	for(int i=0; i<size; ++i) Av[i] /= Av_norm;

	// recursive until certain threshold
	for(int i=0; i<size; ++i){
		if(abs(v_old[i] - Av[i]) > ep){
			power_ei_vector(A, Av, size, ep);
			break;
		}
	}
	v = Av;
	free(v_old);
	return 0;
}

int power_ei_value(std::complex<double>* A, std::complex<double>* v, int size, std::complex<double> &ei_value){
	std::complex<double> e_val=0;
	for(int i=0; i<size; ++i){
		e_val += A[i] * v[i];
	}
	ei_value = e_val/v[0];
	return 0;
}

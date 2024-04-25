#include <stdio.h>
#include <complex>
#include <iostream>
#include <random>

using namespace std;

std::complex<double>* matrix_mul(std::complex<double> *srcA,
							std::complex<double> *srcB, 
							int M_A, int N_A, 
							int M_B, int N_B);
std::complex<double>* conj_transpose(std::complex<double> *src, int M, int N);
std::complex<double>* cov_matrix(std::complex<double> *src, int M, int N);

int eigen_power(std::complex<double> *src, int M, int N, double ep, std::complex<double>* &ei_value, std::complex<double>* &ei_vector);
int power_iter(std::complex<double>* &v, int size, double ep);

double vec_norm(std::complex<double> *src, int size);
double random_num(double min=1.0, double max=10.0);


int main(){

	int M = 5;
	int N = 3;
	int L = 1;

	std::complex<double> *a;
//	std::complex<double> *b(1,2);
	std::complex<double> *b;
	std::complex<double> *c;
	a = (std::complex<double> *)malloc(sizeof(std::complex<double>)*M*N);
	b = (std::complex<double> *)malloc(sizeof(std::complex<double>)*N*L);
	//b = (std::complex<double> *)malloc(sizeof(std::complex<double>)*4);
	for(int i=0; i<M*N; ++i){
		a[i]= std::complex<double>(random_num(),random_num());
	}

	for(int i=0; i<L*N; ++i){
		b[i]= std::complex<double>(random_num(),random_num());
	}


	printf("Matrix Multiplication\n");
	c = matrix_mul(a, b, M, N, N, L);
	for(int i = 0; i < M; ++i) {
		for(int j = 0; j < L; ++j) {
			cout << c[i*M +j] << " ";
		}
		cout << endl;
	}

	printf("Matrix tranpose\n");
	std::complex<double> *cT = conj_transpose(c, M, L);
	for(int i = 0; i < M; ++i) {
		for(int j = 0; j < L; ++j) {
			cout << cT[i*M +j] << " ";
		}
		cout << endl;
	}

	printf("Matrix cov\n");
	std::complex<double> *c_matrix = cov_matrix(c, M, L);
	for(int i = 0; i < M; ++i) {
		for(int j = 0; j < L; ++j) {
			cout << c_matrix[i*M +j] << " ";
		}
		cout << endl;
	}


	std::complex<double>* ei_value;
	std::complex<double>* ei_vector;
	eigen_power(c_matrix, M, M, 0.0001, ei_value, ei_vector);
	for(int i = 0; i < M; ++i) {
		cout << ei_value[i] << endl;
	}

	free(ei_value);
	free(c_matrix);
	free(a);
	free(b);
	free(c);
}

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

std::complex<double>* conj_transpose(std::complex<double> *src, int M, int N){
	std::complex<double> *result = (std::complex<double> *)malloc(sizeof(std::complex<double>) * N * M);
	for(int i=0; i<N; ++i){
		for(int j=0; j<M; ++j){
			result[i*N+j] = conj(src[j*M+i]);
		}
	}
	return result;
}

std::complex<double>* cov_matrix(std::complex<double> *src, int M, int N){
	std::complex<double> *srcT = conj_transpose(src, N, M);
	std::complex<double> *result = matrix_mul(src, srcT, M, N, N, M);
	free(srcT);
	return result;
}

double vec_norm(std::complex<double> *src, int size){
	double result;
	for(int i=0; i<size; ++i){
		result += norm(src[i]);
	}
	return sqrt(result);
}

int eigen_power(std::complex<double> *src, int M, int N, double ep, std::complex<double>* &ei_value, std::complex<double>* &ei_vector){
	if (M != N) return -1;
	std::complex<double> *v = (std::complex<double> *)malloc(sizeof(std::complex<double>) * N);
	for(int i=0; i<N; ++i){
		v[i] = std::complex<double>(random_num(), random_num());
	}
	power_iter(v, N, ep);
	ei_value = v;
	return 0;
}

int power_iter(std::complex<double>* &v, int size, double ep){
	// copy the vector 
	std::complex<double> *v_old = (std::complex<double> *)malloc(sizeof(std::complex<double>) * size);
	for(int i=0; i<size; ++i) v_old[i] = v[i];

	double v_norm = vec_norm(v, size);

	// estimate eigenvector
	for(int i=0; i<size; ++i) v[i] /= v_norm;

	// recursive until certain threshold
	for(int i=0; i<size; ++i){
		if(abs(v_old[i] - v[i]) > ep){
			power_iter(v, size, ep);
			break;
		}
	}
	return 0;
}

double random_num(double min, double max){
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(min, max);
	return dis(gen);
}
































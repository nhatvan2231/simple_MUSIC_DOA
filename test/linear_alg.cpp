#include <stdio.h>
#include <complex>
#include <iostream>

using namespace std;

std::complex<double>* matrix_mul(std::complex<double> *srcA,
							std::complex<double> *srcB, 
							int M_A, int N_A, 
							int M_B, int N_B);
std::complex<double>* conj_transpose(std::complex<double> *src, int M, int N);
std::complex<double>* cov_matrix(std::complex<double> *src, int M, int N);
int main(){
	std::complex<double> *a;
//	std::complex<double> *b(1,2);
	std::complex<double> *b;
	std::complex<double> *c;
	a = (std::complex<double> *)malloc(sizeof(std::complex<double>)*2);
	b = (std::complex<double> *)malloc(sizeof(std::complex<double>)*2);
	c = (std::complex<double> *)malloc(sizeof(std::complex<double>)*2*2);
	//b = (std::complex<double> *)malloc(sizeof(std::complex<double>)*4);
	a[0]= std::complex<double>(1.0,2.0);
	b[0]= std::complex<double>(3.0,4.0);
	a[1]= std::complex<double>(3.0,4.0);
	b[1]= std::complex<double>(3.0,4.0);

	c = matrix_mul(a, b, 2, 1, 1, 2);
	for(int i = 0; i < 2; ++i) {
		for(int j = 0; j < 2; ++j) {
			cout << c[i*2 +j] << " ";
		}
		cout << endl;
	}

	std::complex<double> *cT = conj_transpose(c, 2, 2);
	for(int i = 0; i < 2; ++i) {
		for(int j = 0; j < 2; ++j) {
			cout << cT[i*2 +j] << " ";
		}
		cout << endl;
	}

	std::complex<double> *c_matrix = cov_matrix(c, 2, 2);
	for(int i = 0; i < 2; ++i) {
		for(int j = 0; j < 2; ++j) {
			cout << c_matrix[i*2 +j] << " ";
		}
		cout << endl;
	}

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

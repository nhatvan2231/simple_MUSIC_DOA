//#include "simple_la.h"
#include "qr_method.h"

using namespace std;

int main(){

	int M = 4;
	int N = 4;
	int L = 2;

	std::complex<double> *a;
//	std::complex<double> *b(1,2);
	std::complex<double> *b;
	std::complex<double> *c;
	a = (std::complex<double> *)malloc(sizeof(std::complex<double>)*M*N);
	b = (std::complex<double> *)malloc(sizeof(std::complex<double>)*N*L);
	//b = (std::complex<double> *)malloc(sizeof(std::complex<double>)*4);
	for(int i=0; i<M*N; ++i){
		a[i]= std::complex<double>(random_num(1,2),random_num(1,2));
	}

	for(int i=0; i<L*N; ++i){
		b[i]= std::complex<double>(random_num(1,2),random_num(1,2));
	}


	printf("Matrix Multiplication\n");
	c = matrix_mul(a, b, M, N, N, L);
	for(int i = 0; i < M; ++i) {
		for(int j = 0; j < L; ++j) {
			cout << c[i*L +j] << " ";
		}
		cout << endl;
	}

	printf("Matrix tranpose\n");
	std::complex<double> *cT = conj_transpose(c, M, L);
	for(int i = 0; i < L; ++i) {
		for(int j = 0; j < M; ++j) {
			cout << cT[i*M +j] << " ";
		}
		cout << endl;
	}

	printf("Matrix cov\n");
	std::complex<double> *c_matrix = cov_matrix(c, M, L);
	for(int i = 0; i < M; ++i) {
		for(int j = 0; j < M; ++j) {
			cout << c_matrix[i*M +j] << " ";
		}
		cout << endl;
	}


	std::complex<double>* R = (std::complex<double>*)malloc(sizeof(std::complex<double>)*M*M);
	householder_reflections(c_matrix, M, R);
	printf("R matrix\n");
	for(int i = 0; i < M; ++i) {
		for(int j = 0; j < M; ++j) {
			cout << R[i*M +j] << " ";
		}
		cout << endl;
	}
	//std::complex<double>* ei_value;
	//std::complex<double>* ei_vector;
	//eigen_power(c_matrix, M, M, 0.0001, ei_vector, ei_value);
	//printf("Power Eigen vector\n");
	//for(int i = 0; i < M; ++i) {
	//	cout << ei_vector[i] << endl;
	//}
	//printf("Power Eigen value\n");
	//cout << ei_value[0] << endl;

	//free(ei_vector);
	//free(ei_value);
	free(R);
	free(c_matrix);
	free(a);
	free(b);
	free(c);
}

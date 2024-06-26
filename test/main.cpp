#include "simple_la.h"
#include "qr_method.h"
#include <iomanip>

using namespace std;

int main(){

	int M = 3;
	int N = 4;
	int L = 3;

	std::complex<double> *a;
//	std::complex<double> *b(1,2);
	std::complex<double> *b;
	std::complex<double> *c;
	a = (std::complex<double> *)malloc(sizeof(std::complex<double>)*M*N);
	b = (std::complex<double> *)malloc(sizeof(std::complex<double>)*N*L);
	c = (std::complex<double> *)malloc(sizeof(std::complex<double>)*M*L);
	//b = (std::complex<double> *)malloc(sizeof(std::complex<double>)*4);
	for(int i=0; i<M*N; ++i){
		a[i]= std::complex<double>(random_num(1,2),random_num(1,2));
		//a[i] = std::complex<double>(i%N, i%N);
	}
	printf("Matrix a\n");
	for(int i = 0; i < M; ++i) {
		for(int j = 0; j < N; ++j) {
			cout << a[i*N +j] << " ";
		}
		cout << endl;
	}

	for(int i=0; i<N*L; ++i){
		b[i]= std::complex<double>(random_num(1,2),random_num(1,2));
		//b[i] = std::complex<double>(i%L, i%L);
	}
	printf("Matrix b\n");
	for(int i = 0; i < N; ++i) {
		for(int j = 0; j < L; ++j) {
			cout << b[i*L +j] << " ";
		}
		cout << endl;
	}

	printf("Matrix Multiplication\n");
	matrix_mul(a, b, M, N, N, L, c);
	for(int i = 0; i < M; ++i) {
		for(int j = 0; j < L; ++j) {
			cout << c[i*L +j] << " ";
		}
		cout << endl;
	}


	printf("Matrix cov\n");
	std::complex<double>* c_matrix = (std::complex<double> *)malloc(sizeof(std::complex<double>)*M*M);
	cov_matrix(c, M, L, c_matrix);
	for(int i = 0; i < M; ++i) {
		for(int j = 0; j < M; ++j) {
			cout << c_matrix[i*M +j] << " ";
		}
		cout << endl;
	}
	c_matrix[0] = std::complex<double>(6,0);
	c_matrix[1] = std::complex<double>(5.5,0);
	c_matrix[2] = std::complex<double>(-1,0);
	c_matrix[3] = std::complex<double>(5.5,0);
	c_matrix[4] = std::complex<double>(1,0);
	c_matrix[5] = std::complex<double>(-2,0);
	c_matrix[6] = std::complex<double>(-1,0);
	c_matrix[7] = std::complex<double>(-2,0);
	c_matrix[8] = std::complex<double>(-3,0);
	for(int i = 0; i < M; ++i) {
		for(int j = 0; j < M; ++j) {
			cout << c_matrix[i*M +j] << " ";
		}
		cout << endl;
	}

	//std::complex<double>* r_test = (std::complex<double> *)calloc(M*M,sizeof(std::complex<double>));
	//std::complex<double>* q_test = (std::complex<double> *)calloc(M*M,sizeof(std::complex<double>));
	//householder_reflections(c_matrix, M, q_test, r_test);
	//printf("Q\n");
	//for(int i = 0; i < M; ++i) {
	//	for(int j = 0; j < M; ++j) {
	//		cout << setw(20) << setprecision(2) << q_test[i*M +j] << " ";
	//	}
	//	cout << endl;
	//}
	//printf("R\n");
	//for(int i = 0; i < M; ++i) {
	//	for(int j = 0; j < M; ++j) {
	//		cout << setw(20) << r_test[i*M +j] << " ";
	//	}
	//	cout << endl;
	//}

	//std::complex<double>* c_test = (std::complex<double> *)calloc(M*M,sizeof(std::complex<double>));
	//matrix_mul(r_test, q_test, M, M, M, M, c_test);
	//printf("c\n");
	//for(int i = 0; i < M; ++i) {
	//	for(int j = 0; j < M; ++j) {
	//		cout << setw(20) << setprecision(5) << c_test[i*M +j] << " ";
	//	}
	//	cout << endl;
	//}
	//free(r_test);
	//free(q_test);

	////c_test[3] = std::complex<double>(0,0);
	////c_test[6] = std::complex<double>(0,0);
	////c_test[7] = std::complex<double>(0,0);

	//r_test = (std::complex<double> *)calloc(M*M,sizeof(std::complex<double>));
	//q_test = (std::complex<double> *)calloc(M*M,sizeof(std::complex<double>));
	//householder_reflections(c_test, M, q_test, r_test);
	//printf("Q\n");
	//for(int i = 0; i < M; ++i) {
	//	for(int j = 0; j < M; ++j) {
	//		cout << setw(20) << setprecision(2) << q_test[i*M +j] << " ";
	//	}
	//	cout << endl;
	//}
	//printf("R\n");
	//for(int i = 0; i < M; ++i) {
	//	for(int j = 0; j < M; ++j) {
	//		cout << setw(20) << r_test[i*M +j] << " ";
	//	}
	//	cout << endl;
	//}

	//free(c_test);
	//c_test = (std::complex<double> *)calloc(M*M,sizeof(std::complex<double>));

	//matrix_mul(q_test, r_test, M, M, M, M, c_test);
	//printf("c\n");
	//for(int i = 0; i < M; ++i) {
	//	for(int j = 0; j < M; ++j) {
	//		cout << setw(20) << setprecision(5) << c_test[i*M +j] << " ";
	//	}
	//	cout << endl;
	//}

	//free(r_test);
	//free(q_test);

	std::complex<double>* ei_eve_test = (std::complex<double> *)malloc(M*M*sizeof(std::complex<double>));
	std::complex<double>* ei_val_test = (std::complex<double> *)malloc(M*sizeof(std::complex<double>));
	eigen_vv(c_matrix, M, ei_eve_test, ei_val_test);
	printf("eigen vectors\n");
	for(int i = 0; i < M; ++i) {
		for(int j = 0; j < M; ++j) {
			cout << setw(20) << setprecision(2) << ei_eve_test[i*M +j] << " ";
		}
		cout << endl;
	}

	printf("eigen values\n");
	for(int i=0; i<M; ++i){
		cout << ei_val_test[i] << endl;
	}
	free(ei_eve_test);
	free(ei_val_test);




	free(c_matrix);
	free(a);
	free(b);
	free(c);
}

#include "qr_method.h"

using namespace std;

/********************************************
 ***********QR decomposition method**********
********************************************/

//int qr_decomposition(std::complex<double>* A);

int householder_reflections(std::complex<double>* A, int size, std::complex<double>* &R){
	//std::complex<double>* R = (std::complex<double> *)malloc(sizeof(std::complex<double>) * size * size);
	std::complex<double>* tmp_A = matrix_copy(A, size, size);
	for(int i=0; i<size-1; ++i){
		std::complex<double>* ai = get_col(A, size, size, i);
		std::complex<double>* P = (std::complex<double> *)malloc(sizeof(std::complex<double>) * size * size);
		compute_P(ai, size, i, P);
		R = matrix_mul(P, tmp_A, size, size, size, size);
		tmp_A = matrix_copy(R, size, size);
		free(ai);
		free(P);
	}
	free(tmp_A);
	return 0;
}

int compute_P(std::complex<double>* colA, int size, int col, std::complex<double>* &mat_P){
	// P = I - 2 * n * nT
	// I is identity matrix
	int size_i = size - col;
	mat_P = (std::complex<double>*)identity_mat(size);
	double* mat_i = identity_mat(size_i);
	std::complex<double>* ai = colA + col; // column offset in pointer
	std::complex<double>* vec_n = (std::complex<double> *)malloc(sizeof(std::complex<double>) * size_i);
	compute_n(ai, size_i, vec_n);
	for(int i=0; i<size_i; ++i){
		for(int j=0; j<size_i; ++j){
			int ii = i+col;
			int jj = j+col;
			mat_P[ii*size+jj] = mat_i[i*size_i+j] - (2.0 * vec_n[i*size_i+j] * conj(vec_n[j*size_i+i]));
		}
	}
	return 0;
}

int compute_n(std::complex<double>* colA, int size, std::complex<double>* &vec_n){
	// vector n = vector n/ norm(n)
	std::complex<double>* vec_u = (std::complex<double> *)malloc(sizeof(std::complex<double>) * size);
	double u_norm = vec_norm(vec_u, size); // get the norm of vector u
	compute_u(colA, size, vec_u); // get u vector
	
	for(int i=0; i<size; ++i){
		vec_n[i] = vec_u[i]/u_norm;
	}
	return 0;

}
int compute_u(std::complex<double>* colA, int size, std::complex<double>* &vec_u){
	// vector u = vector Ai - (sign * norm(Ai) * vector b)
	// vector b is the vector for reflection -> b = [1, 0, 0]
	double A_norm = vec_norm(colA, size);
	int sign = -1 * real(colA[0]) / real(colA[0]); // opposite sign to first element of column 
	
	// vector for reflection
	double vec_b[size] = {};
	vec_b[0] = 1;
	
	printf("vector u\n");
	for(int i=0; i<size; ++i){
		vec_u[i] = colA[i] - (sign * A_norm * vec_b[i]);
		cout << vec_u[i] << endl;
	}
	return 0;
}

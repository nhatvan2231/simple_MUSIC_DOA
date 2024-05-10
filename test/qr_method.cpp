#include "qr_method.h"
#include <iomanip>


/********************************************
 ***********QR decomposition method**********
********************************************/

//int qr_decomposition(std::complex<double>* A);

int householder_reflections(std::complex<double>* A, int size, std::complex<double>* &R){
	//std::complex<double>* R = (std::complex<double> *)malloc(sizeof(std::complex<double>) * size * size);
	//std::complex<double>* tmp_A = (std::complex<double> *)malloc(sizeof(std::complex<double>) * size * size); 
	//matrix_copy(A, size, size, tmp_A);
	matrix_copy(A, size, size, R);
	for(int i=0; i<size-1; ++i){
		std::complex<double>* tmp_A = (std::complex<double> *)calloc(size*size, sizeof(std::complex<double>));
		std::complex<double>* subA = (std::complex<double> *)malloc(sizeof(std::complex<double>) * size);
		std::complex<double>* P = (std::complex<double> *)malloc(sizeof(std::complex<double>) * size * size);
		get_col(R, size, size, i, subA);
		identity_mat(size, P); // match the original size of A
		compute_P(subA, size, i, P);
		matrix_mul(P, R, size, size, size, size, tmp_A);
		matrix_copy(tmp_A, size, size, R);
		free(P);
		free(subA);
		free(tmp_A);
	}
	return 0;
}

int compute_P(std::complex<double>* colA, int size, int col, std::complex<double>* &mat_P){
	// P = I - 2 * n * nT
	// I is identity matrix
	int size_i = size - col;
	std::complex<double>* mat_i = (std::complex<double> *)malloc(sizeof(std::complex<double>) * size_i * size_i);
	std::complex<double>* vec_n = (std::complex<double> *)malloc(sizeof(std::complex<double>) * size_i);
	std::complex<double>* vec_ncov = (std::complex<double> *)malloc(sizeof(std::complex<double>) * size_i * size_i);

	identity_mat(size_i, mat_i); // size of sub A
	std::complex<double>* ai = colA + col; // column offset in pointer
	//printf("col A\n");
	//for(int i=0; i<size_i; ++i){
	//	std::cout << ai[i] << std::endl;
	//}
	compute_n(ai, size_i, vec_n);
	cov_matrix(vec_n, size_i ,1 , vec_ncov);
	for(int i=0; i<size_i; ++i){
		for(int j=0; j<size_i; ++j){
			int ii = i+col;
			int jj = j+col;
			//mat_P[ii*size+jj] = mat_i[i*size_i+j] - (2.0 * vec_n[i*size_i+j] * conj(vec_n[j*size_i+i]));
			mat_P[ii*size+jj] = mat_i[i*size_i+j] - (2.0 * vec_ncov[i*size_i+j]);
		}
		//std::cout << std::endl;
	}
	//printf("Compute P %d\n", col);
	//for(int i = 0; i < size; ++i) {
	//	for(int j = 0; j < size; ++j) {
	//		std::cout << std::setw(25) << std::setprecision(5) << mat_P[i*size +j] << " ";
	//	}
	//	std::cout << std::endl;
	//}
	free(mat_i);
	free(vec_ncov);
	//free(ai);
	free(vec_n);
	//printf("DONE P\n");
	return 0;
}

int compute_n(std::complex<double>* colA, int size, std::complex<double>* &vec_n){
	// vector n = vector u/ norm(u)
	std::complex<double>* vec_u = (std::complex<double> *)malloc(sizeof(std::complex<double>) * size);
	compute_u(colA, size, vec_u); // get u vector
	double u_norm = vec_norm(vec_u, size); // get the norm of vector u
	
	//printf("Compute n\n");
	for(int i=0; i<size; ++i){
		vec_n[i] = vec_u[i]/u_norm;
		//std::cout << vec_n[i] << std::endl;
	}
	free(vec_u);
	//printf("DONE n\n");
	return 0;

}
int compute_u(std::complex<double>* colA, int size, std::complex<double>* &vec_u){
	// vector u = vector Ai - (sign * norm(Ai) * vector b)
	// vector b is the vector for reflection -> b = [1, 0, 0]
	double A_norm = vec_norm(colA, size);
	std::complex<double> sign;
	if (real(colA[0]) > 0) sign = std::complex<double>(-1,0);
	else sign = std::complex<double>(1,0); // opposite sign to first element of column 
	
	// vector for reflection
	double* vec_b = (double*)malloc(sizeof(double)*size);
	//double vec_b[size] = {};
	vec_b[0] = 1;
	
	for(int i=0; i<size; ++i){
		std::complex<double> tmp = sign *   A_norm * vec_b[i];
		vec_u[i] = colA[i] - tmp;
		//std::cout << vec_u[i] << " " << colA[i] <<  std::endl;
	}
	//printf("DONE u\n");
	free(vec_b);
	return 0;
}

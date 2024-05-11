#include "qr_method.h"
#include <iomanip>


/********************************************
 ***********QR decomposition method**********
********************************************/
int isRowEchelon(std::complex<double>* A, int size, double eps){
	double sum=0;
	for(int i=0; i<size;++i){
		for(int j=0; j<i; ++j){
			sum += std::abs(A[i*size+j]);
		}
	}
	return eps > sum? 1:0;
}

int eigen_vv(std::complex<double>* A, int size, std::complex<double>* &ei_vectors, std::complex<double>* &ei_values){
	//copy the orignal matrix
	std::complex<double>* mat_A = (std::complex<double> *)malloc(size * size * sizeof(std::complex<double>));
	matrix_copy(A, size, size, mat_A);

	identity_mat(size, ei_vectors);

	int iteration = 5;
	int i_counter = 0;
	bool flagstop = false;
	while((i_counter < iteration) && !flagstop){
		std::complex<double>* Q = (std::complex<double> *)calloc(size*size,sizeof(std::complex<double>));
		std::complex<double>* R = (std::complex<double> *)calloc(size*size,sizeof(std::complex<double>));
		std::complex<double>* tmp_A = (std::complex<double> *)calloc(size * size, sizeof(std::complex<double>));
		std::complex<double>* tmp_eve = (std::complex<double> *)calloc(size * size, sizeof(std::complex<double>));

		printf("fuk\n");
		for(int i = 0; i < size; ++i) {
			for(int j = 0; j < size; ++j) {
				std::cout << std::setw(25) << mat_A[i*size +j] << " ";
			}
			std::cout << std::endl;
		}

		// compute QR
		householder_reflections(mat_A, size, Q, R);

		//compute A
		matrix_mul(Q, R, size, size, size, size, tmp_A);

		printf("A tmp\n");
		for(int i = 0; i < size; ++i) {
			for(int j = 0; j < size; ++j) {
				std::cout << std::setw(25) << tmp_A[i*size +j] << " ";
			}
			std::cout << std::endl;
		}
		matrix_copy(tmp_A, size, size, mat_A);

		//compute ei_vectors = Q0*Q1*Q2*Q3...
		//matrix_mul(ei_vectors, Q, size, size, size, size, tmp_eve);
		//matrix_copy(tmp_eve, size, size,0, ei_vectors);

		printf("Q temp\n");
		for(int i = 0; i < size; ++i) {
			for(int j = 0; j < size; ++j) {
				std::cout << std::setw(25) << Q[i*size +j] << " ";
			}
			std::cout << std::endl;
		}
		printf("R temp\n");
		for(int i = 0; i < size; ++i) {
			for(int j = 0; j < size; ++j) {
				std::cout << std::setw(25) << R[i*size +j] << " ";
			}
			std::cout << std::endl;
		}

		printf("A mat\n");
		for(int i = 0; i < size; ++i) {
			for(int j = 0; j < size; ++j) {
				std::cout << std::setw(25) << mat_A[i*size +j] << " ";
			}
			std::cout << std::endl;
		}

		// check if close enough
		//flagstop = isRowEchelon(mat_A, size);
		
		i_counter++;
		free(Q);
		free(R);
		free(tmp_A);
		free(tmp_eve);
	}
	// the eigen values is the diagnol matrix
	//for(int i=0; i<size; ++i) ei_values[i] = mat_A[i*size+i];
	free(mat_A);
	printf("DONE: %d\n", i_counter);
	return 0;
}

// QR decomposition using householder reflections algorithm
int householder_reflections(std::complex<double>* A, int size, std::complex<double>* &Q, std::complex<double>* &R){
	identity_mat(size, Q);
	matrix_copy(A, size, size, R);
	
	for(int i=0; i<size; ++i){
		std::complex<double>* tmp_Q = (std::complex<double> *)calloc(size*size, sizeof(std::complex<double>));
		std::complex<double>* tmp_R = (std::complex<double> *)calloc(size*size, sizeof(std::complex<double>));
		std::complex<double>* tmp_H = (std::complex<double> *)malloc(sizeof(std::complex<double>) * size * size);
		std::complex<double>* tmp_HT = (std::complex<double> *)malloc(sizeof(std::complex<double>) * size * size);
		std::complex<double>* subA = (std::complex<double> *)malloc(sizeof(std::complex<double>) * size);

		// Compue the H
		get_col(R, size, size, i, subA);
		identity_mat(size, tmp_H); // match the original size of A
		compute_H(subA, size, i, tmp_H);

		//printf("H temp\n");
		//for(int i = 0; i < size; ++i) {
		//	for(int j = 0; j < size; ++j) {
		//		std::cout << std::setw(20) << tmp_H[i*size +j] << " ";
		//	}
		//	std::cout << std::endl;
		//}

		// Compute Q 
		conj_transpose(tmp_H, size, size, tmp_HT);
		matrix_mul(Q, tmp_H, size, size, size, size, tmp_Q);
		matrix_copy(tmp_Q, size, size, Q);

		//printf("Q temp\n");
		//for(int i = 0; i < size; ++i) {
		//	for(int j = 0; j < size; ++j) {
		//		std::cout << std::setw(20) << Q[i*size +j] << " ";
		//	}
		//	std::cout << std::endl;
		//}

		// Compute R
		matrix_mul(tmp_H, R, size, size, size, size, tmp_R);
		matrix_copy(tmp_R, size, size, R);
		//printf("R temp\n");
		//for(int i = 0; i < size; ++i) {
		//	for(int j = 0; j < size; ++j) {
		//		std::cout << std::setw(20) << R[i*size +j] << " ";
		//	}
		//	std::cout << std::endl;
		//}

		free(subA);
		free(tmp_HT);
		free(tmp_H);
		free(tmp_R);
		free(tmp_Q);
	}
	return 0;
}

int compute_H(std::complex<double>* colA, int size, int col, std::complex<double>* &mat_H){
	// H = I - 2 * n * nT
	// I is identity matrix
	int size_i = size - col;
	std::complex<double>* mat_i = (std::complex<double> *)malloc(sizeof(std::complex<double>) * size_i * size_i);
	std::complex<double>* vec_n = (std::complex<double> *)malloc(sizeof(std::complex<double>) * size_i);
	std::complex<double>* vec_ncov = (std::complex<double> *)malloc(sizeof(std::complex<double>) * size_i * size_i);

	identity_mat(size_i, mat_i); // size of sub A
	std::complex<double>* ai = colA + col; // column offset in pointer
	compute_n(ai, size_i, vec_n);
	cov_matrix(vec_n, size_i ,1 , vec_ncov);
	for(int i=0; i<size_i; ++i){
		for(int j=0; j<size_i; ++j){
			int ii = i+col;
			int jj = j+col;
			mat_H[ii*size+jj] = mat_i[i*size_i+j] - (2.0 * vec_ncov[i*size_i+j]);
		}
	}
	free(mat_i);
	free(vec_ncov);
	free(vec_n);
	return 0;
}

int compute_n(std::complex<double>* colA, int size, std::complex<double>* &vec_n){
	// vector n = vector u/ norm(u)
	//printf("subA1\n");
	//for(int i=0; i<size;++i) std::cout << colA[i] << std::endl;
	std::complex<double>* vec_u = (std::complex<double> *)malloc(sizeof(std::complex<double>) * size);
	compute_u(colA, size, vec_u); // get u vector
	double u_norm = vec_norm(vec_u, size); // get the norm of vector u
	
	for(int i=0; i<size; ++i){
		vec_n[i] = vec_u[i]/u_norm;
	}
	free(vec_u);
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
	vec_b[0] = 1;
	
	for(int i=0; i<size; ++i){
		std::complex<double> tmp = sign *   A_norm * vec_b[i];
		vec_u[i] = colA[i] - tmp;
	}
	free(vec_b);
	return 0;
}

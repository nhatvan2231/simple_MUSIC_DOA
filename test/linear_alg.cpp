#include "linear_alg.h"

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

// copy matrix/vector
std::complex<double>* matrix_copy(std::complex<double> *src, int M, int N){
	std::complex<double> *result = (std::complex<double> *)malloc(sizeof(std::complex<double>) * M * N);
	for(int i=0; i<M*N; ++i){
		result[i] = src[i];
	}
	return result;
}

// calculate matrix multiplication
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

// calculate the conjugate transpose of a matrix or vector
std::complex<double>* conj_transpose(std::complex<double> *src, int M, int N){
	std::complex<double> *result = (std::complex<double> *)malloc(sizeof(std::complex<double>) * N * M);
	for(int i=0; i<N; ++i){
		for(int j=0; j<M; ++j){
			result[i*M+j] = conj(src[j*N+i]);
		}
	}
	return result;
}

// calculate the covariance matrix without normalize
std::complex<double>* cov_matrix(std::complex<double> *src, int M, int N){
	std::complex<double> *srcT = conj_transpose(src, N, M);
	std::complex<double> *result = matrix_mul(src, srcT, M, N, N, M);
	free(srcT);
	return result;
}

// get the i column of the matrix
std::complex<double>* get_col(std::complex<double>* src, int M, int N, int col_i){
	std::complex<double>* col = (std::complex<double> *)malloc(sizeof(std::complex<double>) * M);
	for(int i=0; i<M; ++i){
		col[i] = src[i*N + col_i];
	}
	return col;
}

// get the identity matrix with size x size
double* identity_mat(int size){
	double* mat_I = (double*)malloc(sizeof(double) * size * size);
	for(int i=0; i<size; ++i){
		for(int j=0; j<size; ++j){
			mat_I[i*size +j] = 0;
			if(i==j) mat_I[i*size +j] = 1;
		}
	}
	return mat_I;
}

// calculate the norm of vector
double vec_norm(std::complex<double> *src, int size){
	double result;
	for(int i=0; i<size; ++i){
		result += norm(src[i]);
	}
	return sqrt(result);
}

// generate random number
double random_num(double min, double max){
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(min, max);
	return dis(gen);
}

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



















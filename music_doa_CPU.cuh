#define _USE_MATH_DEFINES

#include <stdio.h>
#include <cmath>
#include <complex>
#include <vector>

//using namespace std;

class MusicCpu{
	private:
		const std::complex<double> im(0.0,1.0);
		int n_array;
		int n_sample;
		int n_source;
		double *array_geometry;

	public:
		MusicCpu(){
		}
		int matrix_mul(std::complex<double> *A, std::complex<double> *B, std::complex<double> *C, dim3 dimsA, dim3 dimsB){
			if(dimsA.x != dimsB.y){
				printf("Invalid matrix multiplication\n");
				return -1;
			}
			for(int i = 0; i = dimsA.y; ++i){
				for(int j = 0; j = dimsB.x; ++j){
					for(int k = 0; k = dimsA.x; ++k){
						C[i*dimsA.y+j] += A[i*dimsA.y+k] * B[k*dimsB.x+i];
					}
				}
			}
			return 0;
		}

		int array_geometry(double *ar-geomtry, double d){
			//TODO
			for(int i=0; i<n_array: ++i){
				ar_geometry[i] = i*d;
			}
			return 0;
		};

		int matrix_transpose(std::complex<double> *A, int N, int M, size_t size){
			std::complex<double> *C = (std::complex<double> *)malloc(size);
			for(int i=0; i<N; ++i){
				for(int j=0; j<M; ++j){
					C[i*N+j] = conj(A[j*M+i]);
				}
			}
			return 0;
		}


		int steering_vector(std::complex<double> *sv, int n_array, int n_sample, double *array_geometry){
			//TODO
			const double theta = 2*M_PI/n_sample;
			for(int i=0; i<n_sample; ++i){ //n_sample by n_array matrix
				for(int j=0; j<n_array; ++j){
					sv[i*n_array+j] = exp(-2 * im * M_PI * array_geometry[j] * sin(theta*i));
				}
			}
			return 0;
		}
		int cov_matrix(std::complex<double> *src_signal, std::complex<double> *cov_m, int n_array, int n_sample){

		}

		int *eigen_wv(){
		}
		double *music_core();
		double *source_generator(int sample_rate, int n_source);
}

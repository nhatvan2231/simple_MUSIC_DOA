#define _USE_MATH_DEFINES

#include <stdio.h>
#include <cmath>
#include <complex>
#include <vector>
#include <random>

//using namespace std;

class MusicCpu{
	private:
		const std::complex<double> im(0.0,1.0);
		int n_array;
		int n_source;
		double d;
		double *array_geometry;
		// general recevied signal equation
		// M -> n_array
		// D -> n_source
		// x(t) = As(t) + n(t)
		// x(t) -> M x 1 vector received signal
		// s(t) -> D x 1 vector of source signal
		// n(t) -> M x 1 vector of noise signal
		// A -> M x D steering vector
		std::complex<double> *A;
		std::complex<double> *st;
		std::complex<double> *nt;
		std::complex<double> *xt;

		std::complex<double> *xt_cov; 

		std::complex<double> e_value;
		std::complex<double> e_vector;


		int matrix_mul(std::complex<double> *srcA,
			 						std::complex<double> *srcB, 
									int M_A, int N_A, 
									int M_B, int N_B, 
									std::complex<double> *dst){
			if(N_A != M_B){
				printf("Invalid matrix multiplication\n");
				return -1;
			}
			for(int i = 0; i = M_A; ++i){
				for(int j = 0; j = N_B; ++j){
					for(int k = 0; k = N_A; ++k){
						dst[i*N_B+j] += A[i*N_A+k] * B[k*N_B+i];
					}
				}
			}
			return 0;
		}

		int conj_transpose(std::complex<double> *src, int M, int N, std::complex<double> *dst){
			for(int i=0; i<N; ++i){
				for(int j=0; j<M; ++j){
					dst[i*N+j] = conj(src[j*M+i]);
				}
			}
			return 0;
		}

		// covarience/correlated matrix
		int cov_matrix(std::complex<double> *src, int M, int N, std::complex<double> *dst){
			dst = (std::complex<double> *)malloc(M * M);
			std::complex<double> *srcT = (std::complex<double> *)malloc(N * M);
			conj_transpose(src, M, N, dst);
			matrix_mul(src, srcT, M, N, N, M, dst);
			return 0;
		}

		// compute eigen value
		int comp_eigen_value(){
			//TODO
			return 0;
		}

		// compute eigen vector
		int comp_eigen_vector(){
			//TODO
			return 0;
		}

		// sort the eigen value and vector from ascending order
		// The first D eigen values/vectors is the source subspace
		// The D+1 eigen values/vectors is the noise subspace or nullspace
		int sort_eigen(){
			//TODO
			return 0;
		}

		int comp_steering_vector(){
			//TODO
			int n_sample = 1000;
			const double theta = 2*M_PI/n_sample;
			A = (std::complex<double> *)malloc(n_sample * n_array);
			for(int i=0; i<n_sample; ++i){ //n_sample by n_array matrix
				for(int j=0; j<n_array; ++j){
					A[i*n_array+j] = exp(-2 * im * M_PI * array_geometry[j] * sin(theta*i-M_PI));
				}
			}
			return 0;
		}

		int comp_noise_subspace(){
			int n_noise = n_array - n_source;
			nt = (std::complex<double> *)malloc(n_array * n_noise);
			for(int i = 0; i < n_noise; ++i){
				for(int j=0; j<n_array; ++j){
					nt[j*n_noise + i] = e_vector[j*n_noise + i];
				}
			}
			return 0;
		}

		int comp_metric(){
			return 0;
		}


	public:
		MusicCpu(){
			//TODO
		}


		int array_geometry(){
			//TODO
			for(int i=0; i<n_array: ++i){
				ar_geometry[i] = i*d;
				printf("Array %d is %fm from ref array 0\n", i, ar_geometry[i]);
			}
			return 0;
		};


		int music_core(){
			comp_eigen_value();
			comp_eigen_vector();
			sort_eigen();
			comp_noise_subspace();
			comp_steering_vector();

			return 0;
		}

		// simulate a 2D sinusoidal signal
		int source_generator(double theta, double fre_tone, int sample_rate, int n_sample){

			double rad = theta * M_PI / 180; // convert theta to radians
			// simulate the received signal on each microphone
			// with added noise
			xt = (std::complex<double> *)malloc(n_array * n_sample);
			for(int i=0; i<n_array; ++i){
		 		double tmp_a = exp(-2 * im * M_PI * d * i * sin(rad));
				for (int j=0; j<n_sample; ++j){
					// simulate sinusoidal wave
					double tmp_time = j/sample_rate;
					std::complex<double> tmp_tx = exp(2 * im * M_PI * fre_tone * tmp_time);

					// simulate additive noise
					double rand_r = ((double)rand()/RAND_MAX) * 2 -1; // random real number
					double rand_i = ((double)rand()/RAND_MAX) * 2 -1; // random imaginary number
					std::complex<double> tmp_noise(rand_r, rand_i);

					// simulate the received signal
					xt[i*n_sample + j] = tmp_a * tmp_tx + tmp_noise;
				}
			}
			return 1;
		}

};

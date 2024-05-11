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
		int n_noise;
		int n_scan;
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

		std::complex<double> *e_value;
		std::complex<double> *e_vector;
		
		std::complex<double> doa;



		int comp_steering_vector(){
			n_scan = 1000;
			const double theta = 2*M_PI/n_scan;
			A = (std::complex<double> *)malloc(sizeof(std::complex<double> * n_scan * n_array);
			for(int i=0; i<n_scan; ++i){ //n_sample by n_array matrix
				for(int j=0; j<n_array; ++j){
					A[i*n_array+j] = exp(-2 * im * M_PI * array_geometry[j] * sin(theta*i-M_PI));
				}
			}
			return 0;
		}

		int comp_noise_subspace(){
			n_noise = n_array - n_source;
			nt = (std::complex<double> *)malloc(sizeof(std::complex<double> * n_array * n_noise);
			for(int i = 0; i < n_noise; ++i){
				for(int j=0; j<n_array; ++j){
					nt[j*n_noise + i] = e_vector[j*n_noise + i];
				}
			}
			return 0;
		}

		int comp_metric(){
			std::complex<double> metric;
			std::complex<double> *result_matrix;
			std::complex<double> *nt_T;
			std::complex<double> *A_T;
			std::complex<double> *tmp_matrix;

			conj_transpose(nt, n_noise, n_array, nt_T);
			conj_transpose(A, n_scan, n_array, A_T);


			matrix_mul(nt_T, A, n_noise, n_array, n_array, n_scan, result_matrix);

			return 0;
		}

		int array_geometry(){
			//TODO
			for(int i=0; i<n_array: ++i){
				ar_geometry[i] = i*d;
				printf("Array %d is %fm from ref array 0\n", i, ar_geometry[i]);
			}
			return 0;
		};


	public:
		MusicCpu(int num_array, int num_source, double distant){
			//TODO
			n_array = num_array;
			n_source = num_source;
			d = distant;
			aray_geometry();
		}

		~MusicCpu(){
			free(array_geometry);
			free(A);
			free(st);
			free(nt);
			free(xt);
			free(xt_cov); 
			free(e_value);
			free(e_vector);
		}




		int music_core(){
			comp_eigen_value();
			comp_eigen_vector();
			sort_eigen();
			comp_noise_subspace();
			comp_steering_vector();
			return 0;
		}

		int set_signal(std::complex<double> *signal){
			//TODO
			return 0;
		}

		// simulate a 2D sinusoidal signal
		int signal_generator(double theta, double fre_tone, int sample_rate, int n_sample){

			double rad = theta * M_PI / 180; // convert theta to radians
			// simulate the received signal on each microphone
			// with added noise
			xt = (std::complex<double> *)malloc(sizeof(std::complex<double> * n_array * n_sample);
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

//#include "simple_la.cuh"
#include <thrust/complex.h>
#include <stdio.h>
#include <iostream>

using namespace std;

inline
cudaError_t checkCuda(cudaError_t result)
{
#if defined(DEBUG) || defined(_DEBUG)
	if (result != cudaSuccess) {
		fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
		assert(result == cudaSuccess);
	}
#endif
	return result;
}


// copy matrix/vector 
__global__ void matrixCopy(thrust::complex<double> *src, thrust::complex<double>* dst, int N){
        int idx = blockIdx.x * blockDim.x + threadIdx.x;
        int idy = blockIdx.y * blockDim.y + threadIdx.y;
        dst[idx*N + idy] = src[idx*N + idy];
}

// matrix multiplcation
__global__ void matrixMul(thrust::complex<double>* srcA, 
			thrust::complex<double>* srcB, 
			thrust::complex<double>* dstC,
			int N){
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int idy = threadIdx.y + blockDim.y * blockIdx.y;

	thrust::complex<double>  c_temp = 0.0f;
	if (idx < N && idy < N){
		for(int i = 0; i < N; ++i){
			c_temp += srcA[idy*N + i] * srcB[i*N + idx];// + 0.0f;
		}
		dstC[idy * N + idx] += c_temp;
	}
	__syncthreads();

}

// conjugate tranpose
__global__ void conj_transpose(thrust::complex<double>* src, thrust::complex<double>* dst){
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int idy = threadIdx.y + blockDim.y * blockIdx.y;

	dst[idy*blockDim.y+idx] = thrust::conj(src[idx*blockDim.x+idy]);
}
	
void cov_matrix(thrust::complex<double>* src, int M, int N, thrust::complex<double>* &dst){
	size_t size = (M*N)*sizeof(thrust::complex<double>);
	thrust::complex<double>* srcT = (thrust::complex<double> *)malloc(size);
	if(M*N <= 32*32)
		dim3 grid_block(1,1);
	else if(M*N/2 <= 32*32)
		dim3 grid_block(2,2);
		
	dim3 block_thread(N/grid_block.x, N/grid_block.y);


	checkCuda(cudaFreeHost(srcT));
}

int main(void){
	int N = 2;
	size_t size = (N*N)*sizeof(thrust::complex<double>);

	thrust::complex<double>* h_A = (thrust::complex<double> *)malloc(size);
	thrust::complex<double>* h_C = (thrust::complex<double> *)malloc(size);
	for(int i = 0; i < N*N; ++i){
	//	h_A[i] = thrust::complex<double>(rand()/(float)RAND_MAX, rand()/(float)RAND_MAX);
		h_A[i] = thrust::complex<double>(i, 0);
	}

	thrust::complex<double>* d_A = NULL;
	thrust::complex<double>* d_C = NULL;

	checkCuda(cudaMalloc((void**)&d_A, size));
	checkCuda(cudaMalloc((void**)&d_C, size));

	checkCuda(cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice));

	dim3 grid_block(1,1);
	dim3 block_thread(N/grid_block.x, N/grid_block.y);
	matrixCopy<<<grid_block, block_thread>>>(d_A, d_C, N);

	checkCuda(cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost));

	for(int i = 0; i < N*N; ++i){
		if(h_C[i] != h_A[i]) 
			printf("FUCK\n");
	}
	printf("DONE\n");

	matrixMul<<<grid_block, block_thread>>>(d_A,d_A, d_C, N);
	checkCuda(cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost));

	for(int i=0;i<N;++i){
		for(int j=0;j<N;++j){
			cout << h_C[i*N+j] << " ";
		}
		cout << endl;
	}

	conj_transpose<<<grid_block, block_thread>>>(d_A, d_C);
	checkCuda(cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost));
	for(int i=0;i<N;++i){
		for(int j=0;j<N;++j){
			cout << h_C[i*N+j] << " ";
		}
		cout << endl;
	}
			

	checkCuda(cudaFreeHost(h_A));
	checkCuda(cudaFreeHost(h_C));
	checkCuda(cudaFree(d_A));
	checkCuda(cudaFree(d_C));
}

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include "EXECUTE_HELPER.h"

__device__ int getGlobalIdx_1D_2D(){
    return blockIdx.x * blockDim.x * blockDim.y + threadIdx.y * blockDim.x + threadIdx.x;
}

__global__ void cuda_Hello_World(){
    int tid = getGlobalIdx_1D_2D();
    printf("Hello World from thread ID:%d\n",tid);
}
//op_nop - no operator
//op_mat_inv -  matrix inversion
__device__ void op_mat_inv(){

}
//op_mat_mul - matrix multiplication
__device__ void op_mat_mul(){

}
//op mat add - matrix addition
__device__ void op_mat_add(){

}
//op_vec_mul - multiply two vectors
__device__ void op_vec_mul(){

}
//op_vec_add - add two vectors
__global__ void op_vec_add(double * result,double * vector0,double * vector1,double scalar,int size){
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid < size){
        result[tid] = vector0[tid] + scalar * vector1[tid];
    }
}

void cuda_kernel_vector_addition(double * result,double * vector0,double * vector1, double scalar,int size){
    
    double * device_vector0;
    double * device_vector1;
    double * device_result;

    cudaMalloc((void**)&device_vector0,size * sizeof(double) );
    cudaMalloc((void**)&device_vector1,size * sizeof(double) );
    cudaMalloc((void**)&device_result,size * sizeof(double));
    //copy the input vectors into the device
    cudaMemcpy(device_vector0,vector0,size * sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(device_vector1,vector1,size * sizeof(double),cudaMemcpyHostToDevice);

    op_vec_add<<<1,size>>>(device_result,device_vector0,device_vector1,scalar,size);

    cudaDeviceSynchronize();

    cudaMemcpy(result,device_result,size * sizeof(double),cudaMemcpyDeviceToHost);

    cudaFree(device_vector0);
    cudaFree(device_vector0);
    cudaFree(device_result);
}



void Execute_Helper_Kernel(int){
    //create cublas Handle here
    cublasStatus_t status;
    cublasHandle_t handle;
    status = cublasCreate(&handle);
    
    if (status != CUBLAS_STATUS_SUCCESS) {
        fprintf(stderr, "!!!! CUBLAS initialization error\n");
        return;
    }else{
        printf("CUBLAS initialized succesfully");
    }
    //create and call the schedhuler here

    //return the result back into the parameter
    cudaDeviceSynchronize();
} 

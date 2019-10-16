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
__device__ void op_vec_add(){

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

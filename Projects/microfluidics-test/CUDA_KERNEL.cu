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

void cuda_dgemm_helper(bool at,bool bt,int m, int n, int k,const double * alpha,const double * A,int lda,const double * B,int ldb,const double * beta,double * C,double ldc){
    cublasStatus_t status;
    cublasHandle_t handle;
    status = cublasCreate(&handle);
    status  = cublasDgemm(handle,at?CUBLAS_OP_T:CUBLAS_OP_N,bt?CUBLAS_OP_T:CUBLAS_OP_N,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
    
    if (status != CUBLAS_STATUS_SUCCESS) {
        printf("TIMES MM not working");
    }

    cublasDestroy(handle);

    return;
}

void cuda_dgemm(bool at,bool bt,int m, int n, int k,const double * alpha,const double * A,int lda,const double * B,int ldb,const double * beta,double * C,double ldc){
    //allocate device variables here
    double * device_A;
    double * device_B;
    double * device_C;

    cudaMalloc(&device_A,m*k*sizeof(double));
    cudaMalloc(&device_B,n*k*sizeof(double));
    cudaMalloc(&device_C,m*n*sizeof(double));

    cudaMemcpy(device_A,A,m * k * sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(device_B,B,n* k * sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(device_C,C,m * n * sizeof(double),cudaMemcpyHostToDevice);

    cuda_dgemm_helper(at,bt,m,n,k,alpha,device_A,lda,device_B,ldb,beta,device_C,ldc);
    
    cudaMemcpy(C,device_C,m * n * sizeof(double),cudaMemcpyDeviceToHost);
    
    cudaDeviceSynchronize();

    cudaFree(device_A);
    cudaFree(device_B);
    cudaFree(device_C);

    return;
}

void cuda_gemv_helper(bool t,int m, int n,const double * alpha,const double * A, int lda, const double * x, int cx, const double * beta,double * y,int incy)
{
    cublasStatus_t status;
    cublasHandle_t handle;
    
    status = cublasCreate(&handle);
    status = cublasDgemv(handle, t?CUBLAS_OP_T:CUBLAS_OP_N,m,n,alpha,A,lda,x,cx,beta,y,incy);
    
    if (status != CUBLAS_STATUS_SUCCESS) {
        printf("TIMES MM not working");
    }

    cublasDestroy(handle);

    return;
}

void cuda_dgemv(bool t,int m, int n,const double * alpha,const double * A, int lda, const double * x, int cx, const double * beta,double * y,int incy)
{
    double * device_x; //vector
    double * device_A; //matrix
    double * device_y; //vector 

    cudaMalloc(&device_x,(t?m:n)*sizeof(double));
    cudaMalloc(&device_A,m*n*sizeof(double));
    cudaMalloc(&device_y,(t?n:m)*sizeof(double));

    cudaMemcpy(device_x,x,(t?m:n)*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(device_A,A,m* n * sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(device_y,y,(t?n:m)*sizeof(double),cudaMemcpyHostToDevice);

    cuda_gemv_helper(t,m,n,alpha,device_A,lda,device_x,cx,beta,device_y,incy);

    cudaMemcpy(y,device_y,(t?n:m)*sizeof(double),cudaMemcpyDeviceToHost);

    cudaDeviceSynchronize();

    cudaFree(device_x);
    cudaFree(device_A);
    cudaFree(device_y);

    return;
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

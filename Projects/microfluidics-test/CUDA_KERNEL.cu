#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

__device__ int getGlobalIdx_1D_2D(){
    return blockIdx.x * blockDim.x * blockDim.y + threadIdx.y * blockDim.x + threadIdx.x;
}

__global__ void cuda_Hello_World(){
    int tid = getGlobalIdx_1D_2D();
    printf("Hello World from thread ID:%d\n",tid);
}

void Hello_World(int){
    dim3 block(3,3);
    dim3 grid(1);
    cuda_Hello_World<<<grid,block>>>();
    cudaDeviceSynchronize();
} 

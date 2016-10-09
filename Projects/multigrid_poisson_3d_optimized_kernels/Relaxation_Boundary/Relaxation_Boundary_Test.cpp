//#####################################################################
// Copyright 2009, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#ifdef ENABLE_TIMING
#include <x86intrin.h>
#endif

#include "../Thread_Queueing/PTHREAD_QUEUE.h"
#include "Relaxation_Boundary_Helper.h"
using namespace PhysBAM;
extern PhysBAM::PTHREAD_QUEUE* pthread_queue;

struct timeval starttime,stoptime;
void start_timer(){gettimeofday(&starttime,NULL);__rdtsc();}
void stop_timer(){gettimeofday(&stoptime,NULL);__rdtsc();}
double get_time(){return (double)stoptime.tv_sec-(double)starttime.tv_sec+(double)1e-6*(double)stoptime.tv_usec-(double)1e-6*(double)starttime.tv_usec;}

int main(int argc,char *argv[]){
    
    typedef float T;
    static const int size=512;
    T* const u=new T[size];
    T* const b=new T[size];
    T* const one_over_diagonal_part=new T[size];
    int* const boundary_index=new int[size];
    int* const block_start=new int [size];
    int* const block_end=new int[size];
    int number_of_red_blocks=size/2;
    int number_of_black_blocks=size/2;
    int loops=1;
    bool reverse_order=false;


    if(argc!=2){printf("Must specify number of threads\n");exit(1);}
    int number_of_threads=atoi(argv[1]);
    printf("Using %d threads\n",number_of_threads);
    pthread_queue=new PhysBAM::PTHREAD_QUEUE(number_of_threads);  

    // Allocate data
    printf("Allocation");start_timer();

//    Relaxation_Boundary_Size_Specific_Helper<T,size>::Allocate_Data(u,b,r,bit_writemask);

    stop_timer();printf(" [Seconds: %g]\n",get_time());

    // Initialize data
    printf("Initialization");start_timer();

//    Relaxation_Boundary_Size_Specific_Helper<T,size>::Initialize_Data(u,b,r,bit_writemask);
    Relaxation_Boundary_Helper<T> test(size,size,size,u,b,one_over_diagonal_part,boundary_index,block_start,block_end,number_of_red_blocks,number_of_black_blocks,loops,reverse_order);

    stop_timer();printf(" [Seconds: %g]\n",get_time());

    while(1){
        // Running benchmark
        printf("Running benchmark");start_timer();

        test.Relax_Parallel(number_of_threads);
        
        stop_timer();printf(" [Seconds: %g]\n",get_time());
    }

    return 0;
}

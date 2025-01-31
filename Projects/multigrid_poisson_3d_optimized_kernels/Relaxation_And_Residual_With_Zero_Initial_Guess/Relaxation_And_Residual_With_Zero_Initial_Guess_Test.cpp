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
#include "Relaxation_And_Residual_With_Zero_Initial_Guess_Helper.h"
using namespace PhysBAM;
extern PhysBAM::PTHREAD_QUEUE* pthread_queue;

struct timeval starttime,stoptime;
void start_timer(){gettimeofday(&starttime,NULL);__rdtsc();}
void stop_timer(){gettimeofday(&stoptime,NULL);__rdtsc();}
double get_time(){return (double)stoptime.tv_sec-(double)starttime.tv_sec+(double)1e-6*(double)stoptime.tv_usec-(double)1e-6*(double)starttime.tv_usec;}

int main(int argc,char *argv[]){
    
    typedef float T;
    static const int size=512;

    T* u;
    T* b;
    T* r;
    unsigned char* bit_writemask;
    unsigned char* bit_interiormask;

    if(argc!=2){printf("Must specify number of threads\n");exit(1);}
    int number_of_threads=atoi(argv[1]);
    printf("Using %d threads\n",number_of_threads);
    pthread_queue=new PhysBAM::PTHREAD_QUEUE(number_of_threads);  

    // Allocate data
    printf("Allocation");start_timer();

    Relaxation_And_Residual_With_Zero_Initial_Guess_Size_Specific_Helper<T,size>::Allocate_Data(u,b,r,bit_writemask,bit_interiormask);

    stop_timer();printf(" [Seconds: %g]\n",get_time());

    // Initialize data
    printf("Initialization");start_timer();

    Relaxation_And_Residual_With_Zero_Initial_Guess_Size_Specific_Helper<T,size>::Initialize_Data(u,b,r,bit_writemask,bit_interiormask);
    Relaxation_And_Residual_With_Zero_Initial_Guess_Helper<T> test(size,size,size,u,b,r,bit_writemask,bit_interiormask,1);

    stop_timer();printf(" [Seconds: %g]\n",get_time());

    while(1){
        // Running benchmark
        printf("Running benchmark");start_timer();

        test.Run_Parallel(number_of_threads);
        
        stop_timer();printf(" [Seconds: %g]\n",get_time());
    }

    return 0;
}

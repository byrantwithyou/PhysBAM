//#####################################################################
// Copyright 2009, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "Multiplication_And_Dot_Product_Helper.h"
#include "PTHREAD_QUEUE.h"
using namespace PhysBAM;
extern PhysBAM::PTHREAD_QUEUE* pthread_queue;

unsigned long long read_tsc(){__asm__("rdtsc");}
struct timeval starttime,stoptime;
void start_timer(){gettimeofday(&starttime,NULL);read_tsc();}
void stop_timer(){gettimeofday(&stoptime,NULL);read_tsc();}
double get_time(){return (double)stoptime.tv_sec-(double)starttime.tv_sec+(double)1e-6*(double)stoptime.tv_usec-(double)1e-6*(double)starttime.tv_usec;}

int main(int argc,char *argv[]){
    
    typedef float T;
    static const int size=512;
    T* x;
    T* y;
    T* diagonal_part;

    if(argc!=2){printf("Must specify number of threads\n");exit(1);}
    int number_of_threads=atoi(argv[1]);
    printf("Using %d threads\n",number_of_threads);
    pthread_queue=new PhysBAM::PTHREAD_QUEUE(number_of_threads);  

    // Allocate data
    printf("Allocation");start_timer();

    Multiplication_And_Dot_Product_Size_Specific_Helper<T,size>::Allocate_Data(x,y,diagonal_part);

    stop_timer();printf(" [Seconds: %g]\n",get_time());

    // Initialize data
    printf("Initialization");start_timer();

    Multiplication_And_Dot_Product_Size_Specific_Helper<T,size>::Initialize_Data(x,y,diagonal_part);
    Multiplication_And_Dot_Product_Helper<T> test(size,size,size,x,y,diagonal_part,(T)1/(T)size);

    stop_timer();printf(" [Seconds: %g]\n",get_time());

    while(1){
        // Running benchmark
        printf("Running benchmark");start_timer();
	test.Run_Parallel(number_of_threads);
        
        stop_timer();printf(" [Seconds: %g]\n",get_time());
    }

    return 0;
}

//#####################################################################
// Copyright 2009, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.
//#####################################################################

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include "../multigrid_poisson_3d_optimized_kernels/Thread_Queueing/PTHREAD_QUEUE.h"
//#include "MG_PRECONDITIONED_CONJUGATE_GRADIENT.h"
//#include "MULTIGRID_POISSON.h"
//#include "MULTIGRID_POISSON_REFINEMENT.h"
//#include "MULTIGRID_POISSON_SOLVER.h"
#include "MULTIGRID_POISSON_TESTS.h"


using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;
int main(int argc,char* argv[])
{ 
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    static const int d=3;

    if(argc<3 || argc >4){
	std::cout<<"Usage: "<<argv[0]<<" <test_number> <number_of_threads> <resolution>(optional)"<<std::endl; return 1;
    }

    const int test_number=atoi(argv[1]);
    const int number_of_threads=atoi(argv[2]);
    int resolution =256;
    if(argc==4)
	resolution=atoi(argv[3]);


    pthread_queue=new PhysBAM::PTHREAD_QUEUE(number_of_threads);
    MULTIGRID_POISSON_TESTS<T,d> tests(test_number,number_of_threads,resolution);
    tests.Run();
    return 0;
}

#include "EXECUTE_HELPER.h"
#include "CACHED_ELIMINATION_MATRIX.h"
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Log/LOG.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/FACE_RANGE_ITERATOR.h>
#include "CACHED_ELIMINATION_MATRIX.h"
#include "FLUID_LAYOUT.h"
#include "FREQUENCY_TRACKER.h"
#include "JOB_SCHEDULER.h"
#include "EXECUTE_HELPER.h"
#include <stdio.h>

//these are the prototypes for the cuda_kernels we want to run
void Execute_Helper_Kernel(int);

namespace PhysBAM
{
  //variables copies over from Cached Elimination Matrix
  enum op_type
{
    op_nop,
    op_mat_inv,op_mat_mul,op_mat_add,
    op_vec_mul,op_vec_add,
    op_last
};

enum {pseudo_inv=1<<30,raw_op_mask=~pseudo_inv};

  //create structure for CUDA JOB
  struct CUDA_JOB
  {
    /* data */
    int op;
    int a[3],o;
    void cuda_execute(CACHED_ELIMINATION_MATRIX<double>*cem);
  };

//Define the struct cuda_execute based off the execute in CACHED_ELIMINATION MATRIX
void CUDA_JOB::cuda_execute(CACHED_ELIMINATION_MATRIX<double>*cem){
  //here we want to make our calls to cuda_kernels 
  //auto &vl= cem->vector_list;
  //int raw_op = this->op&raw_op_mask;
  MATRIX_MXN<double>* MO=0,*MA[3]={};

}

//Define the cuda_execute_jobs function based off the execute_job functions in CACHED_ELIMINATION_MATRIX
void cuda_execute_jobs(CACHED_ELIMINATION_MATRIX<double>*cem,int num_threads){
    cem->Compute_Job_Deps();
    cem->Eliminate_Trans();
    cem->Combine_Ops();
    cem->Simplify_Jobs();
    cem->Relabel();
    JOB_SCHEDULER<CACHED_ELIMINATION_MATRIX<double>::JOB,CACHED_ELIMINATION_MATRIX<double>> cuda_scheduler(cem);
    for(int i = 0; i < cem->jobs.m;i++){
        int id = cuda_scheduler.Add_Job(&cem->jobs(i),0);
        PHYSBAM_ASSERT(i==id);
    }
    for(auto a:cem->dep_list ){
        cuda_scheduler.Register_Dependency(a.x,a.y);
    }
    cuda_scheduler.Compute_Priority_By_Paths();
    cuda_scheduler.Execute_Jobs(num_threads);
}

void Execute_Helper(CACHED_ELIMINATION_MATRIX<double>* cem,int num_threads){
    // Schedule CUDA kernels here
    
    Execute_Helper_Kernel(2);
    //flaten the array here

    //call the cuda kernel 

    //retrieve array
  }
}

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
#include <cblas.h>
#include <lapacke.h>

//these are the prototypes for the cuda_kernels we want to run
void Execute_Helper_Kernel(int);
void cuda_kernel_vector_addition(double* result,double * vector1,double * vector2,double scalar,int size);
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

int cuda_arg_type[op_last][4] =
{
    [op_nop]={0,0,0,0},
    [op_mat_inv]={1,0,0,1},
    [op_mat_mul]={1,1,1,1},
    [op_mat_add]={1,1,0,1},
    [op_vec_mul]={2,1,2,2},
    [op_vec_add]={2,2,0,2}
};

  //create structure for CUDA JOB
  //struct CUDA_JOB:
  //{
    /* data */
  //  int op;
  //   int a[3],o;
  //  void Execute(CACHED_ELIMINATION_MATRIX<double>*cem);
  //};

void Execute(CACHED_ELIMINATION_MATRIX<double>::JOB* job,CACHED_ELIMINATION_MATRIX<double>*cem);

template<class JOB,class DATA>
struct CUDA_JOB_SCHEDULER:public JOB_SCHEDULER_CORE
{
    CUDA_JOB_SCHEDULER(DATA* user_data)
    {
        execute_job=[](void* job,void* data){Execute((JOB*)job,(DATA*)data);};
        data=user_data;
    }

    int Add_Job(JOB* job,int priority)
    {return JOB_SCHEDULER_CORE::Add_Job(job,priority);}
};

//declare function static to make it only available to this file
static void Pseudo_Inverse(MATRIX_MXN<double>& A)
{
    ARRAY<double> s(A.m),u(A.m*A.m),vt(A.m*A.m),superb(A.m);
    int ret=LAPACKE_dgesvd(LAPACK_COL_MAJOR,'S','S',A.m,A.n,A.x.Get_Array_Pointer(),A.m,
        s.Get_Array_Pointer(),u.Get_Array_Pointer(),A.m,vt.Get_Array_Pointer(),A.m,superb.Get_Array_Pointer());
    PHYSBAM_ASSERT(!ret);
    double tol=1e-16*A.m*s.Max();
    for(int i=0;i<A.m;i++)
    {
       double ss=s(i);
       if(s(i)>tol) ss=1.0/s(i);
       cblas_dscal(A.m,ss,&vt(i*A.m),1);
    }
    cblas_dgemm(CblasColMajor,CblasTrans,CblasTrans,
        A.m,A.m,A.m,1.0,vt.Get_Array_Pointer(),A.m,u.Get_Array_Pointer(),A.m,0.0,A.x.Get_Array_Pointer(),A.m);
}

//Declare function as static so it's only availble to this file
static void Inverse(MATRIX_MXN<double>& A)
{
    ARRAY<int> idiv(A.m);
    int ret_f = LAPACKE_dsytrf( LAPACK_COL_MAJOR, 'U', A.m, A.x.Get_Array_Pointer(), A.m, idiv.Get_Array_Pointer() );
    PHYSBAM_ASSERT(!ret_f);
    int ret_i = LAPACKE_dsytri( LAPACK_COL_MAJOR, 'U', A.m, A.x.Get_Array_Pointer(), A.m, idiv.Get_Array_Pointer() );
    PHYSBAM_ASSERT(!ret_i);
    for(int r=0;r<A.m;r++)
        for(int c=0;c<r;c++)
            A(r,c)=A(c,r);
}

//Declare function as static so it's only availble to this file
static void Times_MV(ARRAY<double>& v,double a,const MATRIX_MXN<double>& M,bool t,const ARRAY<double>& u,double b)
{
    cblas_dgemv(CblasColMajor,t?CblasTrans:CblasNoTrans,M.m,M.n,
        a, M.x.Get_Array_Pointer(), M.m, u.Get_Array_Pointer(), 1, b,
        v.Get_Array_Pointer(), 1);
}

static void Times_MM(MATRIX_MXN<double>& A,double sa,const MATRIX_MXN<double>& B,bool bt,const MATRIX_MXN<double>& C,bool ct,double sbc)
{
    int m=bt?B.n:B.m;
    int k=bt?B.m:B.n;
    int n=ct?C.m:C.n;
    if(A.m!=m || A.n!=n) A.Resize(m,n);
    
    cblas_dgemm( CblasColMajor, bt?CblasTrans:CblasNoTrans, ct?CblasTrans:CblasNoTrans,
        m,n,k,sbc,B.x.Get_Array_Pointer(),
        B.m, C.x.Get_Array_Pointer(), C.m,
        sa, A.x.Get_Array_Pointer(), A.m);
}

//Define the struct cuda_execute based off the execute in CACHED_ELIMINATION MATRIX
void Execute(CACHED_ELIMINATION_MATRIX<double>::JOB* job,CACHED_ELIMINATION_MATRIX<double>*cem){
  //here we want to make our calls to cuda_kernels 
 
  auto &vl= cem->vector_list;
  int raw_op = job->op&raw_op_mask;
  MATRIX_MXN<double>* MO=0,*MA[3]={};
  for(int l = 0; l < 3; l++){
    if(cuda_arg_type[raw_op][l] == 1 && job->a[l] > 1)
        MA[l] = &cem->matrix_cache.Use(job->a[l]&raw_mask);
  }
  
  if(cuda_arg_type[raw_op][3] == 1 && job->o > 1)
      MO = &cem->matrix_cache.Use(job->o);
   
    switch (raw_op)
    {
      case op_nop:break;
      case op_mat_inv: 
      {

         *MO=*MA[0];
          if(job->op&pseudo_inv) Pseudo_Inverse(*MO);
          else Inverse(*MO);
          break;
          
      }
      break;
      case op_mat_mul: 
      {
        int s0=0,s1=((job->a[1]^job->a[2])&use_neg)?-1:1;
        if(job->a[0]>=0) s0=(job->a[0]&use_neg)?-1:1;
        if(job->a[0]<0) MO->Resize(job->a[1]&use_trans?MA[1]->n:MA[1]->m,job->a[2]&use_trans?MA[2]->m:MA[2]->n);
        else if((job->a[0]&raw_mask)!=job->o) *MO=*MA[0];
        Times_MM(*MO,s0,*MA[1],job->a[1]&use_trans,*MA[2],job->a[2]&use_trans,s1);
      }
      break;
      case op_mat_add: 
      {
         if(job->a[1]&use_trans)
         {
          if(job->a[1]&use_neg) *MO=*MA[0]-MA[1]->Transposed();
          else *MO=*MA[0]+MA[1]->Transposed();
         }
         else
         {
            if(job->a[1]&use_neg) *MO=*MA[0]-*MA[1];
            else *MO=*MA[0]+*MA[1];
         }
      }
      break;
      case op_vec_mul: 
      {
        int s0=0,s1=((job->a[1]^job->a[2])&use_neg)?-1:1;
        if(job->a[0]>=0) s0=(job->a[0]&use_neg)?-1:1;
        auto& M=*MA[1];
        if(job->a[0]<0) vl(job->o).Resize(job->a[1]&use_trans?M.n:M.m);
        else if(job->a[0]!=job->o) vl(job->o)=vl(job->a[0]&raw_mask);
        Times_MV(vl(job->o),s1,M,job->a[1]&use_trans,vl(job->a[2]&raw_mask),s0);
      }
      break;
      case op_vec_add:
       {
        PHYSBAM_ASSERT(vl(job->a[0]).m==vl(job->a[1]).m);  
        if(job->a[1] & use_neg){
          //subtraction between two vectors
          LOG::printf("Vector subtraction");
          cuda_kernel_vector_addition(vl(job->o).Get_Array_Pointer(),vl(job->a[0]).Get_Array_Pointer(),vl(job->a[1]).Get_Array_Pointer(),-1,vl(job->a[0]).m);
        }
        else{
          LOG::printf("Vector addition");
          cuda_kernel_vector_addition(vl(job->o).Get_Array_Pointer(),vl(job->a[0]).Get_Array_Pointer(),vl(job->a[1]).Get_Array_Pointer(),1,vl(job->a[0]).m);
        }
       }
        break;
      default: PHYSBAM_FATAL_ERROR();
    }
  
    for(int l = 0; l < 3; l++){
        if(cuda_arg_type[raw_op][l] > 0 && job->a[l] > 1){
          bool rem=!--cem->data_refs[cuda_arg_type[raw_op][l]-1][job->a[l]&raw_mask];
          if(cuda_arg_type[raw_op][l]==1) cem->matrix_cache.Release(job->a[l]&raw_mask,rem);
          else if(rem) vl(job->a[l]&raw_mask).Clean_Memory();
        }
    }
    if(cuda_arg_type[raw_op][3]==1 && job->o>1)
        cem->matrix_cache.Release(job->o,false);
}

//Define the cuda_execute_jobs function based off the execute_job functions in CACHED_ELIMINATION_MATRIX
void cuda_execute_jobs(CACHED_ELIMINATION_MATRIX<double>*cem,int num_threads){
    cem->Compute_Job_Deps();
    cem->Eliminate_Trans();
    cem->Combine_Ops();
    cem->Simplify_Jobs();
    cem->Relabel();
    CUDA_JOB_SCHEDULER<CACHED_ELIMINATION_MATRIX<double>::JOB,CACHED_ELIMINATION_MATRIX<double>> cuda_scheduler(cem);
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
    
    cuda_execute_jobs(cem,num_threads);
  }
}

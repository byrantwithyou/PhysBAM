//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __CACHED_ELIMINATION_MATRIX__
#define __CACHED_ELIMINATION_MATRIX__
#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include <Core/Vectors/VECTOR.h>
#include "COMMON.h"
#include <mutex>

namespace PhysBAM{

template<class T> class MATRIX_MXN;
template<class T>
struct CACHED_ELIMINATION_MATRIX
{
    struct MATRIX_INFO
    {
        MATRIX_MXN<T> M;
        bool sym;
        ARRAY<int> prod_list;
    };

    ARRAY<MATRIX_INFO> block_list;
    ARRAY<ARRAY<T> > vector_list;
    ARRAY<int> orig_sizes;
    HASHTABLE<VECTOR<int,2>,int> blocks_to_canonical_block_id;

    HASHTABLE<VECTOR<int,3>,int> cached_ops;
    HASHTABLE<ARRAY<int>,int> prod_lookup;

    struct MATRIX_BLOCK
    {
        int c;
        int matrix_id;
    };

    ARRAY<ARRAY<MATRIX_BLOCK> > rows;
    ARRAY<ARRAY<T> > orig_rhs;
    ARRAY<int> rhs;
    ARRAY<ARRAY<T> > test_sol;

    ARRAY<bool> valid_row;
    ARRAY<int> elimination_order;
    bool quiet,pinv_last_blk=false;

    // Note: Jobs are constructed in SSA form
    struct JOB
    {
        int op;
        int a[3],o;

        void Execute(CACHED_ELIMINATION_MATRIX<T>* cem);
    };
    ARRAY<JOB> jobs;
    ARRAY<int> provider[2];
    ARRAY<ARRAY<int> > user[2];
    int num_orig_blocks;
    int num_orig_vectors;
    ARRAY<VECTOR<int,2> > dep_list;
    std::mutex mtx;
    ARRAY<int> data_refs[2];

    template<class F>
    void Foreach_Single_User_Pair(F f)
    {
        for(int i=0;i<2;i++)
            for(int j=0;j<user[i].m;j++)
                if(user[i](j).m==1 && provider[i](j)>=0)
                    f(provider[i](j),user[i](j)(0));
    }
    
    void Fill_Blocks(ARRAY<VECTOR<int,2>,DOF_ID>& dof_map,const ARRAY<TRIPLE<DOF_ID,DOF_ID,CODE_ID> >& coded_entries,
        const ARRAY<T,CODE_ID>& code_values,const ARRAY<T,DOF_ID>& rhs_vector);
    void Begin_Fill_Blocks();
    void End_Fill_Blocks();
    void Fill_Orig_Rows();
    void Eliminate_Row(int r);
    int& Get_Block(int r,int c);
    int Get_Block_Lazy(int r,int c) const;
    int Compute_Inv(int m,bool pinv);
    int Compute_Mul(int a,int b);
    int Compute_Add(int a,int b);
    int Compute_Elim(int a,int b,int c);
    void Print_Full() const;
    void Print_Current() const;
    void Back_Solve();
    void Unpack_Vector(ARRAY<VECTOR<int,2>,DOF_ID>& dof_map,ARRAY<ARRAY<T> >& u,const ARRAY<T,DOF_ID>& v);
    void Pack_Vector(ARRAY<VECTOR<int,2>,DOF_ID>& dof_map,ARRAY<T,DOF_ID>& v,const ARRAY<ARRAY<T> >& u);
    void Pack_Vector(ARRAY<VECTOR<int,2>,DOF_ID>& dof_map,ARRAY<T,DOF_ID>& v,const ARRAY<int>& u);
    int Matrix_Times(int m,int in);
    int Sub_Times(int out,int m,int in);
    int Transposed(int a) const;
    int Negate(int a) const;
    bool Symmetric(int a) const;
    ARRAY<int> Transposed(const ARRAY<int>& a) const;
    ARRAY<int> Prod_List(int a) const;
    void Reduce_Rows_By_Frequency(int begin,int end,int fill_limit);
    void Full_Reordered_Elimination();
    void Execute_Jobs(int num_threads);
    void Compute_Job_Deps();
    void Eliminate_Trans();
    void Simplify_Jobs();
    void Combine_Ops();
    void Relabel();
    void Combine_With_Next_Job(int j);
    int Provides_Argument(const JOB& j,int a) const;
    const ARRAY<int>& Uses_Output(const JOB& j) const;
    void Remove_Job_Deps(int j);
    void Add_Job_Deps(int j);
    int Create_Matrix_Block(bool sym);
    void Add_Block_Matrix_Entry(int r,int c,int matrix_id);
};
}
#endif

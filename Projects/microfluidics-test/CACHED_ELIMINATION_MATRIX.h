//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __CACHED_ELIMINATION_MATRIX__
#define __CACHED_ELIMINATION_MATRIX__
#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Vectors/VECTOR.h>

namespace PhysBAM{

template<class T> class MATRIX_MXN;
template<class T>
struct CACHED_ELIMINATION_MATRIX
{
    enum {zero_block=0,id_block=1,invalid_block=-1,use_trans=1<<30,is_vec=1<<29};

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
    
    enum {op_inv,op_mul,op_sub,op_Av,op_vec_neg,op_vec_sub,op_u_sub_Av,free_mat,free_vec};

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
    bool quiet;
    
    struct JOB
    {
        int op;
        int a,b,c,o;
        int is_vec_mask;
        int priority;
        int job_id;
        ARRAY<int> users;

        void Execute(CACHED_ELIMINATION_MATRIX<T>* cem);
    };
    ARRAY<JOB> jobs;

    void Fill_Blocks(ARRAY<VECTOR<int,2> >& dof_map,const ARRAY<VECTOR<int,3> >& coded_entries,
        const ARRAY<T>& code_values,const ARRAY<T>& rhs_vector);
    void Fill_Orig_Rows();
    void Eliminate_Row(int r);
    int& Get_Block(int r,int c);
    int Get_Block_Lazy(int r,int c) const;
    int Compute_Inv(int m);
    int Compute_Mul(int a,int b);
    int Compute_Sub(int a,int b);
    int Compute_Elim(int a,int b,int c);
    void Print_Full() const;
    void Print_Current() const;
    void Back_Solve();
    void Unpack_Vector(ARRAY<VECTOR<int,2> >& dof_map,ARRAY<ARRAY<T> >& u,const ARRAY<T>& v);
    void Pack_Vector(ARRAY<VECTOR<int,2> >& dof_map,ARRAY<T>& v,const ARRAY<ARRAY<T> >& u);
    void Pack_Vector(ARRAY<VECTOR<int,2> >& dof_map,ARRAY<T>& v,const ARRAY<int>& u);
    int Matrix_Times(int m,int in);
    int Sub_Times(int out,int m,int in);
    int Transposed(int a) const;
    bool Symmetric(int a) const;
    ARRAY<int> Transposed(const ARRAY<int>& a) const;
    ARRAY<int> Prod_List(int a) const;
    void Reduce_Rows_By_Frequency(int begin,int end,int fill_limit);
    void Full_Reordered_Elimination();
    void Execute_Jobs(int num_threads);
};
}
#endif

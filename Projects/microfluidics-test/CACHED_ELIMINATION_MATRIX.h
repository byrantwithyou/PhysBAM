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
template<class TV> struct FLUID_LAYOUT;
template<class T>
struct CACHED_ELIMINATION_MATRIX
{
    enum {zero_block=0,id_block=1,invalid_block=-1,use_trans=1<<30};

    struct MATRIX_INFO
    {
        MATRIX_MXN<T> M;
        // symmetric if left-multiplied by this matrix
        // id_block if symmetric
        // zero_block if no special properties
        bool sym;
        ARRAY<int> prod_list;
    };

    ARRAY<MATRIX_INFO> block_list;
    ARRAY<int> orig_sizes;
    HASHTABLE<VECTOR<int,2>,int> blocks_to_canonical_block_id;
    
    enum {op_inv,op_mul,op_sub};

    HASHTABLE<VECTOR<int,3>,int> cached_ops;
    HASHTABLE<ARRAY<int>,int> prod_lookup;

    struct MATRIX_BLOCK
    {
        int c;
        int matrix_id;
    };

    ARRAY<ARRAY<MATRIX_BLOCK> > rows;
    ARRAY<ARRAY<T> > rhs;
    ARRAY<ARRAY<T> > test_sol;

    ARRAY<bool> valid_row;
    ARRAY<int> elimination_order;
    
    MATRIX_MXN<T>& Get_Orig_By_Blocks(int b0,int b1);

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
    void Test_State(const char* str) const;
    void Unpack_Vector(ARRAY<VECTOR<int,2> >& dof_map,ARRAY<ARRAY<T> >& u,const ARRAY<T>& v);
    void Pack_Vector(ARRAY<VECTOR<int,2> >& dof_map,ARRAY<T>& v,const ARRAY<ARRAY<T> >& u);
    void Add_Times(ARRAY<ARRAY<T> >& out,const ARRAY<ARRAY<T> >& in) const;
    void Add_Times(ARRAY<T>& out,T a,int m,const ARRAY<T>& in,T b) const;
    int Transposed(int a) const;
    bool Symmetric(int a) const;
    ARRAY<int> Transposed(const ARRAY<int>& a) const;
    ARRAY<int> Prod_List(int a) const;
};

template<class T,class TV>
void Setup_Block_Map(CACHED_ELIMINATION_MATRIX<T>& cem,const FLUID_LAYOUT<TV>& fl);
}
#endif

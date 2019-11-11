//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MATRIX_CONSTRUCTION_FEM__
#define __MATRIX_CONSTRUCTION_FEM__
#include "BLOCK_MATRIX.h"
#include "BLOCK_VECTOR.h"
#include "COMPONENT_LAYOUT_FEM.h"
#include "MANAGED_MATRIX.h"
#include <functional>

namespace PhysBAM{

template<class T> struct CACHED_ELIMINATION_MATRIX;
template<class T> struct SYSTEM_MATRIX_HELPER;

template<class TV>
struct MATRIX_CONSTRUCTION_FEM
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<T,2> TV2;
    typedef VECTOR<int,2> IV2;
    typedef VECTOR<int,3> IV3;

    COMPONENT_LAYOUT_FEM<T>& cl;
    T mu=1;

    MATRIX_CONSTRUCTION_FEM(COMPONENT_LAYOUT_FEM<T>& cl,const std::string& cache_pattern,int cache_size=-1);
    ~MATRIX_CONSTRUCTION_FEM();

    MANAGED_MATRIX<T,BLOCK_MATRIX<TV> > canonical_matrix_cache;
    ARRAY<int,REFERENCE_BLOCK_ID> reference_matrix;
    ARRAY<std::function<void(MATRIX_MXN<T>&)>,REFERENCE_BLOCK_ID> diagonal_system_blocks;
    ARRAY<std::function<void(MATRIX_MXN<T>&)>,REFERENCE_CONNECTION_ID> regular_system_blocks;
    ARRAY<ARRAY<std::function<void(MATRIX_MXN<T>&)>,RID_ID>,REFERENCE_IRREGULAR_ID> irregular_system_blocks;

    HASHTABLE<CANONICAL_BLOCK<T>*,int> canonical_block_matrices;
    ARRAY<BLOCK_VECTOR<TV>,BLOCK_ID> rhs_block_list;

    const BLOCK_MATRIX<TV>& Canonical_Matrix(BLOCK_ID b);
    void Release_Canonical_Matrix(BLOCK_ID b);
    void Compute_Matrix_Blocks();
    void Compute_RHS();
    void Copy_To_CEM(CACHED_ELIMINATION_MATRIX<T>& cem);
    void Fill_Canonical_Block_Matrix(BLOCK_MATRIX<TV>& mat,const REFERENCE_BLOCK_DATA& rb);
    void Fill_Block_Matrix(BLOCK_MATRIX<TV>& M,const REFERENCE_BLOCK_DATA& rd);
    void Fill_Connection_Matrix(BLOCK_MATRIX<TV>& M,const REFERENCE_CONNECTION_DATA& cd);
    void Fill_Irregular_Connection_Matrix(BLOCK_MATRIX<TV>& M,const REFERENCE_IRREGULAR_DATA& ri,RID_ID j);
    void Copy_Matrix_Data(BLOCK_MATRIX<TV>& A,BLOCK_ID b,
        const DOF_PAIRS& dpa,const DOF_PAIRS& dpb,BLOCK_ID ar,BLOCK_ID ac);
    void Copy_Vector_Data(const BLOCK_VECTOR<TV>& B,BLOCK_ID a,BLOCK_ID b,const DOF_PAIRS& dp);
    void Init_Block_Matrix(BLOCK_MATRIX<TV>& M,BLOCK_ID a,BLOCK_ID b,bool compressed) const;
    void Init_Block_Vector(BLOCK_VECTOR<TV>& M,BLOCK_ID b,bool compressed) const;
    void Times_U_Dot_V(BLOCK_ID b,BLOCK_VECTOR<TV>& v,const BLOCK_VECTOR<TV>& u) const;
    void Times_P_U(BLOCK_ID b,BLOCK_VECTOR<TV>& w,const ARRAY<T>& div_v,const ARRAY<T>& div_e) const;
    void Times_Line_Integral_U_Dot_V(BLOCK_ID b,INTERVAL<int> bc_e,BLOCK_VECTOR<TV>& w,const BLOCK_VECTOR<TV>& u) const;
    void Apply_To_RHS(BLOCK_ID b,const BLOCK_VECTOR<TV>& w);
    void Transform_To_World_Space(BLOCK_MATRIX<TV>& M,const BLOCK_MATRIX<TV>& B,BLOCK_ID a,BLOCK_ID b) const;
    void Transform_Solution(const CACHED_ELIMINATION_MATRIX<T>& cem,bool inverse,bool transpose);
    void Dump_World_Space_System(ARRAY<int,BLOCK_ID> first[3],int size,SPARSE_MATRIX_FLAT_MXN<T>& SM) const;
    void Dump_World_Space_Vector(ARRAY<int,BLOCK_ID> first[3],int size,ARRAY<T>& vec) const;
    void Dump_World_Space_System() const;
    void Dump_World_Space_Vector(const char* name) const;
    void Dump_Matrix_Block(SYSTEM_MATRIX_HELPER<T>& h,ARRAY<int,BLOCK_ID> first[3],const BLOCK_MATRIX<TV>& M,BLOCK_ID b0,BLOCK_ID b1) const;
    int Compute_Global_Dof_Mapping(ARRAY<int,BLOCK_ID> first[3]) const;
    std::tuple<BLOCK_ID,int,int,int> Inverse_DOF_Lookup(int global_dof) const; // BLOCK_ID,vep,dof,dim
    void Print_Statistics() const;
};

}
#endif

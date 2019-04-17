//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MATRIX_CONSTRUCTION_FEM__
#define __MATRIX_CONSTRUCTION_FEM__
#include "BLOCK_MATRIX.h"
#include "BLOCK_VECTOR.h"
#include "COMPONENT_LAYOUT_FEM.h"

namespace PhysBAM{

template<class T> struct CACHED_ELIMINATION_MATRIX;
template<class T> struct SYSTEM_MATRIX_HELPER;

template<class TV>
struct MATRIX_CONSTRUCTION_FEM
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,2> IV2;
    typedef VECTOR<int,3> IV3;

    COMPONENT_LAYOUT_FEM<T>& cl;
    T mu=1;

    MATRIX_CONSTRUCTION_FEM(COMPONENT_LAYOUT_FEM<T>& cl);
    ~MATRIX_CONSTRUCTION_FEM();
    
    ARRAY<BLOCK_MATRIX<T>,REFERENCE_BLOCK_ID> diagonal_system_blocks;
    ARRAY<BLOCK_MATRIX<T>,REFERENCE_CONNECTION_ID> regular_system_blocks;
    ARRAY<ARRAY<BLOCK_MATRIX<T>,RID_ID>,REFERENCE_IRREGULAR_ID> irregular_system_blocks;

    HASHTABLE<CANONICAL_BLOCK<T>*,BLOCK_MATRIX<T> > canonical_block_matrices;
    ARRAY<BLOCK_VECTOR<T>,BLOCK_ID> rhs_block_list;

    void Compute_Matrix_Blocks();
    void Compute_RHS();
    void Copy_To_CEM(CACHED_ELIMINATION_MATRIX<T>& cem);
    void Fill_Canonical_Block_Matrix(BLOCK_MATRIX<T>& mat,const CANONICAL_BLOCK<T>* cb);
    void Fill_Block_Matrix(BLOCK_MATRIX<T>& M,const REFERENCE_BLOCK_DATA& rd);
    void Fill_Connection_Matrix(BLOCK_MATRIX<T>& M,const REFERENCE_CONNECTION_DATA& cd);
    void Fill_Irregular_Connection_Matrix(ARRAY<BLOCK_MATRIX<T>,RID_ID>& M,const REFERENCE_IRREGULAR_DATA& ri);
    void Copy_Matrix_Data(BLOCK_MATRIX<T>& A,BLOCK_ID b,
        const DOF_PAIRS& dpa,const DOF_PAIRS& dpb,BLOCK_ID ar,BLOCK_ID ac) const;
    void Copy_Vector_Data(const BLOCK_VECTOR<T>& B,BLOCK_ID b,const DOF_PAIRS& dp);
    void Init_Block_Matrix(BLOCK_MATRIX<T>& M,BLOCK_ID a,BLOCK_ID b) const;
    void Init_Block_Vector(BLOCK_VECTOR<T>& M,BLOCK_ID b) const;
    void Init_Block_Vector(BLOCK_VECTOR<T>& M,const CANONICAL_BLOCK<T>* cb) const;
    void Times_U_Dot_V(BLOCK_ID b,BLOCK_VECTOR<T>& v,const BLOCK_VECTOR<T>& u) const;
    void Times_P_U(BLOCK_ID b,BLOCK_VECTOR<T>& w,const ARRAY<T>& div_v,const ARRAY<T>& div_e) const;
    void Times_Line_Integral_U_Dot_V(BLOCK_ID b,INTERVAL<int> bc_e,BLOCK_VECTOR<T>& w,const BLOCK_VECTOR<T>& u) const;
    void Apply_To_RHS(BLOCK_ID b,const BLOCK_VECTOR<T>& w);
    void Transform_To_World_Space(BLOCK_MATRIX<T>& M,const BLOCK_MATRIX<T>& B,BLOCK_ID a,BLOCK_ID b) const;
    void Transform_Solution(const CACHED_ELIMINATION_MATRIX<T>& cem,bool inverse,bool transpose);
    void Dump_World_Space_System() const;
    void Dump_World_Space_Vector(const char* name) const;
    void Dump_Matrix_Block(SYSTEM_MATRIX_HELPER<T>& h,ARRAY<int,BLOCK_ID> first[3],const BLOCK_MATRIX<T>& M,BLOCK_ID b0,BLOCK_ID b1) const;
    int Compute_Global_Dof_Mapping(ARRAY<int,BLOCK_ID> first[3]) const;
};

}
#endif

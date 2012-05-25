//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> template<int d> void SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<TV>::
Initialize(const BASIS_STENCIL_UNIFORM<TV,d>& s,CELL_MANAGER_COLOR<TV>& cm_input,CELL_DOMAIN_INTERFACE_COLOR<TV> &cdi_input)
{
    cm=&cm_input;
    cdi=&cdi_input;

    for(int i=0;i<s.diced.m;i++) flat_diff.Append(cdi->Flatten_Diff(s.diced(i).index_offset));

    flat_diff.Prune_Duplicates();
    flat_diff.Sort();

    data.Resize(cdi->colors);
}
//#####################################################################
// Function Mark_Active_Cells
//#####################################################################
template<class TV> void SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<TV>::
Mark_Active_Cells(T tol)
{
    for(int c=0;c<cdi->colors;c++){
        MATRIX_MXN<T>& d=data(c);
        for(int l=0;l<d.m;l++)
            for(int k=0;k<d.n;k++)
                if(abs(d(l,k))>tol)
                    cm->Set_Active(cdi->flat_base_scalar(l)+flat_diff(k),c);
                else d(l,k)=0;}
}
//#####################################################################
// Function Build_Matrix
//#####################################################################
template<class TV> void SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<TV>::
Build_Matrix(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& matrix)
{
    if(!cdi->wrap) PHYSBAM_FATAL_ERROR();
    matrix.Resize(cdi->colors);

    for(int c=0;c<cdi->colors;c++){
        SPARSE_MATRIX_FLAT_MXN<T>& M=matrix(c);
        ARRAY<int>& comp_n=cm->compressed(c);
        int m=cdi->constraint_base_scalar;
        int n=cm->dofs(c);
        M.Reset(n);
        M.offsets.Resize(m+1);
        M.m=m;

        MATRIX_MXN<T>& d=data(c);
        for(int row=0;row<d.m;row++,row++){
            ARRAY<SPARSE_MATRIX_ENTRY<T> > entries;
            for(int j=0;j<d.n;j++){
                T value=d(row,j);
                if(value){
                    int column=comp_n(cdi->flat_base_scalar(row)+flat_diff(j));
                    M.offsets(row+1)++;
                    entries.Append(SPARSE_MATRIX_ENTRY<T>(column,value));}}
            if(cdi->Is_Boundary_Constraint_Scalar(row)) entries.Sort();
            M.A.Append_Elements(entries);}

        for(int i=0;i<M.offsets.m-1;i++) M.offsets(i+1)+=M.offsets(i);}
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<TV>::
Resize()
{
    for(int c=0;c<cdi->colors;c++)
        data(c).Resize(cdi->flat_base_scalar.m,flat_diff.m);
}
template class SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<VECTOR<float,2> >;
template class SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<VECTOR<float,3> >;
template void SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<VECTOR<float,2> >::Initialize<1>(BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,
    CELL_MANAGER_COLOR<VECTOR<float,2> >&,CELL_DOMAIN_INTERFACE_COLOR<VECTOR<float,2> >&);
template void SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<VECTOR<float,3> >::Initialize<1>(BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,
    CELL_MANAGER_COLOR<VECTOR<float,3> >&,CELL_DOMAIN_INTERFACE_COLOR<VECTOR<float,3> >&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<VECTOR<double,2> >;
template class SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<VECTOR<double,3> >;
template void SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<VECTOR<double,2> >::Initialize<1>(BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,
    CELL_MANAGER_COLOR<VECTOR<double,2> >&,CELL_DOMAIN_INTERFACE_COLOR<VECTOR<double,2> >&);
template void SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<VECTOR<double,3> >::Initialize<1>(BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,
    CELL_MANAGER_COLOR<VECTOR<double,3> >&,CELL_DOMAIN_INTERFACE_COLOR<VECTOR<double,3> >&);
#endif

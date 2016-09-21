//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Matrices/MATRIX_MXN.h>
#include <Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <Geometry/Finite_Elements/CELL_MANAGER_COLOR.h>
#include <Geometry/Finite_Elements/SYSTEM_SURFACE_BLOCK_HELPER_COLOR.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> template<int d> void SYSTEM_SURFACE_BLOCK_HELPER_COLOR<TV>::
Initialize(CELL_DOMAIN_INTERFACE_COLOR<TV> &cdi_input,const BASIS_STENCIL_UNIFORM<TV,d>& s,CELL_MANAGER_COLOR<TV>& cm_input)
{
    cm=&cm_input;
    cdi=&cdi_input;

    for(int i=0;i<s.diced.m;i++) flat_diff.Append(cdi->Flatten_Diff(s.diced(i).index_offset));

    Prune_Duplicates(flat_diff);
    flat_diff.Sort();

    for(int i=0;i<TV::m;i++){
        data(i).Resize(cdi->colors);
        rhs_data(i).Resize(cdi->colors);}
}
//#####################################################################
// Function Mark_Active_Cells
//#####################################################################
template<class TV> void SYSTEM_SURFACE_BLOCK_HELPER_COLOR<TV>::
Mark_Active_Cells(T tol)
{
    for(int i=0;i<TV::m;i++)
        for(int c=0;c<cdi->colors;c++){
            MATRIX_MXN<T>& d=data(i)(c);
            for(int l=0;l<d.m;l++)
                for(int k=0;k<d.n;k++)
                    if(abs(d(l,k))>tol)
                        cm->Set_Active((*cdi->flat_base(i))(l)+flat_diff(k),c);
                    else d(l,k)=0;}
}
//#####################################################################
// Function Build_Matrix
//#####################################################################
template<class TV> void SYSTEM_SURFACE_BLOCK_HELPER_COLOR<TV>::
Build_Matrix(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& matrix,ARRAY<T>& constraint_rhs)
{
    matrix.Resize(cdi->colors);

    for(int c=0;c<cdi->colors;c++){
        SPARSE_MATRIX_FLAT_MXN<T>& M=matrix(c);
        const ARRAY<int>& comp_n=cm->compressed(c);
        int m=cdi->total_number_of_surface_constraints;
        int n=cm->dofs(c);
        M.Reset(n);
        M.offsets.Resize(m+1);
        M.m=m;
        constraint_rhs.Resize(m);

        int row=0;
        for(int orientation=0;orientation<TV::m;orientation++){
            const MATRIX_MXN<T>& d=data(orientation)(c);
            const ARRAY<T>& rd=rhs_data(orientation)(c);
            for(int i=0;i<d.m;i++,row++){
                constraint_rhs(row)+=rd(i);
                ARRAY<SPARSE_MATRIX_ENTRY<T> > entries;
                for(int j=0;j<d.n;j++){
                    T value=d(i,j);
                    if(value){
                        int column=comp_n((*cdi->flat_base(orientation))(i)+flat_diff(j));
                        M.offsets(row+1)++;
                        entries.Append(SPARSE_MATRIX_ENTRY<T>(column,value));}}
                if(cdi->Is_Boundary_Constraint(i,orientation)) entries.Sort();
                M.A.Append_Elements(entries);}}

        for(int i=0;i<M.offsets.m-1;i++) M.offsets(i+1)+=M.offsets(i);}
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void SYSTEM_SURFACE_BLOCK_HELPER_COLOR<TV>::
Resize()
{
    for(int i=0;i<TV::m;i++)
        for(int c=0;c<cdi->colors;c++){
            data(i)(c).Resize(cdi->flat_base(i)->m,flat_diff.m);
            rhs_data(i)(c).Resize(cdi->flat_base(i)->m);}
}
namespace PhysBAM{
template class SYSTEM_SURFACE_BLOCK_HELPER_COLOR<VECTOR<float,2> >;
template class SYSTEM_SURFACE_BLOCK_HELPER_COLOR<VECTOR<float,3> >;
template class SYSTEM_SURFACE_BLOCK_HELPER_COLOR<VECTOR<double,2> >;
template class SYSTEM_SURFACE_BLOCK_HELPER_COLOR<VECTOR<double,3> >;
template void SYSTEM_SURFACE_BLOCK_HELPER_COLOR<VECTOR<double,2> >::Initialize<1>(CELL_DOMAIN_INTERFACE_COLOR<VECTOR<double,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,CELL_MANAGER_COLOR<VECTOR<double,2> >&);
template void SYSTEM_SURFACE_BLOCK_HELPER_COLOR<VECTOR<double,3> >::Initialize<1>(CELL_DOMAIN_INTERFACE_COLOR<VECTOR<double,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,CELL_MANAGER_COLOR<VECTOR<double,3> >&);
template void SYSTEM_SURFACE_BLOCK_HELPER_COLOR<VECTOR<float,2> >::Initialize<1>(CELL_DOMAIN_INTERFACE_COLOR<VECTOR<float,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,CELL_MANAGER_COLOR<VECTOR<float,2> >&);
template void SYSTEM_SURFACE_BLOCK_HELPER_COLOR<VECTOR<float,3> >::Initialize<1>(CELL_DOMAIN_INTERFACE_COLOR<VECTOR<float,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,CELL_MANAGER_COLOR<VECTOR<float,3> >&);
}

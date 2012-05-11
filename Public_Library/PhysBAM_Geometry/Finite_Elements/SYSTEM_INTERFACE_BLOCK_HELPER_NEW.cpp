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
#include <PhysBAM_Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_NEW.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER_NEW.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_INTERFACE_BLOCK_HELPER_NEW.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> template<int d> void SYSTEM_INTERFACE_BLOCK_HELPER_NEW<TV>::
Initialize(const BASIS_STENCIL_UNIFORM<TV,d>& s,CELL_MANAGER_NEW<TV>& cm_input,CELL_DOMAIN_INTERFACE_NEW<TV> &cdi_input)
{
    cm=&cm_input;
    cdi=&cdi_input;

    for(int i=0;i<s.diced.m;i++) flat_diff.Append(cdi->Flatten_Diff(s.diced(i).index_offset));

    flat_diff.Sort();
    flat_diff.Prune_Duplicates();

    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            data(i)[s].Resize(cdi->interface_dofs*TV::m,flat_diff.m);
}
//#####################################################################
// Function Mark_Active_Cells
//#####################################################################
template<class TV> void SYSTEM_INTERFACE_BLOCK_HELPER_NEW<TV>::
Mark_Active_Cells(T tol)
{
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++){
            MATRIX_MXN<T>& d=data(i)[s];
            for(int l=0;l<d.m;l++)
                for(int k=0;k<d.n;k++)
                    if(abs(d(l,k))>tol)
                        cm->Set_Active(cdi->flat_base(l)+flat_diff(k),s);
                    else d(l,k)=0;}
}
//#####################################################################
// Function Build_Matrix
//#####################################################################
template<class TV> void SYSTEM_INTERFACE_BLOCK_HELPER_NEW<TV>::
Build_Matrix(VECTOR<SPARSE_MATRIX_FLAT_MXN<T>,2>& matrix)
{
    if(!cdi->periodic_bc) PHYSBAM_FATAL_ERROR();

    for(int s=0;s<2;s++){
        SPARSE_MATRIX_FLAT_MXN<T>& M=matrix(s);
        ARRAY<int>& comp_n=cm->compressed[s];
        int m=cdi->interface_dofs*TV::m;
        int n=cm->dofs[s];
        M.Reset(n);
        M.offsets.Resize(m+1);
        M.m=m;

        int row=0;
        for(int orientation=0;orientation<TV::m;orientation++){
            MATRIX_MXN<T>& d=data(orientation)[s];
            for(int i=0;i<cdi->interface_dofs;i++,row++){
                ARRAY<SPARSE_MATRIX_ENTRY<T> > entries;
                for(int j=0;j<d.n;j++){
                    T value=d(i,j);
                    if(value){
                        int column=comp_n(cdi->flat_base(i)+flat_diff(j));
                        M.offsets(row+1)++;
                        entries.Append(SPARSE_MATRIX_ENTRY<T>(column,value));}}
                if(cdi->Is_Boundary_Cut_Cell(i)) entries.Sort();
                M.A.Append_Elements(entries);}}

        for(int i=0;i<M.offsets.m-1;i++) M.offsets(i+1)+=M.offsets(i);}
}
template class SYSTEM_INTERFACE_BLOCK_HELPER_NEW<VECTOR<float,2> >;
template class SYSTEM_INTERFACE_BLOCK_HELPER_NEW<VECTOR<float,3> >;
template void SYSTEM_INTERFACE_BLOCK_HELPER_NEW<VECTOR<float,2> >::Initialize<1>(BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,
    CELL_MANAGER_NEW<VECTOR<float,2> >&,CELL_DOMAIN_INTERFACE_NEW<VECTOR<float,2> >&);
template void SYSTEM_INTERFACE_BLOCK_HELPER_NEW<VECTOR<float,3> >::Initialize<1>(BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,
    CELL_MANAGER_NEW<VECTOR<float,3> >&,CELL_DOMAIN_INTERFACE_NEW<VECTOR<float,3> >&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SYSTEM_INTERFACE_BLOCK_HELPER_NEW<VECTOR<double,2> >;
template class SYSTEM_INTERFACE_BLOCK_HELPER_NEW<VECTOR<double,3> >;
template void SYSTEM_INTERFACE_BLOCK_HELPER_NEW<VECTOR<double,2> >::Initialize<1>(BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,
    CELL_MANAGER_NEW<VECTOR<double,2> >&,CELL_DOMAIN_INTERFACE_NEW<VECTOR<double,2> >&);
template void SYSTEM_INTERFACE_BLOCK_HELPER_NEW<VECTOR<double,3> >::Initialize<1>(BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,
    CELL_MANAGER_NEW<VECTOR<double,3> >&,CELL_DOMAIN_INTERFACE_NEW<VECTOR<double,3> >&);
#endif

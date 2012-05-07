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
#include <PhysBAM_Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK_HELPER.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> template<int d0,int d1> void SYSTEM_VOLUME_BLOCK_HELPER<TV>::
Initialize(const BASIS_STENCIL_UNIFORM<TV,d0>& s0,const BASIS_STENCIL_UNIFORM<TV,d1>& s1,
    CELL_MANAGER<TV>& cm0_input,CELL_MANAGER<TV>& cm1_input,CELL_DOMAIN_INTERFACE<TV> &cdi_input)
{
    cm0=&cm0_input;
    cm1=&cm1_input;
    cdi=&cdi_input;

    for(int i=0;i<s0.diced.m;i++)
        for(int j=0;j<s1.diced.m;j++){
            int overlap=s0.diced(i).subcell&s1.diced(j).subcell;
            if(overlap){
                const typename BASIS_STENCIL_UNIFORM<TV,d0>::DICED& diced0=s0.diced(i);
                const typename BASIS_STENCIL_UNIFORM<TV,d1>::DICED& diced1=s1.diced(j);
                flat_diff.Append(cdi->Flatten_Diff(diced1.index_offset-diced0.index_offset));}}

    flat_diff.Sort();
    flat_diff.Prune_Duplicates();

    for(int s=0;s<2;s++) data[s].Resize(cdi->flat_size,flat_diff.m);
}
//#####################################################################
// Function Mark_Active_Cells
//#####################################################################
template<class TV> void SYSTEM_VOLUME_BLOCK_HELPER<TV>::
Mark_Active_Cells(T tol)
{
    if(cdi->periodic_bc) // add ghost rows to material ones 
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(cdi->grid,cdi->padding,GRID<TV>::GHOST_REGION,-1);it.Valid();it.Next()){
            int i=cdi->Flatten(it.index);
            int r=cdi->remap(i);
            for(int s=0;s<2;s++)
                for(int k=0;k<data[s].n;k++){
                    data[s](r,k)+=data[s](i,k);
                    data[s](i,k)=0;}}
    for(int s=0;s<2;s++)
        for(int l=0;l<data[s].m;l++)
            for(int k=0;k<data[s].n;k++)
                if(abs(data[s](l,k))>tol){
                    cm0->Set_Active(l,s);
                    cm1->Set_Active(l+flat_diff(k),s);}
                else data[s](l,k)=0;
}
//#####################################################################
// Function Build_Matrix
//#####################################################################
template<class TV> void SYSTEM_VOLUME_BLOCK_HELPER<TV>::
Build_Matrix(VECTOR<SPARSE_MATRIX_FLAT_MXN<T>,2>& matrix)
{
    if(!cdi->periodic_bc) PHYSBAM_FATAL_ERROR();

    for(int s=0;s<2;s++){
        SPARSE_MATRIX_FLAT_MXN<T>& M=matrix(s);
        MATRIX_MXN<T>& d=data[s];
        ARRAY<int>& comp_m=cm0->compressed[s];
        ARRAY<int>& comp_n=cm1->compressed[s];
        int m=cm0->dofs[s];
        int n=cm1->dofs[s];
            
        M.Reset(n);
        M.offsets.Resize(m+1);
        M.m=m;
            
        for(int i=0;i<d.m;i++){
            if(cdi->Is_Outside_Cell(i)) continue;
            int row=comp_m(i);
            if(row>=0){
                ARRAY<SPARSE_MATRIX_ENTRY<T> > entries;
                for(int j=0;j<d.n;j++){
                    T value=d(i,j);
                    if(value){
                        int column=comp_n(i+flat_diff(j));
                        M.offsets(row+1)++;
                        entries.Append(SPARSE_MATRIX_ENTRY<T>(column,value));}}
                if(cdi->Is_Boundary_Cell(i)) entries.Sort();
                M.A.Append_Elements(entries);}}

        for(int i=0;i<M.offsets.m-1;i++) M.offsets(i+1)+=M.offsets(i);}
}
template class SYSTEM_VOLUME_BLOCK_HELPER<VECTOR<float,2> >;
template class SYSTEM_VOLUME_BLOCK_HELPER<VECTOR<float,3> >;
template void SYSTEM_VOLUME_BLOCK_HELPER<VECTOR<float,2> >::Initialize<0,1>(BASIS_STENCIL_UNIFORM<VECTOR<float,2>,0> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,CELL_MANAGER<VECTOR<float,2> >&,CELL_MANAGER<VECTOR<float,2> >&,
    CELL_DOMAIN_INTERFACE<VECTOR<float,2> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER<VECTOR<float,2> >::Initialize<1,1>(BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,CELL_MANAGER<VECTOR<float,2> >&,CELL_MANAGER<VECTOR<float,2> >&,
    CELL_DOMAIN_INTERFACE<VECTOR<float,2> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER<VECTOR<float,3> >::Initialize<0,1>(BASIS_STENCIL_UNIFORM<VECTOR<float,3>,0> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,CELL_MANAGER<VECTOR<float,3> >&,CELL_MANAGER<VECTOR<float,3> >&,
    CELL_DOMAIN_INTERFACE<VECTOR<float,3> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER<VECTOR<float,3> >::Initialize<1,1>(BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,CELL_MANAGER<VECTOR<float,3> >&,CELL_MANAGER<VECTOR<float,3> >&,
    CELL_DOMAIN_INTERFACE<VECTOR<float,3> >&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SYSTEM_VOLUME_BLOCK_HELPER<VECTOR<double,2> >;
template class SYSTEM_VOLUME_BLOCK_HELPER<VECTOR<double,3> >;
template void SYSTEM_VOLUME_BLOCK_HELPER<VECTOR<double,2> >::Initialize<0,1>(BASIS_STENCIL_UNIFORM<VECTOR<double,2>,0> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,CELL_MANAGER<VECTOR<double,2> >&,CELL_MANAGER<VECTOR<double,2> >&,
    CELL_DOMAIN_INTERFACE<VECTOR<double,2> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER<VECTOR<double,2> >::Initialize<1,1>(BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,CELL_MANAGER<VECTOR<double,2> >&,CELL_MANAGER<VECTOR<double,2> >&,
    CELL_DOMAIN_INTERFACE<VECTOR<double,2> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER<VECTOR<double,3> >::Initialize<0,1>(BASIS_STENCIL_UNIFORM<VECTOR<double,3>,0> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,CELL_MANAGER<VECTOR<double,3> >&,CELL_MANAGER<VECTOR<double,3> >&,
    CELL_DOMAIN_INTERFACE<VECTOR<double,3> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER<VECTOR<double,3> >::Initialize<1,1>(BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,CELL_MANAGER<VECTOR<double,3> >&,CELL_MANAGER<VECTOR<double,3> >&,
    CELL_DOMAIN_INTERFACE<VECTOR<double,3> >&);
#endif

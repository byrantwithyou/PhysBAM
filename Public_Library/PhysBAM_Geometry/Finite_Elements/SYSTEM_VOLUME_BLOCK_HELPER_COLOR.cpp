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
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK_HELPER_COLOR.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> template<int d0,int d1> void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>::
Initialize(const BASIS_STENCIL_UNIFORM<TV,d0>& s0,const BASIS_STENCIL_UNIFORM<TV,d1>& s1,
    CELL_MANAGER_COLOR<TV>& cm0_input,CELL_MANAGER_COLOR<TV>& cm1_input,CELL_DOMAIN_INTERFACE_COLOR<TV> &cdi_input)
{
    cm0=&cm0_input;
    cm1=&cm1_input;
    cdi=&cdi_input;

    for(int i=0;i<s0.diced.m;i++)
        for(int j=0;j<s1.diced.m;j++)
            if(s0.diced(i).subcell&s1.diced(j).subcell)
                flat_diff.Append(cdi->Flatten_Diff(s1.diced(j).index_offset-s0.diced(i).index_offset));

    flat_diff.Prune_Duplicates();
    flat_diff.Sort();
    
    data.Resize(cdi->colors);
    for(int c=0;c<cdi->colors;c++) data(c).Resize(cdi->flat_size,flat_diff.m);
}
//#####################################################################
// Function Mark_Active_Cells
//#####################################################################
template<class TV> void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>::
Mark_Active_Cells(T tol)
{
    if(cdi->wrap) // add ghost rows to material ones 
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(cdi->grid,cdi->padding,GRID<TV>::GHOST_REGION,-1);it.Valid();it.Next()){
            int i=cdi->Flatten(it.index);
            int r=cdi->remap(i);
            for(int c=0;c<cdi->colors;c++)
                for(int k=0;k<data(c).n;k++){
                    data(c)(r,k)+=data(c)(i,k);
                    data(c)(i,k)=0;}}
    for(int c=0;c<cdi->colors;c++)
        for(int l=0;l<data(c).m;l++)
            for(int k=0;k<data(c).n;k++)
                if(abs(data(c)(l,k))>tol){
                    cm0->Set_Active(l,c);
                    cm1->Set_Active(l+flat_diff(k),c);}
                else data(c)(l,k)=0;
}
//#####################################################################
// Function Build_Matrix
//#####################################################################
template<class TV> void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>::
Build_Matrix(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& matrix)
{
    if(!cdi->wrap) PHYSBAM_FATAL_ERROR();
    matrix.Resize(cdi->colors);

    for(int c=0;c<cdi->colors;c++){
        SPARSE_MATRIX_FLAT_MXN<T>& M=matrix(c);
        MATRIX_MXN<T>& d=data(c);
        ARRAY<int>& comp_m=cm0->compressed(c);
        ARRAY<int>& comp_n=cm1->compressed(c);
        int m=cm0->dofs(c);
        int n=cm1->dofs(c);
            
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
namespace PhysBAM{
template class SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,2> >;
template class SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,3> >;
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,2> >::Initialize<0,1>(BASIS_STENCIL_UNIFORM<VECTOR<float,2>,0> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,CELL_MANAGER_COLOR<VECTOR<float,2> >&,CELL_MANAGER_COLOR<VECTOR<float,2> >&,
    CELL_DOMAIN_INTERFACE_COLOR<VECTOR<float,2> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,2> >::Initialize<1,1>(BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,CELL_MANAGER_COLOR<VECTOR<float,2> >&,CELL_MANAGER_COLOR<VECTOR<float,2> >&,
    CELL_DOMAIN_INTERFACE_COLOR<VECTOR<float,2> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,3> >::Initialize<0,1>(BASIS_STENCIL_UNIFORM<VECTOR<float,3>,0> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,CELL_MANAGER_COLOR<VECTOR<float,3> >&,CELL_MANAGER_COLOR<VECTOR<float,3> >&,
    CELL_DOMAIN_INTERFACE_COLOR<VECTOR<float,3> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,3> >::Initialize<1,1>(BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,CELL_MANAGER_COLOR<VECTOR<float,3> >&,CELL_MANAGER_COLOR<VECTOR<float,3> >&,
    CELL_DOMAIN_INTERFACE_COLOR<VECTOR<float,3> >&);
template class SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,2> >;
template class SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,3> >;
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,2> >::Initialize<0,1>(BASIS_STENCIL_UNIFORM<VECTOR<double,2>,0> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,CELL_MANAGER_COLOR<VECTOR<double,2> >&,CELL_MANAGER_COLOR<VECTOR<double,2> >&,
    CELL_DOMAIN_INTERFACE_COLOR<VECTOR<double,2> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,2> >::Initialize<1,1>(BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,CELL_MANAGER_COLOR<VECTOR<double,2> >&,CELL_MANAGER_COLOR<VECTOR<double,2> >&,
    CELL_DOMAIN_INTERFACE_COLOR<VECTOR<double,2> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,3> >::Initialize<0,1>(BASIS_STENCIL_UNIFORM<VECTOR<double,3>,0> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,CELL_MANAGER_COLOR<VECTOR<double,3> >&,CELL_MANAGER_COLOR<VECTOR<double,3> >&,
    CELL_DOMAIN_INTERFACE_COLOR<VECTOR<double,3> >&);
template void SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,3> >::Initialize<1,1>(BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,CELL_MANAGER_COLOR<VECTOR<double,3> >&,CELL_MANAGER_COLOR<VECTOR<double,3> >&,
    CELL_DOMAIN_INTERFACE_COLOR<VECTOR<double,3> >&);
}

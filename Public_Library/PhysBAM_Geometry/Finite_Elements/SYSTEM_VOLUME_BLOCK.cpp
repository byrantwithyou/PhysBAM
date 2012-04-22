//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Classes SYSTEM_VOLUME_BLOCK
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV,int static_degree> template<int d0,int d1> void SYSTEM_VOLUME_BLOCK<TV,static_degree>::
Initialize(const BASIS_STENCIL_UNIFORM<TV,d0>& s0,const BASIS_STENCIL_UNIFORM<TV,d1>& s1,CELL_MANAGER<TV>& cm0_input,CELL_MANAGER<TV>& cm1_input,
    CELL_DOMAIN_INTERFACE<TV> &cdi_input,const VECTOR<T,2>& scale_input)
{
    scale=scale_input;
    cm0=&cm0_input;
    cm1=&cm1_input;
    cdi=&cdi_input;
    RANGE<TV_INT> overlap_range(TV_INT(),TV_INT()+1);
    HASHTABLE<int,int> ht;
    for(int i=0;i<s0.diced.m;i++)
        for(int j=0;j<s1.diced.m;j++){
            int overlap=s0.diced(i).subcell&s1.diced(j).subcell;
            if(overlap){
                const typename BASIS_STENCIL_UNIFORM<TV,d0>::DICED& diced0=s0.diced(i);
                const typename BASIS_STENCIL_UNIFORM<TV,d1>::DICED& diced1=s1.diced(j);
                OVERLAP_POLYNOMIAL op={diced0.index_offset,diced1.index_offset,overlap};
                op.polynomial=diced0.polynomial*diced1.polynomial;
                int fd=cdi->Flatten_Diff(diced1.index_offset-diced0.index_offset);
                if(!ht.Get(fd,op.flat_diff_index)){
                    ht.Insert(fd,flat_diff.m);
                    op.flat_diff_index=flat_diff.m;
                    flat_diff.Append(fd);}
                overlap_polynomials.Append(op);}}
    for(int s=0;s<2;s++) data[s].Resize(cdi->flat_size,flat_diff.m);
}
//#####################################################################
// Function Add_Open_Entries
//#####################################################################
template<class TV,int static_degree> void SYSTEM_VOLUME_BLOCK<TV,static_degree>::
Add_Open_Entries(const TV_INT& cell,int inside)
{
    for(int j=0;j<open_entries.m;j++) Add_Open_Entry(cell,inside,open_entries(j));
}
//#####################################################################
// Function Add_Open_Subcell_Entries
//#####################################################################
template<class TV,int static_degree> void SYSTEM_VOLUME_BLOCK<TV,static_degree>::
Add_Open_Subcell_Entries(const TV_INT& cell,int block,int inside)
{
    for(int j=0;j<open_entries.m;j++) Add_Open_Entry(cell,block,open_subcell_entries[block](j));
}
template class SYSTEM_VOLUME_BLOCK<VECTOR<float,2>,2>;
template class SYSTEM_VOLUME_BLOCK<VECTOR<float,3>,2>;
template void SYSTEM_VOLUME_BLOCK<VECTOR<float,3>,2>::Initialize<1,0>(
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,0> const&,
    CELL_MANAGER<VECTOR<float,3> >&,CELL_MANAGER<VECTOR<float,3> >&,
    CELL_DOMAIN_INTERFACE<VECTOR<float,3> >&,const VECTOR<float,2>&);
template void SYSTEM_VOLUME_BLOCK<VECTOR<float,3>,2>::Initialize<1,1>(
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,
    CELL_MANAGER<VECTOR<float,3> >&,CELL_MANAGER<VECTOR<float,3> >&,
    CELL_DOMAIN_INTERFACE<VECTOR<float,3> >&,const VECTOR<float,2>&);
template void SYSTEM_VOLUME_BLOCK<VECTOR<float,2>,2>::Initialize<1,0>(
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,0> const&,
    CELL_MANAGER<VECTOR<float,2> >&,CELL_MANAGER<VECTOR<float,2> >&,
    CELL_DOMAIN_INTERFACE<VECTOR<float,2> >&,const VECTOR<float,2>&);
template void SYSTEM_VOLUME_BLOCK<VECTOR<float,2>,2>::Initialize<1,1>(
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,
    CELL_MANAGER<VECTOR<float,2> >&,CELL_MANAGER<VECTOR<float,2> >&,
    CELL_DOMAIN_INTERFACE<VECTOR<float,2> >&,const VECTOR<float,2>&);
#ifndef COMPILATE_WITHOUT_DOUBLE_SUPPORT
template class SYSTEM_VOLUME_BLOCK<VECTOR<double,2>,2>;
template class SYSTEM_VOLUME_BLOCK<VECTOR<double,3>,2>;
template void SYSTEM_VOLUME_BLOCK<VECTOR<double,3>,2>::Initialize<1,0>(
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,0> const&,
    CELL_MANAGER<VECTOR<double,3> >&,CELL_MANAGER<VECTOR<double,3> >&,
    CELL_DOMAIN_INTERFACE<VECTOR<double,3> >&,const VECTOR<double,2>&);
template void SYSTEM_VOLUME_BLOCK<VECTOR<double,3>,2>::Initialize<1,1>(
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,
    CELL_MANAGER<VECTOR<double,3> >&,CELL_MANAGER<VECTOR<double,3> >&,
    CELL_DOMAIN_INTERFACE<VECTOR<double,3> >&,const VECTOR<double,2>&);
template void SYSTEM_VOLUME_BLOCK<VECTOR<double,2>,2>::Initialize<1,0>(
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,0> const&,
    CELL_MANAGER<VECTOR<double,2> >&,CELL_MANAGER<VECTOR<double,2> >&,
    CELL_DOMAIN_INTERFACE<VECTOR<double,2> >&,const VECTOR<double,2>&);
template void SYSTEM_VOLUME_BLOCK<VECTOR<double,2>,2>::Initialize<1,1>(
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,
    CELL_MANAGER<VECTOR<double,2> >&,CELL_MANAGER<VECTOR<double,2> >&,
    CELL_DOMAIN_INTERFACE<VECTOR<double,2> >&,const VECTOR<double,2>&);
#endif

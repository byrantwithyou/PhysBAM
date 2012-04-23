//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Classes SYSTEM_INTERFACE_BLOCK
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_INTERFACE_BLOCK.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV,int static_degree> template<int d> void SYSTEM_INTERFACE_BLOCK<TV,static_degree>::
Initialize(const BASIS_STENCIL_UNIFORM<TV,d>& s,CELL_MANAGER<TV>& cm_input,CELL_DOMAIN_INTERFACE<TV>& cdi_input,T scale_input,bool ignore_orientation_input)
{
    scale=scale_input;
    ignore_orientation=ignore_orientation_input;
    cm=&cm_input;
    cdi=&cdi_input;

    ARRAY<int> flat_diffs;
    overlap_polynomials.Resize(s.diced.m);
    for(int i=0;i<overlap_polynomials.m;i++){
        OVERLAP_POLYNOMIAL& op=overlap_polynomials(i);
        const typename BASIS_STENCIL_UNIFORM<TV,d>::DICED& diced=s.diced(i);
        op.index_offset=diced.index_offset;
        op.subcell=diced.subcell;
        op.polynomial=diced.polynomial;
        op.flat_diff_index.Resize(cdi->coarse_range);
        for(RANGE_ITERATOR<TV::m> it(cdi->coarse_range);it.Valid();it.Next())
            flat_diffs.Append(cdi->Flatten_Diff(it.index+diced.index_offset));}

    Sort(flat_diffs);
    flat_diffs.Prune_Duplicates();

    for(int i=0;i<overlap_polynomials.m;i++){
        OVERLAP_POLYNOMIAL& op=overlap_polynomials(i);
        for(RANGE_ITERATOR<TV::m> it(cdi->coarse_range);it.Valid();it.Next())
            op.flat_diff_index(it.index)=flat_diffs.Binary_Search(cdi->Flatten_Diff(it.index+op.index_offset));}
    
    for(int s=0;s<2;s++) data[s].Resize(cdi->interface_elements,flat_diff.m);
}
template class SYSTEM_INTERFACE_BLOCK<VECTOR<float,2>,2>;
template class SYSTEM_INTERFACE_BLOCK<VECTOR<float,3>,2>;
template void SYSTEM_INTERFACE_BLOCK<VECTOR<float,3>,2>::Initialize<1>(
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,CELL_MANAGER<VECTOR<float,3> >&,
    CELL_DOMAIN_INTERFACE<VECTOR<float,3> >&,float,bool);
template void SYSTEM_INTERFACE_BLOCK<VECTOR<float,2>,2>::Initialize<1>(
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,CELL_MANAGER<VECTOR<float,2> >&,
    CELL_DOMAIN_INTERFACE<VECTOR<float,2> >&,float,bool);
#ifndef COMPILATE_WITHOUT_DOUBLE_SUPPORT
template class SYSTEM_INTERFACE_BLOCK<VECTOR<double,2>,2>;
template class SYSTEM_INTERFACE_BLOCK<VECTOR<double,3>,2>;
template void SYSTEM_INTERFACE_BLOCK<VECTOR<double,3>,2>::Initialize<1>(
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,CELL_MANAGER<VECTOR<double,3> >&,
    CELL_DOMAIN_INTERFACE<VECTOR<double,3> >&,double,bool);
template void SYSTEM_INTERFACE_BLOCK<VECTOR<double,2>,2>::Initialize<1>(
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,CELL_MANAGER<VECTOR<double,2> >&,
    CELL_DOMAIN_INTERFACE<VECTOR<double,2> >&,double,bool);
#endif

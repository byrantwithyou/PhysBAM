//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK_HELPER_COLOR.h>
using namespace PhysBAM;
template<class TV,int static_degree> template<int d0,int d1> void SYSTEM_VOLUME_BLOCK_COLOR<TV,static_degree>::
Initialize(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>& helper_input,const BASIS_STENCIL_UNIFORM<TV,d0>& s0,
    const BASIS_STENCIL_UNIFORM<TV,d1>& s1,const ARRAY<T>& scale_input)
{
    scale=scale_input;
    helper=&helper_input;

    for(int i=0;i<s0.diced.m;i++)
        for(int j=0;j<s1.diced.m;j++){
            int overlap=s0.diced(i).subcell&s1.diced(j).subcell;
            if(overlap){
                const typename BASIS_STENCIL_UNIFORM<TV,d0>::DICED& diced0=s0.diced(i);
                const typename BASIS_STENCIL_UNIFORM<TV,d1>::DICED& diced1=s1.diced(j);
                OVERLAP_POLYNOMIAL op;
                op.flat_index_offset=helper->cdi->Flatten_Diff(diced0.index_offset);
                op.flat_index_diff_ref=helper->flat_diff.Binary_Search(helper->cdi->Flatten_Diff(diced1.index_offset-diced0.index_offset));
                op.subcell=overlap;
                op.polynomial=diced0.polynomial*diced1.polynomial;
                overlap_polynomials.Append(op);}}
}
template<class TV,int static_degree> void SYSTEM_VOLUME_BLOCK_COLOR<TV,static_degree>::
Add_Open_Entries(int flat_index,int color)
{
    for(int j=0;j<open_entries.m;j++) Add_Open_Entry(flat_index,color,open_entries(j));
}
template<class TV,int static_degree> void SYSTEM_VOLUME_BLOCK_COLOR<TV,static_degree>::
Add_Open_Subcell_Entries(int flat_index,int subcell,int color)
{
    for(int j=0;j<open_subcell_entries[subcell].m;j++) Add_Open_Entry(flat_index,color,open_subcell_entries[subcell](j));
}
template class SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<float,2>,2>;
template class SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<float,3>,2>;
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<float,2>,2>::Initialize<0,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,ARRAY<float> const&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<float,2>,2>::Initialize<1,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,ARRAY<float> const&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<float,3>,2>::Initialize<0,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,ARRAY<float> const&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<float,3>,2>::Initialize<1,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,ARRAY<float> const&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<double,2>,2>;
template class SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<double,3>,2>;
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<double,2>,2>::Initialize<0,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,ARRAY<double> const&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<double,2>,2>::Initialize<1,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,ARRAY<double> const&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<double,3>,2>::Initialize<0,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,ARRAY<double> const&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<double,3>,2>::Initialize<1,1>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,ARRAY<double> const&);
#endif

//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Utilities/NONCOPYABLE.h>
#include <Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK_COLOR.h>
#include <Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK_HELPER_COLOR.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV,int static_degree> template<int d0,class ...Args> void SYSTEM_VOLUME_BLOCK_COLOR<TV,static_degree>::
Initialize(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>& helper_input,const ARRAY<T>& scale_input,const BASIS_STENCIL_UNIFORM<TV,d0>& s0,Args&& ...args)
{
    scale=scale_input;
    helper=&helper_input;

    ARRAY<int> diffs;
    for(int i=0;i<s0.diced.m;i++){
        const typename BASIS_STENCIL_UNIFORM<TV,d0>::DICED& diced0=s0.diced(i);
        Initialize_Helper(diced0.subcell,helper->cdi->Flatten_Diff(diced0.index_offset),diced0.index_offset,diced0.polynomial,diffs,args...);}
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV,int static_degree> template<int d1,int pd,class ...Args> void SYSTEM_VOLUME_BLOCK_COLOR<TV,static_degree>::
Initialize_Helper(int subcell,int flat_index_offset,const TV_INT& index_offset,const STATIC_POLYNOMIAL<T,TV::m,pd>& poly,ARRAY<int>& diffs,const BASIS_STENCIL_UNIFORM<TV,d1>& s1,Args&& ...args)
{
    for(int j=0;j<s1.diced.m;j++)
        if(int overlap=subcell&s1.diced(j).subcell){
            const typename BASIS_STENCIL_UNIFORM<TV,d1>::DICED& diced1=s1.diced(j);
            diffs.Append(helper->cdi->Flatten_Diff(diced1.index_offset-index_offset));
            Initialize_Helper(overlap,flat_index_offset,index_offset,poly*diced1.polynomial,diffs,args...);
            diffs.Pop();}
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV,int static_degree> template<int pd> void SYSTEM_VOLUME_BLOCK_COLOR<TV,static_degree>::
Initialize_Helper(int subcell,int flat_index_offset,const TV_INT& index_offset,const STATIC_POLYNOMIAL<T,TV::m,pd>& poly,ARRAY<int>& diffs)
{
    OVERLAP_POLYNOMIAL op;
    op.flat_index_offset=flat_index_offset;
    op.flat_index_diff_ref=helper->flat_diff.Binary_Search(diffs,LEXICOGRAPHIC_COMPARE());
    op.subcell=subcell;
    op.polynomial=poly;
    overlap_polynomials.Append(op);
}
//#####################################################################
// Function Add_Open_Entries
//#####################################################################
template<class TV,int static_degree> void SYSTEM_VOLUME_BLOCK_COLOR<TV,static_degree>::
Add_Open_Entries(int flat_index,int color)
{
    for(int j=0;j<open_entries.m;j++) Add_Open_Entry(flat_index,color,open_entries(j));
}
//#####################################################################
// Function Add_Open_Subcell_Entries
//#####################################################################
template<class TV,int static_degree> void SYSTEM_VOLUME_BLOCK_COLOR<TV,static_degree>::
Add_Open_Subcell_Entries(int flat_index,int subcell,int color)
{
    for(int j=0;j<open_subcell_entries[subcell].m;j++) Add_Open_Entry(flat_index,color,open_subcell_entries[subcell](j));
}
//#####################################################################
// Function Add_Entry
//#####################################################################
template<class TV,int static_degree> void SYSTEM_VOLUME_BLOCK_COLOR<TV,static_degree>::
Add_Entry(int flat_index,int flat_index_diff_ref,int color,T value)
{
    helper->data(color)(flat_index,flat_index_diff_ref)+=value*scale(color);
}
//#####################################################################
// Function Add_Open_Entry
//#####################################################################
template<class TV,int static_degree> void SYSTEM_VOLUME_BLOCK_COLOR<TV,static_degree>::
Add_Open_Entry(int flat_index,int color,OPEN_ENTRY& oe)
{
    Add_Entry(flat_index+oe.flat_index_offset,oe.flat_index_diff_ref,color,oe.x);
}
namespace PhysBAM{
template class SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<float,2>,2>;
template class SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<float,3>,2>;
template class SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<double,2>,2>;
template class SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<double,3>,2>;
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<double,2>,2>::Initialize<0,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&>
    (SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,2> >&,ARRAY<double,int> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<double,2>,2>::Initialize<1,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&>
    (SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,2> >&,ARRAY<double,int> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<double,3>,2>::Initialize<0,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&>
    (SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,3> >&,ARRAY<double,int> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<double,3>,2>::Initialize<1,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&>
    (SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,3> >&,ARRAY<double,int> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<float,2>,2>::Initialize<0,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&>
    (SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,2> >&,ARRAY<float,int> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<float,2>,2>::Initialize<1,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&>
    (SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,2> >&,ARRAY<float,int> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<float,3>,2>::Initialize<0,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&>
    (SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,3> >&,ARRAY<float,int> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<float,3>,2>::Initialize<1,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&>
    (SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,3> >&,ARRAY<float,int> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<double,2>,2>::Initialize<0,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1>&>
    (SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,2> >&,ARRAY<double,int> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1>&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<double,2>,2>::Initialize<1,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1>&>
    (SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,2> >&,ARRAY<double,int> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1>&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<double,3>,2>::Initialize<0,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1>&>
    (SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,3> >&,ARRAY<double,int> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1>&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<double,3>,2>::Initialize<1,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1>&>
    (SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,3> >&,ARRAY<double,int> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1>&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<float,2>,2>::Initialize<0,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1>&>
    (SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,2> >&,ARRAY<float,int> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1>&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<float,2>,2>::Initialize<1,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1>&>
    (SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,2> >&,ARRAY<float,int> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1>&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<float,3>,2>::Initialize<0,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1>&>
    (SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,3> >&,ARRAY<float,int> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,0> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1>&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<float,3>,2>::Initialize<1,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1>&>
    (SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,3> >&,ARRAY<float,int> const&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1>&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<double,2>,2>::Initialize<1,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1>&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1>&>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,2> >&,
    ARRAY<double,int> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1>&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1>&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<double,3>,2>::Initialize<1,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1>&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1>&>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<double,3> >&,
    ARRAY<double,int> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1>&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1>&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<float,2>,2>::Initialize<1,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1>&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1>&>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,2> >&,
    ARRAY<float,int> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1>&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1>&);
template void SYSTEM_VOLUME_BLOCK_COLOR<VECTOR<float,3>,2>::Initialize<1,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1>&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1>&>(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<VECTOR<float,3> >&,
    ARRAY<float,int> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1>&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1>&);
}

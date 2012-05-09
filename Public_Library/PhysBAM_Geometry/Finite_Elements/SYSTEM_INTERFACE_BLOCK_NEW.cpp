//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_NEW.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_INTERFACE_BLOCK_NEW.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_INTERFACE_BLOCK_HELPER_NEW.h>
using namespace PhysBAM;
template<class TV,int static_degree> template<int d> void SYSTEM_INTERFACE_BLOCK_NEW<TV,static_degree>::
Initialize(SYSTEM_INTERFACE_BLOCK_HELPER_NEW<TV>& helper_input,const BASIS_STENCIL_UNIFORM<TV,d>& s,T scale_input,bool ignore_orientation_input)
{
    scale=scale_input;
    ignore_orientation=ignore_orientation_input;
    helper=&helper_input;
    cdi=helper->cdi;

    overlap_polynomials.Resize(s.diced.m);
    for(int i=0;i<overlap_polynomials.m;i++){
        OVERLAP_POLYNOMIAL& op=overlap_polynomials(i);
        const typename BASIS_STENCIL_UNIFORM<TV,d>::DICED& diced=s.diced(i);
        op.flat_index_offset=cdi->Flatten_Diff(diced.index_offset);
        op.flat_index_diff_ref=helper->flat_diff.Binary_Search(op.flat_index_offset);
        op.subcell=diced.subcell;
        op.polynomial=diced.polynomial;}
}
template void SYSTEM_INTERFACE_BLOCK_NEW<VECTOR<float,2>,2>::Initialize<1>(SYSTEM_INTERFACE_BLOCK_HELPER_NEW<VECTOR<float,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,float,bool);
template void SYSTEM_INTERFACE_BLOCK_NEW<VECTOR<float,3>,2>::Initialize<1>(SYSTEM_INTERFACE_BLOCK_HELPER_NEW<VECTOR<float,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,float,bool);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void SYSTEM_INTERFACE_BLOCK_NEW<VECTOR<double,2>,2>::Initialize<1>(SYSTEM_INTERFACE_BLOCK_HELPER_NEW<VECTOR<double,2> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,double,bool);
template void SYSTEM_INTERFACE_BLOCK_NEW<VECTOR<double,3>,2>::Initialize<1>(SYSTEM_INTERFACE_BLOCK_HELPER_NEW<VECTOR<double,3> >&,BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,double,bool);
#endif

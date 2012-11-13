//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_SURFACE_BLOCK_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_SURFACE_BLOCK_HELPER_COLOR.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV,int static_degree> template<int d> void SYSTEM_SURFACE_BLOCK_COLOR<TV,static_degree>::
Initialize(SYSTEM_SURFACE_BLOCK_HELPER_COLOR<TV>& helper_input,const BASIS_STENCIL_UNIFORM<TV,d>& s,BOUNDARY_CONDITIONS_COLOR<TV>* bc_input,
    ARRAY<ARRAY<T> >& rhs_input,int axis_input,T scale_input)
{
    bc=bc_input;
    axis=axis_input;
    scale=scale_input;
    helper=&helper_input;
    rhs=&rhs_input;
    
    overlap_polynomials.Resize(s.diced.m);
    for(int i=0;i<overlap_polynomials.m;i++){
        OVERLAP_POLYNOMIAL& op=overlap_polynomials(i);
        const typename BASIS_STENCIL_UNIFORM<TV,d>::DICED& diced=s.diced(i);
        op.flat_index_offset=helper->cdi->Flatten_Diff(diced.index_offset);
        op.flat_index_diff_ref=helper->flat_diff.Binary_Search(op.flat_index_offset);
        op.subcell=diced.subcell;
        op.polynomial=diced.polynomial;}
}
template void SYSTEM_SURFACE_BLOCK_COLOR<VECTOR<float,2>,2>::Initialize<1>(SYSTEM_SURFACE_BLOCK_HELPER_COLOR<VECTOR<float,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,2>,1> const&,BOUNDARY_CONDITIONS_COLOR<VECTOR<float,2> >*,ARRAY<ARRAY<float> >&,int,float);
template void SYSTEM_SURFACE_BLOCK_COLOR<VECTOR<float,3>,2>::Initialize<1>(SYSTEM_SURFACE_BLOCK_HELPER_COLOR<VECTOR<float,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<float,3>,1> const&,BOUNDARY_CONDITIONS_COLOR<VECTOR<float,3> >*,ARRAY<ARRAY<float> >&,int,float);
template void SYSTEM_SURFACE_BLOCK_COLOR<VECTOR<double,2>,2>::Initialize<1>(SYSTEM_SURFACE_BLOCK_HELPER_COLOR<VECTOR<double,2> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,2>,1> const&,BOUNDARY_CONDITIONS_COLOR<VECTOR<double,2> >*,ARRAY<ARRAY<double> >&,int,double);
template void SYSTEM_SURFACE_BLOCK_COLOR<VECTOR<double,3>,2>::Initialize<1>(SYSTEM_SURFACE_BLOCK_HELPER_COLOR<VECTOR<double,3> >&,
    BASIS_STENCIL_UNIFORM<VECTOR<double,3>,1> const&,BOUNDARY_CONDITIONS_COLOR<VECTOR<double,3> >*,ARRAY<ARRAY<double> >&,int,double);

//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/UNIFORM_ARRAY_ITERATOR.h>
#include <PhysBAM_Tools/Symbolics/MULTIVARIATE_POLYNOMIAL.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MAPPING.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BASIS_INTEGRATION_UNIFORM<TV>::
BASIS_INTEGRATION_UNIFORM(const GRID<TV>& grid_input,const GRID<TV>& phi_grid_input)
    :grid(grid_input),phi_grid(phi_grid_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BASIS_INTEGRATION_UNIFORM<TV>::
~BASIS_INTEGRATION_UNIFORM()
{
}
//#####################################################################
// Function Compute_Matrix
//#####################################################################
template<class TV> void BASIS_INTEGRATION_UNIFORM<TV>::
Compute_Matrix(SYSTEM_MATRIX_HELPER<T>& helper,const BASIS_STENCIL_UNIFORM<TV>& s0, const BASIS_STENCIL_UNIFORM<TV>& s1, CELL_MAPPING<TV>& cm0, CELL_MAPPING<TV>& cm1)
{
    ARRAY<OVERLAP_POLYNOMIALS> overlap_polynomials;
    for(int i=0;i<s0.diced.m;i++)
        for(int j=0;j<s1.diced.m;j++){
            RANGE<TV_INT> overlap=RANGE<TV_INT>::Intersect(s0.diced(i).range,s1.diced(j).range);
            if(!overlap.Empty_Half_Open()){
                OVERLAP_POLYNOMIALS op = {s0.diced(i).index_offset, s1.diced(j).index_offset, overlap};
                op.polynomial=s0.diced(i).polynomial;
                op.polynomial*=s1.diced(j).polynomial;
                overlap_polynomials.Append(op);}}

    ARRAY<MATRIX_ENTRY> open_entries; // stencil with no boundary conditions
    for(int i=0;i<overlap_polynomials.m;i++){
        RANGE<TV> box=RANGE<TV>(overlap_polynomials(i).range)/2;
        MATRIX_ENTRY me = {overlap_polynomials(i).index_offset0, overlap_polynomials(i).index_offset1, overlap_polynomials(i).polynomial.Definite_Integral(box)};
        open_entries.Append(me);}
    open_entries.Coalesce();

    helper.New_Block();
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        for(int i=0;i<open_entries.m;i++){
            MATRIX_ENTRY me=open_entries(i);
            me.index0+=it.index;
            me.index1+=it.index;

            int i0=cm0.Get_Index(me.index0,0);
            int i1=cm1.Get_Index(me.index1,0);
            assert(i0>=0 && i1>=0);

            helper.data.Append(TRIPLE<int,int,T>(i0,i1,me.x));}}
}
template class BASIS_INTEGRATION_UNIFORM<VECTOR<float,1> >;
template class BASIS_INTEGRATION_UNIFORM<VECTOR<float,2> >;
template class BASIS_INTEGRATION_UNIFORM<VECTOR<float,3> >;
#ifndef COMPILATE_WITHOUT_DOUBLE_SUPPORT
template class BASIS_INTEGRATION_UNIFORM<VECTOR<double,1> >;
template class BASIS_INTEGRATION_UNIFORM<VECTOR<double,2> >;
template class BASIS_INTEGRATION_UNIFORM<VECTOR<double,3> >;
#endif

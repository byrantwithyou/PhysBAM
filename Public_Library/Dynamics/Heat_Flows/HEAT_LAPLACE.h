//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HEAT_LAPLACE
//#####################################################################
//
// Helper class to convert a Laplace-type matrix of the form
//   -Eu = -u^* + b
// where E is an elliptic operator, to a heat equation matrix of the form
//   -1 - cEu = -u^* + cb
// Here the -u^* term is due to f, and b is due to Dirichlet boundary conditions.
//
//#####################################################################
#ifndef __HEAT_LAPLACE__
#define __HEAT_LAPLACE__

#include <Core/Arrays/ARRAY.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Math_Tools/RANGE.h>
#include <Grid_Tools/Arrays/ARRAYS_UNIFORM_FORWARD.h>
namespace PhysBAM{

template<class TV> struct GRID_ARRAYS_POLICY;
template<class TV> class GRID;
template<class T> class SPARSE_MATRIX_FLAT_MXN;

template<class T_LAPLACE>
class HEAT_LAPLACE:public T_LAPLACE
{
    typedef typename T_LAPLACE::VECTOR_T TV;typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef T_LAPLACE BASE;
    using BASE::filled_region_touches_dirichlet;using BASE::Solve;using BASE::Find_Solution_Regions;using BASE::f;

    T coefficient;

    HEAT_LAPLACE(GRID<TV>& grid_input,ARRAY<T,TV_INT>& u_input)
        :T_LAPLACE(grid_input,u_input,true,false,false),coefficient(1)
    {}

    HEAT_LAPLACE(GRID<TV>& grid_input,ARRAY<T,TV_INT>& u_input,const bool initialize_grid_input,const bool multiphase_input)
        :T_LAPLACE(grid_input,u_input,initialize_grid_input,multiphase_input,false),coefficient(1)
    {}

    void Solve(const T time,const bool solution_regions_already_computed=false)
    {if(!solution_regions_already_computed) Find_Solution_Regions();
    filled_region_touches_dirichlet.Fill(true); // don't need to worry about Neumann regions in heat flow
    T_LAPLACE::Solve(time,true);}

    void Find_A_Part_Two(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& A_array,ARRAY<ARRAY<T> >& b_array,ARRAY<int,TV_INT>& cell_index_to_matrix_index)
    {if(!coefficient) PHYSBAM_FATAL_ERROR();
    if(coefficient!=1) f/=coefficient; // this will cancel with the multiplication below, leaving only the Dirichlet part multiplied by coefficient
    T_LAPLACE::Find_A_Part_Two(A_array,b_array,cell_index_to_matrix_index);
    A_array*=(T)3;
    if(coefficient!=1){A_array*=coefficient;b_array*=coefficient;}
    for(int i=0;i<A_array.m;i++) A_array(i)-=(T)1;}

//#####################################################################
};
}
#endif

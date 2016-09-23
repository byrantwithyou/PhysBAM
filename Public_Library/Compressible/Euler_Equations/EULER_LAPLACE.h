//#####################################################################
// Copyright 2007, Jon Gretarsson, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_LAPLACE
//#####################################################################
#ifndef __EULER_LAPLACE__
#define __EULER_LAPLACE__

#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Grid_Tools/Arrays/ARRAYS_UNIFORM_FORWARD.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
namespace PhysBAM{

template<class TV> struct GRID_ARRAYS_POLICY;

template<class T_LAPLACE>
class EULER_LAPLACE:public T_LAPLACE
{
    typedef typename T_LAPLACE::GRID_T::VECTOR_T TV;typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef TV_INT INDEX;
public:
    typedef T_LAPLACE BASE;
    using BASE::filled_region_touches_dirichlet;using BASE::Solve;using BASE::Find_Solution_Regions;using BASE::f;
    using BASE::grid;using BASE::solve_neumann_regions;using BASE::dt;using BASE::dt_is_set;using BASE::Initialize_Grid;

    const ARRAY<T,TV_INT>& one_over_rho_c_squared;

    EULER_LAPLACE(const GRID<TV>& grid_input,ARRAY<T,TV_INT>& u_input,const ARRAY<T,TV_INT>& one_over_rho_c_squared_input)
        :T_LAPLACE(grid_input,u_input,true,false,false),one_over_rho_c_squared(one_over_rho_c_squared_input)
    {}

    EULER_LAPLACE(const GRID<TV>& grid_input,ARRAY<T,TV_INT>& u_input,const ARRAY<T,TV_INT>& one_over_rho_c_squared_input,const bool flag)
        :T_LAPLACE(grid_input,u_input,flag,flag,false),one_over_rho_c_squared(one_over_rho_c_squared_input)
    {}

    void Solve(const T time,const bool solution_regions_already_computed=false)
    {if(!solution_regions_already_computed) Find_Solution_Regions();
    T_LAPLACE::Solve(time,true);}

    void Find_A(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& A_array,ARRAY<ARRAY<T> >& b_array,const ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,
            ARRAY<int,TV_INT>& cell_index_to_matrix_index)
    {assert(dt_is_set);dt_is_set=false;

    T_LAPLACE::Find_A(domain,A_array,b_array,filled_region_cell_count,cell_index_to_matrix_index);
    for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next()){INDEX cell_index=iterator.Cell_Index();
        int color=this->filled_region_colors(cell_index);
        if(color!=-2 && (filled_region_touches_dirichlet(color)||solve_neumann_regions)){
            int matrix_index=cell_index_to_matrix_index(cell_index);
            SPARSE_MATRIX_FLAT_MXN<T>& A=A_array(this->filled_region_colors(cell_index));
            A(matrix_index,matrix_index)-=one_over_rho_c_squared(cell_index)/(dt*dt);}}
    }

//#####################################################################
};
}
#endif

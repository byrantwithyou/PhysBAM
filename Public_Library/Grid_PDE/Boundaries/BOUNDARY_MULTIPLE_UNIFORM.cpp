//#####################################################################
// Copyright 2008, Jon Gretarsson, Nipun Kwatra
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Grid_PDE/Boundaries/BOUNDARY_MULTIPLE_UNIFORM.h>
namespace PhysBAM{
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class TV,class T2> void BOUNDARY_MULTIPLE_UNIFORM<TV,T2>::
Fill_Ghost_Cells(const GRID<TV>& grid,const T_ARRAYS_DIMENSION_T2& u,T_ARRAYS_DIMENSION_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells) const
{
    T_ARRAYS_DIMENSION_T2::Put(u,u_ghost);
    VECTOR<RANGE<TV_INT>,2*TV::m> regions;Find_Ghost_Regions(grid,regions,number_of_ghost_cells);
    for(int side=0;side<2*TV::m;side++) boundaries[side]->Fill_Single_Ghost_Region(grid,u_ghost,regions(side),side,dt,time,number_of_ghost_cells); 
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class TV,class T2> void BOUNDARY_MULTIPLE_UNIFORM<TV,T2>::
Apply_Boundary_Condition(const GRID<TV>& grid,T_ARRAYS_DIMENSION_T2& u,const T time) const
{
    for(int side=0;side<2*TV::m;side++) boundaries[side]->Apply_Boundary_Condition(grid,u,time);
}
template class BOUNDARY_MULTIPLE_UNIFORM<VECTOR<float,2>,VECTOR<float,4> >;
template class BOUNDARY_MULTIPLE_UNIFORM<VECTOR<double,2>,VECTOR<double,4> >;
}

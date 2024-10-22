//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_REFLECTION_ATTENUATION
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_PDE/Boundaries/BOUNDARY_REFLECTION_ATTENUATION.h>
using namespace PhysBAM;
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class TV,class T2> void BOUNDARY_REFLECTION_ATTENUATION<TV,T2>::
Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells) const
{
    ARRAYS_ND_BASE<T2,TV_INT>::Put(u,u_ghost); // interior
    VECTOR<RANGE<TV_INT>,2*TV::m> regions;Find_Ghost_Regions(grid,regions,number_of_ghost_cells);
    for(int side=0;side<GRID<TV>::number_of_faces_per_cell;side++) Fill_Single_Ghost_Region(grid,u_ghost,regions(side),side,dt,time,number_of_ghost_cells);
}
//#####################################################################
// Function Fill_Single_Ghost_Region
//#####################################################################
template<class TV,class T2> void BOUNDARY_REFLECTION_ATTENUATION<TV,T2>::
Fill_Single_Ghost_Region(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const RANGE<TV_INT>& region,const int side,const T dt,const T time,const int number_of_ghost_cells) const
{
    if(Constant_Extrapolation(side)){
        int axis=side/2,boundary=Boundary(side,region);
        for(CELL_ITERATOR<TV> iterator(grid,region);iterator.Valid();iterator.Next()){
            TV_INT cell_index=iterator.Cell_Index(),boundary_cell=cell_index;boundary_cell[axis]=boundary;
            T2 boundary_value=u_ghost(boundary_cell);
            u_ghost(cell_index)=Attenuate_To_Far_Field_Value(boundary_value,dt);}}
    else{BASE::Fill_Single_Ghost_Region(grid,u_ghost,region,side,dt,time,number_of_ghost_cells);}
}
//#####################################################################
// Function Attenuate_To_Far_Field_Value
//#####################################################################
template<class TV,class T2> T2 BOUNDARY_REFLECTION_ATTENUATION<TV,T2>::
Attenuate_To_Far_Field_Value(const T2 boundary_value,const T dt) const
{
    //TODO: Attenuate only when inflow. Will need velocity info too.
    return boundary_value+linear_attenuation*(fixed_boundary_value-boundary_value);
}
//#####################################################################
namespace PhysBAM{
template class BOUNDARY_REFLECTION_ATTENUATION<VECTOR<float,1>,float>;
template class BOUNDARY_REFLECTION_ATTENUATION<VECTOR<float,2>,float>;
template class BOUNDARY_REFLECTION_ATTENUATION<VECTOR<float,3>,float>;
template class BOUNDARY_REFLECTION_ATTENUATION<VECTOR<double,1>,double>;
template class BOUNDARY_REFLECTION_ATTENUATION<VECTOR<double,2>,double>;
template class BOUNDARY_REFLECTION_ATTENUATION<VECTOR<double,3>,double>;
}

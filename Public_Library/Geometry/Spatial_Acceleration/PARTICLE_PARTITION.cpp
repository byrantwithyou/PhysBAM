//#####################################################################
// Copyright 2003-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Matrices/MATRIX.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Spatial_Acceleration/PARTICLE_PARTITION.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLE_PARTITION<TV>::
PARTICLE_PARTITION(const RANGE<TV>& box,const TV_INT& counts,const GEOMETRY_PARTICLES<TV>& particles,const bool use_radius,const bool is_mac_grid)
    :grid(counts,box,is_mac_grid),partition(grid.Domain_Indices()),use_radius(use_radius)
{
    if(use_radius) radius.Resize(grid.Domain_Indices());
    for(int p=0;p<particles.Size();p++) Add_To_Partition(particles.X(p),p);
}
//#####################################################################
// Function Add_To_Partition
//#####################################################################
template<class TV> void PARTICLE_PARTITION<TV>::
Add_To_Partition(const TV& location,const int particle_id)
    {
        assert(grid.Domain().Lazy_Inside(location));
    TV_INT cell=grid.Clamp_To_Cell(location);
    partition(cell).Append(particle_id);
    if(use_radius) radius(cell)=max(radius(cell),(location-grid.Center(cell)).Magnitude());
}
//#####################################################################
// Function Range
//#####################################################################
template<class TV> auto PARTICLE_PARTITION<TV>::
Range(const RANGE<TV>& box) const -> RANGE<TV_INT>
{
    return grid.Clamp_To_Cell(box);
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV> void PARTICLE_PARTITION<TV>::
Intersection_List(const IMPLICIT_OBJECT<TV>& test_surface,const MATRIX<T,d>& rotation,const TV& translation,ARRAY<TV_INT>& intersection_list,const T contour_value) const
{
    PHYSBAM_ASSERT(use_radius);intersection_list.Remove_All();
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell=iterator.Cell_Index(); 
        if(partition(cell).m && !test_surface.Lazy_Outside_Extended_Levelset(rotation*grid.Center(cell)+translation,radius(cell)+contour_value))
            intersection_list.Append(cell);}
}
//#####################################################################
// Function Proximity_List
//#####################################################################
template<class TV> void PARTICLE_PARTITION<TV>::
Proximity_List(const TV& location,const T proximity,ARRAY<int>& proximity_list)
{
        proximity_list.Remove_All();
        for(CELL_ITERATOR<TV> iterator(grid,grid.Clamp_To_Cell(RANGE<TV>(location).Thickened(proximity),0));iterator.Valid();iterator.Next())
            proximity_list.Append_Unique_Elements(partition(iterator.Cell_Index()));
}
template class PARTICLE_PARTITION<VECTOR<double,1> >;
template class PARTICLE_PARTITION<VECTOR<double,2> >;
template class PARTICLE_PARTITION<VECTOR<double,3> >;
template class PARTICLE_PARTITION<VECTOR<float,1> >;
template class PARTICLE_PARTITION<VECTOR<float,2> >;
template class PARTICLE_PARTITION<VECTOR<float,3> >;
}

//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_COLLISION_BODY_INACCURATE_UNION
//##################################################################### 
#include <Tools/Data_Structures/OPERATION_HASH.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Level_Sets/IMPLICIT_OBJECT_ON_A_RAY.h>
#include <Incompressible/Collisions_And_Interactions/FLUID_COLLISION_BODY_INACCURATE_UNION.h>
using namespace PhysBAM;
//#####################################################################
// Function Constructor
//#####################################################################
template<class TV> FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
FLUID_COLLISION_BODY_INACCURATE_UNION(GRID<TV>& grid_input)
    :collision_bodies(grid_input),contour_value(0),grid(grid_input),levelset(grid_input,phi)
{
    collision_geometries_for_rasterization=&collision_bodies.collision_geometry_collection.bodies;
}
template<class TV> FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
FLUID_COLLISION_BODY_INACCURATE_UNION(GRID<TV>& grid_input,T contour_value_input)
    :contour_value(contour_value_input),grid(grid_input),levelset(grid_input,phi)
{
}    
//#####################################################################
// Function Destructor
//#####################################################################
template<class TV> FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
~FLUID_COLLISION_BODY_INACCURATE_UNION()
{
}
//#####################################################################
// Function Implicit_Geometry_Extended_Value
//#####################################################################
template<class TV> typename TV::SCALAR FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Implicit_Geometry_Extended_Value_Helper(const TV& location,UNIFORM_TAG<TV>) const
{return interpolation.Clamped_To_Array(grid,phi,location);}
template<class TV> typename TV::SCALAR FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Implicit_Geometry_Extended_Value(const TV& location) const
{
    return Implicit_Geometry_Extended_Value_Helper(location,typename GRID<TV>::GRID_TAG());
}
//#####################################################################
// Function Update_Intersection_Acceleration_Structures
//#####################################################################
template<class TV> void FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Update_Intersection_Acceleration_Structures(const bool use_swept_simplex_hierarchy,const int state1,const int state2)
{
    for(COLLISION_GEOMETRY_ID i(0);i<collision_bodies.collision_geometry_collection.bodies.m;i++) if(collision_bodies.Is_Active(i)) collision_bodies.collision_geometry_collection.bodies(i)->Update_Intersection_Acceleration_Structures(use_swept_simplex_hierarchy,state1,state2);
}
//#####################################################################
// Function Restore_State
//#####################################################################
template<class TV> void FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Restore_State(const int state_index)
{
    for(COLLISION_GEOMETRY_ID i(0);i<collision_bodies.collision_geometry_collection.bodies.m;i++) if(collision_bodies.Is_Active(i)) collision_bodies.collision_geometry_collection.bodies(i)->Restore_State(state_index);
}
//#####################################################################
// Function Save_State
//#####################################################################
template<class TV> void FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Save_State(const int state_index,const T time)
{
    for(COLLISION_GEOMETRY_ID i(0);i<collision_bodies.collision_geometry_collection.bodies.m;i++) if(collision_bodies.Is_Active(i)) collision_bodies.collision_geometry_collection.bodies(i)->Save_State(state_index,time);
}
//#####################################################################
// Function Read_State
//#####################################################################
template<class TV> void FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Read_State(TYPED_ISTREAM& input,const int state_index)
{
    for(COLLISION_GEOMETRY_ID i(0);i<collision_bodies.collision_geometry_collection.bodies.m;i++) if(collision_bodies.Is_Active(i)) collision_bodies.collision_geometry_collection.bodies(i)->Read_State(input,state_index);
}
//#####################################################################
// Function Write_State
//#####################################################################
template<class TV> void FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Write_State(TYPED_OSTREAM& output,const int state_index) const
{
    for(COLLISION_GEOMETRY_ID i(0);i<collision_bodies.collision_geometry_collection.bodies.m;i++) if(collision_bodies.Is_Active(i)) collision_bodies.collision_geometry_collection.bodies(i)->Write_State(output,state_index);
}
//#####################################################################
// Function Initialize_Grid_Structures
//#####################################################################
template<class TV> void FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Initialize_Grid_Structures(const GRID<TV>& grid_input,OBJECTS_IN_CELL<TV,COLLISION_GEOMETRY_ID>& objects_in_cell,const COLLISION_GEOMETRY_ID id) const
{
    if(&grid!=&grid_input) PHYSBAM_FATAL_ERROR();
    const_cast<FLUID_COLLISION_BODY_INACCURATE_UNION<TV>&>(*this).Initialize_Grid_Structures_Helper(objects_in_cell,id,typename GRID<TV>::GRID_TAG());
}
//#####################################################################
// Function Initialize_Grid_Structures_Helper
//#####################################################################
template<class TV> void FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Initialize_Grid_Structures_Helper(OBJECTS_IN_CELL<TV,COLLISION_GEOMETRY_ID>& objects_in_cell,const COLLISION_GEOMETRY_ID id,UNIFORM_TAG<TV>)
{
    collision_bodies.collision_geometry_collection.Update_Bounding_Boxes();
    // phi and velocity
    phi.Resize(grid.Cell_Indices(3),false,false);phi.Fill(10*grid.dX.Min());
    face_velocities.Resize(grid,3,false,false);face_velocities.Fill((T)0);
    face_velocities_set.Resize(grid,3,false,false);face_velocities_set.Fill(false);
    T_FACE_ARRAYS_INT face_velocities_count(grid,3);
    T_FACE_ARRAYS_COLLISION_GEOMETRY_ID face_operations(grid,3);
    for(COLLISION_GEOMETRY_ID i(0);i<collision_bodies.collision_geometry_collection.bodies.m;i++)if(collision_bodies.Is_Active(i) && collision_bodies.collision_geometry_collection.bodies(i)->active)
        Initialize_Grid_Structures_Subobject(face_velocities_count,face_operations,i,typename GRID<TV>::GRID_TAG());
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) if(face_velocities_count.Component(iterator.Axis())(iterator.Face_Index())){
        face_velocities.Component(iterator.Axis())(iterator.Face_Index())/=face_velocities_count.Component(iterator.Axis())(iterator.Face_Index());
        face_velocities_set.Component(iterator.Axis())(iterator.Face_Index())=true;}
}
//#####################################################################
// Function Initialize_Grid_Structures_Subobject
//#####################################################################
template<class TV> void FLUID_COLLISION_BODY_INACCURATE_UNION<TV>::
Initialize_Grid_Structures_Subobject(T_FACE_ARRAYS_INT& face_velocities_count,T_FACE_ARRAYS_COLLISION_GEOMETRY_ID& face_operations,const COLLISION_GEOMETRY_ID subobject,UNIFORM_TAG<TV>)
{
    COLLISION_GEOMETRY<TV>& collision_body=*collision_bodies.collision_geometry_collection.bodies(subobject);
    RANGE<TV> bounding_box=collision_body.Axis_Aligned_Bounding_Box();
    RANGE<TV_INT> box=grid.Clamp_To_Cell(bounding_box,2).Thickened(1);
    for(CELL_ITERATOR<TV> iterator(grid,box);iterator.Valid();iterator.Next()){
        T phi_value=collision_body.Implicit_Geometry_Extended_Value(iterator.Location());
        phi(iterator.Cell_Index())=min(phi_value,phi(iterator.Cell_Index()));
        if(phi_value<0) for(int axis=0;axis<TV::m;axis++){
            TV_INT face1=iterator.First_Face_Index(axis),face2=iterator.Second_Face_Index(axis);
            if(face_operations.Component(axis)(face1)!=subobject){face_operations.Component(axis)(face1)=subobject;
                face_velocities.Component(axis)(face1)=collision_body.Pointwise_Object_Velocity(grid.Face(FACE_INDEX<TV::m>(axis,face1)))[axis];face_velocities_count.Component(axis)(face1)++;}
            if(face_operations.Component(axis)(face2)!=subobject){face_operations.Component(axis)(face2)=subobject;
                face_velocities.Component(axis)(face2)=collision_body.Pointwise_Object_Velocity(grid.Face(FACE_INDEX<TV::m>(axis,face2)))[axis];face_velocities_count.Component(axis)(face2)++;}}}
}
//##################################################################### 
#define P(...) __VA_ARGS__
#define INSTANTIATION_HELPERV(TV) \
    template FLUID_COLLISION_BODY_INACCURATE_UNION<P(TV) >::FLUID_COLLISION_BODY_INACCURATE_UNION(P(GRID<TV>)&); \
    template void FLUID_COLLISION_BODY_INACCURATE_UNION<P(TV) >::Initialize_Grid_Structures(P(GRID<TV>) const&,OBJECTS_IN_CELL<P(TV),COLLISION_GEOMETRY_ID>&,COLLISION_GEOMETRY_ID) const;
#define INSTANTIATION_HELPER_T(T) \
    INSTANTIATION_HELPERV(P(VECTOR<T,1>)) INSTANTIATION_HELPERV(P(VECTOR<T,2>)) INSTANTIATION_HELPERV(P(VECTOR<T,3>))

INSTANTIATION_HELPER_T(float);
INSTANTIATION_HELPER_T(double);

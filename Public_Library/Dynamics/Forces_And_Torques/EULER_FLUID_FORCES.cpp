//#####################################################################
// Copyright 2007-2008, Nipun Kwatra, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_FLUID_FORCES
//#####################################################################
#include <Core/Data_Structures/HASHTABLE.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_1D.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_2D.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Dynamics/Forces_And_Torques/EULER_FLUID_FORCES.h>
using ::std::sqrt;
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> EULER_FLUID_FORCES<TV>::
EULER_FLUID_FORCES(const GRID<TV>& grid_input,const ARRAY<T,FACE_INDEX<TV::m> >& pressure_at_faces_input,
    const ARRAY<bool,FACE_INDEX<TV::m> >& solid_fluid_face_input,const ARRAY<bool,TV_INT>& cells_inside_fluid_input,
    const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>* collision_bodies_affecting_fluid_input,DEFORMABLE_PARTICLES<TV>& particles_input,
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input):SOLIDS_FORCES<TV>(particles_input,rigid_body_collection_input),
    grid(grid_input),pressure_at_faces(pressure_at_faces_input),solid_fluid_face(solid_fluid_face_input),
    cells_inside_fluid(cells_inside_fluid_input),collision_bodies_affecting_fluid(collision_bodies_affecting_fluid_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> EULER_FLUID_FORCES<TV>::
~EULER_FLUID_FORCES()
{}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void EULER_FLUID_FORCES<TV>::
Add_Velocity_Independent_Forces(GENERALIZED_VELOCITY<TV>& F,const T time) const
{
    COLLISION_GEOMETRY_ID body_id;int simplex_id;
    T distance,max_distance=grid.dX.Magnitude();
    TV one_over_dx=grid.one_over_dX;T cell_size=grid.Cell_Size();
    TV_INT face_index,first_cell_index,second_cell_index;TV first_cell_location,second_cell_location;int axis;
    bool first_cell_inside_fluid,second_cell_inside_fluid;
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){axis=iterator.Axis();face_index=iterator.Face_Index();
        if(!solid_fluid_face.Component(axis)(face_index)) continue;
        first_cell_index=iterator.First_Cell_Index();second_cell_index=iterator.Second_Cell_Index();
        first_cell_inside_fluid=cells_inside_fluid.Valid_Index(first_cell_index) && cells_inside_fluid(first_cell_index);
        second_cell_inside_fluid=cells_inside_fluid.Valid_Index(second_cell_index) && cells_inside_fluid(second_cell_index);
        int direction=0;
        if(first_cell_inside_fluid && !second_cell_inside_fluid) direction=1; 
        else if(!first_cell_inside_fluid && second_cell_inside_fluid) direction=-1;
        if(direction!=0){ // face on solid-fluid boundary
            TV location=iterator.Location();
            collision_bodies_affecting_fluid->collision_geometry_collection.Closest_Boundary_Point(location,max_distance,distance,body_id,simplex_id);
            const COLLISION_GEOMETRY<TV>& collision_body=collision_bodies_affecting_fluid->collision_geometry_collection(body_id);
            if(const RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_geometry=dynamic_cast<const RIGID_COLLISION_GEOMETRY<TV>*>(&collision_body)){
                // apply force and torque to this body from the pressure flux times area
                const RIGID_BODY<TV>& rigid_body=rigid_collision_geometry->rigid_body;
                int rigid_body_index=rigid_body.particle_index;
                T face_area=cell_size*one_over_dx[axis],face_pressure=pressure_at_faces.Component(axis)(face_index);
                TV center_of_mass=rigid_body.Frame().t,force=face_area*face_pressure*direction*TV::Axis_Vector(axis);
                F.rigid_V.array(rigid_body_index).linear+=force;
                F.rigid_V.array(rigid_body_index).angular+=TV::Cross_Product(location-center_of_mass,force);}
            else PHYSBAM_FATAL_ERROR("deformable part not implemented");}}
}
//#####################################################################
namespace PhysBAM{
template class EULER_FLUID_FORCES<VECTOR<float,1> >;
template class EULER_FLUID_FORCES<VECTOR<float,2> >;
template class EULER_FLUID_FORCES<VECTOR<float,3> >;
template class EULER_FLUID_FORCES<VECTOR<double,1> >;
template class EULER_FLUID_FORCES<VECTOR<double,2> >;
template class EULER_FLUID_FORCES<VECTOR<double,3> >;
}

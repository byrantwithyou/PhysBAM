//#####################################################################
// Copyright 2009, Nipun Kwatra, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_BOUNDARY_CONDITION_COLLISIONS
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION_COLLISIONS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> IMPLICIT_BOUNDARY_CONDITION_COLLISIONS<TV>::
IMPLICIT_BOUNDARY_CONDITION_COLLISIONS(COLLISION_BODY_COLLECTION<TV>& collision_geometry_collection_input,
    const bool use_implicit_geometry_input)
    :collision_geometry_collection(collision_geometry_collection_input),use_implicit_geometry(use_implicit_geometry_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> IMPLICIT_BOUNDARY_CONDITION_COLLISIONS<TV>::
~IMPLICIT_BOUNDARY_CONDITION_COLLISIONS()
{}
//#####################################################################
// Function Update_Boundary_Conditions
//#####################################################################
template<class TV> void IMPLICIT_BOUNDARY_CONDITION_COLLISIONS<TV>::
Update_Boundary_Conditions(const GRID<TV>& grid,ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,ARRAY<T,TV_INT>& p,
    ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T time)
{
    T p_inside_solid=0; // TODO: set to something nasty and make sure uncovered cells get good values for initial guess.

    if(use_implicit_geometry){
        COLLISION_GEOMETRY_ID body_id;
        for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next()){TV location=iterator.Location();
            TV_INT cell_index=iterator.Cell_Index();
            if(collision_geometry_collection.Implicit_Geometry_Lazy_Inside_Any_Body(location,body_id)){
                psi_D(cell_index)=true;p(cell_index)=p_inside_solid;}}
        return;}

    for(COLLISION_GEOMETRY_ID i(0);i<collision_geometry_collection.Size();i++) if(collision_geometry_collection.Is_Active(i)){
        T collision_thickness_over_two=collision_geometry_collection.collision_body_thickness*(T).5;
        if(DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>* object=dynamic_cast<DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>*>(&collision_geometry_collection(i))){
            if(!object->volume_object){
                object->object.Update_Bounding_Box();
                object->object.Refresh_Auxiliary_Structures();
                for(CELL_ITERATOR<TV> iterator(grid,grid.Clamp_To_Cell(*object->object.bounding_box,1));iterator.Valid();iterator.Next())
                    if(object->Inside(iterator.Location(),collision_thickness_over_two)){
                        psi_D(iterator.Cell_Index())=true;p(iterator.Cell_Index())=p_inside_solid;}}
            else{
                for(int i=0;i<object->volume_object->mesh.elements.m;i++){T_SIMPLEX simplex=object->volume_object->Get_Element(i);
                    for(CELL_ITERATOR<TV> iterator(grid,grid.Clamp_To_Cell(simplex.Bounding_Box(),1));iterator.Valid();iterator.Next())
                        if(simplex.Inside(iterator.Location(),collision_thickness_over_two)){
                            psi_D(iterator.Cell_Index())=true;p(iterator.Cell_Index())=p_inside_solid;}}}}
        else if(RIGID_COLLISION_GEOMETRY_BASE<TV>* object=dynamic_cast<RIGID_COLLISION_GEOMETRY_BASE<TV>*>(&collision_geometry_collection(i))){
            RIGID_BODY<TV>& rigid_body=dynamic_cast<RIGID_BODY<TV>&>(object->rigid_body);
            if(!rigid_body.simplicial_object->mesh.incident_elements) rigid_body.simplicial_object->mesh.Initialize_Incident_Elements();
            if(!rigid_body.simplicial_object->mesh.adjacent_elements) rigid_body.simplicial_object->mesh.Initialize_Adjacent_Elements();
            object->Update_Bounding_Box();
            for(CELL_ITERATOR<TV> iterator(grid,grid.Clamp_To_Cell(object->Axis_Aligned_Bounding_Box(),1));iterator.Valid();iterator.Next()){
                if(collision_geometry_collection(i).Inside(iterator.Location(),collision_thickness_over_two)){
                    psi_D(iterator.Cell_Index())=true;p(iterator.Cell_Index())=p_inside_solid;}}}
        else PHYSBAM_FATAL_ERROR("Unrecognized collision body type");}
}
//#####################################################################
namespace PhysBAM{
template class IMPLICIT_BOUNDARY_CONDITION_COLLISIONS<VECTOR<float,1> >;
template class IMPLICIT_BOUNDARY_CONDITION_COLLISIONS<VECTOR<float,2> >;
template class IMPLICIT_BOUNDARY_CONDITION_COLLISIONS<VECTOR<float,3> >;
template class IMPLICIT_BOUNDARY_CONDITION_COLLISIONS<VECTOR<double,1> >;
template class IMPLICIT_BOUNDARY_CONDITION_COLLISIONS<VECTOR<double,2> >;
template class IMPLICIT_BOUNDARY_CONDITION_COLLISIONS<VECTOR<double,3> >;
}

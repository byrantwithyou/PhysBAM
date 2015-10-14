//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Joseph Teran, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Deformables/Collisions_And_Interactions/COLLISION_PENALTY_FORCES.h>
#include <Deformables/Collisions_And_Interactions/TETRAHEDRON_COLLISION_BODY.h>
namespace PhysBAM{
using ::std::exp;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COLLISION_PENALTY_FORCES<TV>::
COLLISION_PENALTY_FORCES(DEFORMABLE_PARTICLES<TV>& particles)
    :DEFORMABLES_FORCES<TV>(particles),collision_body_list_id(0)
{
    Set_Stiffness();Set_Separation_Parameter();Set_Self_Collision_Reciprocity_Factor();Set_Perform_Self_Collision();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COLLISION_PENALTY_FORCES<TV>::
~COLLISION_PENALTY_FORCES()
{
}
//#####################################################################
// Function Set_Collision_Body_List
//#####################################################################
template<class TV> void COLLISION_PENALTY_FORCES<TV>::
Set_Collision_Body_List(COLLISION_BODY_COLLECTION<TV>& collision_body_list_input)
{
    collision_body_list=&collision_body_list_input;
    skip_collision_body.Resize(collision_body_list->bodies.m);
    skip_collision_body.Fill(false);
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void COLLISION_PENALTY_FORCES<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void COLLISION_PENALTY_FORCES<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void COLLISION_PENALTY_FORCES<TV>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
    if(collision_body_list_id>=COLLISION_GEOMETRY_ID(0) && typeid((*collision_body_list)(collision_body_list_id))==typeid(TETRAHEDRON_COLLISION_BODY<T>)){
        TETRAHEDRON_COLLISION_BODY<T>& collision_body=(TETRAHEDRON_COLLISION_BODY<T>&)((*collision_body_list)(collision_body_list_id));
        collision_body.tetrahedralized_volume.hierarchy->Update_Boxes(collision_body.collision_thickness);
        collision_body.tetrahedralized_volume.triangulated_surface->hierarchy->Update_Boxes();
        collision_body.tetrahedralized_volume.triangulated_surface->Update_Vertex_Normals();}
}
//#####################################################################
// Function Update_Forces_And_Derivatives
//#####################################################################
template<class TV> void COLLISION_PENALTY_FORCES<TV>::
Update_Forces_And_Derivatives()
{
    for(int p=0;p<check_collision.m;p++){
        int index=check_collision(p);collision_force(p)=TV();collision_force_derivative(p)=SYMMETRIC_MATRIX<T,TV::m>();
        for(COLLISION_GEOMETRY_ID r(0);r<collision_body_list->bodies.m;r++) if(collision_body_list->Is_Active(r)){
                COLLISION_GEOMETRY<TV>& collision_body=*collision_body_list->bodies(r);
                if(!skip_collision_body(r) && (perform_self_collision || collision_body.collision_geometry_id!=collision_body_list_id)){
                    int collision_body_particle_index=-1;if(collision_body_list_id==collision_body.collision_geometry_id) collision_body_particle_index=index;
                    T phi_value;int aggregate=-1;TV normal=collision_body.Implicit_Geometry_Extended_Normal(particles.X(index),phi_value,aggregate,collision_body_particle_index);
                    T scaled_stiffness=stiffness;if(collision_body_list_id==collision_body.collision_geometry_id) scaled_stiffness*=self_collision_reciprocity_factor;
                    if(phi_value<=0){
                        collision_force(p)+=stiffness*(-phi_value+separation_parameter)*normal;
                        collision_force_derivative(p)-=scaled_stiffness*SYMMETRIC_MATRIX<T,TV::m>::Outer_Product(normal);}
                    else if(phi_value<collision_body.collision_thickness){
                        collision_force(p)+=stiffness*separation_parameter*exp(-phi_value/separation_parameter)*normal;
                        collision_force_derivative(p)-=scaled_stiffness*exp(-phi_value/separation_parameter)*SYMMETRIC_MATRIX<T,TV::m>::Outer_Product(normal);}}}}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void COLLISION_PENALTY_FORCES<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(int p=0;p<check_collision.m;p++) F(check_collision(p))+=collision_force(p);
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void COLLISION_PENALTY_FORCES<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void COLLISION_PENALTY_FORCES<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(int p=0;p<check_collision.m;p++) F(check_collision(p))+=collision_force_derivative(p)*V(check_collision(p));
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> auto COLLISION_PENALTY_FORCES<TV>::
CFL_Strain_Rate() const -> T
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template class COLLISION_PENALTY_FORCES<VECTOR<double,1> >;
template class COLLISION_PENALTY_FORCES<VECTOR<double,2> >;
template class COLLISION_PENALTY_FORCES<VECTOR<double,3> >;
template class COLLISION_PENALTY_FORCES<VECTOR<float,1> >;
template class COLLISION_PENALTY_FORCES<VECTOR<float,2> >;
template class COLLISION_PENALTY_FORCES<VECTOR<float,3> >;
}

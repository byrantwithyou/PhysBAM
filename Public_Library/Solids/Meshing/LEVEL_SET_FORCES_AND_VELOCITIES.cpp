//#####################################################################
// Copyright 2002-2007, Christopher Allocco, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <Solids/Meshing/LEVEL_SET_FORCES_AND_VELOCITIES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LEVEL_SET_FORCES_AND_VELOCITIES<TV>::
LEVEL_SET_FORCES_AND_VELOCITIES(T_MESH_OBJECT& mesh_object,IMPLICIT_OBJECT<TV>& implicit_object)
    :mesh_object(mesh_object),implicit_object(implicit_object)
{
    Set_Force_Attraction_Coefficient();Set_Velocity_Attraction_Coefficient();
    Use_External_Forces();
    Use_External_Velocities_Normal_To_Boundary();
    Allow_Tangential_Velocity_Slip();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LEVEL_SET_FORCES_AND_VELOCITIES<TV>::
~LEVEL_SET_FORCES_AND_VELOCITIES()
{
}
//#####################################################################
// Function Add_External_Forces
//#####################################################################
template<class TV> void LEVEL_SET_FORCES_AND_VELOCITIES<TV>::
Add_External_Forces(ARRAY_VIEW<TV> F,const T time)
{
    if(!use_external_forces) return;

    DEFORMABLE_PARTICLES<TV>& particles=dynamic_cast<DEFORMABLE_PARTICLES<TV>&>(mesh_object.particles);
    const ARRAY<int>& boundary_nodes=*mesh_object.mesh.boundary_nodes;
    T_BOUNDARY_OBJECT& boundary_object=mesh_object.Get_Boundary_Object();

    boundary_object.Update_Vertex_Normals();
    for(int i=0;i<boundary_nodes.m;i++){int p=boundary_nodes(i);
        F(p)-=force_attraction_coefficient*sqrt(particles.mass(p))*implicit_object(particles.X(p))*(*boundary_object.vertex_normals)(p);}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class TV> void LEVEL_SET_FORCES_AND_VELOCITIES<TV>::
Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time)
{
    if(!use_external_velocities) return;

    GEOMETRY_PARTICLES<TV>& particles=mesh_object.particles;
    const ARRAY<int>& boundary_nodes=*mesh_object.mesh.boundary_nodes;

    if(allow_tangential_velocity_slip){
        if(use_external_velocities_normal_to_boundary){
            T_BOUNDARY_OBJECT& boundary_object=mesh_object.Get_Boundary_Object();
            boundary_object.Update_Vertex_Normals();
            for(int i=0;i<boundary_nodes.m;i++){int p=boundary_nodes(i);
                TV N=(*boundary_object.vertex_normals)(p);
                V(p)-=(TV::Dot_Product(V(p),N)+velocity_attraction_coefficient*implicit_object(particles.X(p)))*N;}}
        else // use external velocities normal to level set
            for(int i=0;i<boundary_nodes.m;i++){int p=boundary_nodes(i);
                TV N=implicit_object.Normal(particles.X(p));
                V(p)-=(TV::Dot_Product(V(p),N)+velocity_attraction_coefficient*implicit_object(particles.X(p)))*N;}}
    else{
        if(use_external_velocities_normal_to_boundary){
            T_BOUNDARY_OBJECT& boundary_object=mesh_object.Get_Boundary_Object();
            boundary_object.Update_Vertex_Normals();
            for(int i=0;i<boundary_nodes.m;i++){int p=boundary_nodes(i);
                TV N=(*boundary_object.vertex_normals)(p);
                V(p)=-velocity_attraction_coefficient*implicit_object(particles.X(p))*N;}}
        else // use external velocities normal to level set
            for(int i=0;i<boundary_nodes.m;i++){int p=boundary_nodes(i);
                TV N=implicit_object.Normal(particles.X(p));
                V(p)=-velocity_attraction_coefficient*implicit_object(particles.X(p))*N;}}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
template<class TV> void LEVEL_SET_FORCES_AND_VELOCITIES<TV>::
Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time)
{
    if(!use_external_velocities) return;

    GEOMETRY_PARTICLES<TV>& particles=mesh_object.particles;
    const ARRAY<int>& boundary_nodes=*mesh_object.mesh.boundary_nodes;

    if(allow_tangential_velocity_slip){
        if(use_external_velocities_normal_to_boundary){
            T_BOUNDARY_OBJECT& boundary_object=mesh_object.Get_Boundary_Object();
            boundary_object.Update_Vertex_Normals();
            for(int i=0;i<boundary_nodes.m;i++){int p=boundary_nodes(i);
                TV N=(*boundary_object.vertex_normals)(p);
                V(p)-=TV::Dot_Product(V(p),N)*N;}}
        else // use external velocities normal to level set
            for(int i=0;i<boundary_nodes.m;i++){int p=boundary_nodes(i);
                TV N=implicit_object.Normal(particles.X(p));
                V(p)-=TV::Dot_Product(V(p),N)*N;}}
    else for(int i=0;i<boundary_nodes.m;i++){int p=boundary_nodes(i);
        V(p)=TV();}
}
template class LEVEL_SET_FORCES_AND_VELOCITIES<VECTOR<double,3> >;
template class LEVEL_SET_FORCES_AND_VELOCITIES<VECTOR<float,3> >;
}

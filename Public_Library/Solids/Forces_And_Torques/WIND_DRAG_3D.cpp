//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Grid_PDE/Interpolation/LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Forces_And_Torques/WIND_DRAG_3D.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
using namespace PhysBAM;
//#####################################################################
template<class T> const LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM<VECTOR<T,3>,T> WIND_DRAG_3D<T>::interpolation;
template<class T> const LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM<VECTOR<T,3>,VECTOR<T,3> > WIND_DRAG_3D<T>::vector_interpolation;
//#####################################################################
// Constructor
//#####################################################################
namespace{
template<class T> TRIANGULATED_SURFACE<T>& Triangulated_Surface(TRIANGULATED_SURFACE<T>& triangulated_surface)
{
    return triangulated_surface;
}
template<class T> TRIANGULATED_SURFACE<T>& Triangulated_Surface(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume)
{
    if(!tetrahedralized_volume.triangulated_surface) tetrahedralized_volume.Initialize_Triangulated_Surface();
    return *tetrahedralized_volume.triangulated_surface;
}}
template<class T> template<class T_OBJECT> WIND_DRAG_3D<T>::
WIND_DRAG_3D(T_OBJECT& object,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
    :SOLIDS_FORCES<TV>(dynamic_cast<DEFORMABLE_PARTICLES<TV>&>(object.particles),rigid_body_collection_input),triangulated_surface(&Triangulated_Surface(object)),rigid_body(0),use_constant_wind(false),
    use_spatially_varying_wind(false),spatially_varying_wind(0),wind_density(0),spatially_varying_wind_density(0),spatially_varying_wind_pressure(0),linear_normal_viscosity(0),mpi_solids(0)
{}
template<class T> WIND_DRAG_3D<T>::
WIND_DRAG_3D(RIGID_BODY<TV>& rigid_body_input,DEFORMABLE_PARTICLES<TV>& deformable_body_particles_input)
    :SOLIDS_FORCES<TV>(deformable_body_particles_input,rigid_body_input.rigid_body_collection),triangulated_surface(0),rigid_body(&rigid_body_input),use_constant_wind(false),
    use_spatially_varying_wind(false),spatially_varying_wind(0),wind_density(0),spatially_varying_wind_density(0),spatially_varying_wind_pressure(0),linear_normal_viscosity(0),mpi_solids(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> WIND_DRAG_3D<T>::
~WIND_DRAG_3D()
{}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class T> void WIND_DRAG_3D<T>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    if(triangulated_surface) triangulated_surface->mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class T> void WIND_DRAG_3D<T>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,const ARRAY<bool>& rigid_particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    const auto& e=triangulated_surface->mesh.elements;
    Update_Force_Particles(force_particles,e.Flattened(),particle_is_simulated,true);
    Update_Force_Elements(force_elements,e,particle_is_simulated);
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T> void WIND_DRAG_3D<T>::
Update_Position_Based_State(const T time)
{
    T_SIMPLICIAL_OBJECT& simplicial_object=triangulated_surface?*triangulated_surface:*rigid_body->simplicial_object;
    T wind_viscosity=use_constant_wind?constant_wind_viscosity:spatially_varying_wind_viscosity;
    if(wind_viscosity || spatially_varying_wind_pressure || wind_density || spatially_varying_wind_density){
        optimization.Resize(simplicial_object.mesh.elements.m,no_init);
        if(triangulated_surface) for(int t:force_elements){
            int i,j,k;simplicial_object.mesh.elements(t).Get(i,j,k);
            const TV &Xi=particles.X(i),&Xj=particles.X(j),&Xk=particles.X(k);
            optimization(t).center=((T)1/3)*(Xi+Xj+Xk);optimization(t).inward_normal=PLANE<T>::Normal_Direction(Xi,Xk,Xj);
            optimization(t).one_third_area=((T)1/6)*optimization(t).inward_normal.Normalize();
            if(use_spatially_varying_wind) optimization(t).wind_velocity=Spatially_Varying_Wind_Velocity(optimization(t).center);}
        else for(int t=0;t<rigid_body->simplicial_object->mesh.elements.m;t++){
            T_SIMPLEX world_space_simplex=rigid_body->World_Space_Simplex(t);
            optimization(t).center=world_space_simplex.Center(),optimization(t).inward_normal=-world_space_simplex.Normal();
            optimization(t).one_third_area=world_space_simplex.Size();
            if(use_spatially_varying_wind) optimization(t).wind_velocity=Spatially_Varying_Wind_Velocity(optimization(t).center);}}
    if(linear_normal_viscosity && triangulated_surface){ // compute vertex normals for fragment
        vertex_normals.Resize(particles.Size(),init_all);
        for(int t:force_elements){
            const VECTOR<int,3>& nodes=simplicial_object.mesh.elements(t);
            vertex_normals.Subset(nodes)+=PLANE<T>::Normal_Direction(particles.X.Subset(nodes));}
        for(int p:force_particles) vertex_normals(p).Normalize();}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T> VECTOR<T,3> WIND_DRAG_3D<T>::
Add_Velocity_Independent_Forces_Helper(TV relative_velocity,int t) const
{
    // wind drag pressure - per unit area
    T wind_viscosity=use_constant_wind?constant_wind_viscosity:spatially_varying_wind_viscosity;
    TV force=wind_viscosity*(relative_velocity-TV::Dot_Product(relative_velocity,optimization(t).inward_normal)*optimization(t).inward_normal);
    // wind pressure 
    T pressure;
    if(spatially_varying_wind_pressure) pressure=Spatially_Varying_Wind_Pressure(optimization(t).center);
    else{
        T normal_velocity=TV::Dot_Product(relative_velocity,optimization(t).inward_normal);
        if(normal_velocity>0) pressure=sqr(normal_velocity);else pressure=-sqr(normal_velocity);
        if(wind_density) pressure*=wind_density;
        else pressure*=Spatially_Varying_Wind_Density(optimization(t).center);}
    force+=pressure*optimization(t).inward_normal;
    force*=optimization(t).one_third_area;
    return force;
}

template<class T> void WIND_DRAG_3D<T>::
Add_Velocity_Independent_Forces(GENERALIZED_VELOCITY<TV>& F,const T time) const
{
    if(use_constant_wind==use_spatially_varying_wind) PHYSBAM_FATAL_ERROR();
    T wind_viscosity=use_constant_wind?constant_wind_viscosity:spatially_varying_wind_viscosity;
    if(wind_viscosity || spatially_varying_wind_pressure || wind_density || spatially_varying_wind_density){
        if(!spatially_varying_wind_pressure && !wind_density && !spatially_varying_wind_density) PHYSBAM_FATAL_ERROR();
        if(triangulated_surface) for(int t:force_elements){
            int i,j,k;triangulated_surface->mesh.elements(t).Get(i,j,k);
            TV wind_velocity=use_constant_wind?constant_wind:optimization(t).wind_velocity;
            TV relative_velocity=wind_velocity-((T)1/3)*(particles.V(i)+particles.V(j)+particles.V(k));
            TV triangle_force=Add_Velocity_Independent_Forces_Helper(relative_velocity,t);
            F.V.array(i)+=triangle_force;
            F.V.array(j)+=triangle_force;
            F.V.array(k)+=triangle_force;}
        else{
            if(rigid_body->Has_Infinite_Inertia() || rigid_body->particle_index<0) return;
            TWIST<TV>& wrench=F.rigid_V.array(rigid_body->particle_index);
            for(int t=0;t<rigid_body->simplicial_object->mesh.elements.m;t++){
                TV wind_velocity=use_constant_wind?constant_wind:optimization(t).wind_velocity;
                TV relative_velocity=wind_velocity-rigid_body->Pointwise_Object_Velocity(optimization(t).center);
                TV simplex_force=Add_Velocity_Independent_Forces_Helper(relative_velocity,t);
                wrench.linear+=simplex_force;wrench.angular+=TV::Cross_Product(optimization(t).center-rigid_body->Frame().t,simplex_force);}}}
    if(linear_normal_viscosity && triangulated_surface){
        TV wind_acceleration=linear_normal_viscosity*constant_wind;
        for(int p:force_particles){
            if(use_constant_wind) F.V.array(p)+=particles.mass(p)*TV::Dot_Product(wind_acceleration,vertex_normals(p))*vertex_normals(p);
            else F.V.array(p)+=particles.mass(p)*linear_normal_viscosity*TV::Dot_Product(Spatially_Varying_Wind_Velocity(particles.X(p)),vertex_normals(p))*vertex_normals(p);}}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class T> void WIND_DRAG_3D<T>::
Add_Velocity_Dependent_Forces(const GENERALIZED_VELOCITY<TV>& V,GENERALIZED_VELOCITY<TV>& F,const T time) const
{
    if(use_constant_wind==use_spatially_varying_wind) PHYSBAM_FATAL_ERROR();
    T wind_viscosity=use_constant_wind?constant_wind_viscosity:spatially_varying_wind_viscosity;
    if(wind_viscosity){
        if(triangulated_surface) for(int t:force_elements){
            int i,j,k;triangulated_surface->mesh.elements(t).Get(i,j,k);
            // wind drag pressure - per unit area
            TV negative_V=-((T)1/3)*(V.V.array(i)+V.V.array(j)+V.V.array(k));
            TV triangle_force=wind_viscosity*(negative_V-TV::Dot_Product(negative_V,optimization(t).inward_normal)*optimization(t).inward_normal);
            // total force
            triangle_force*=optimization(t).one_third_area; // one third of the triangle force is distriduted to each node
            F.V.array(i)+=triangle_force;F.V.array(j)+=triangle_force;F.V.array(k)+=triangle_force;}
        else{
            if(rigid_body->Has_Infinite_Inertia() || rigid_body->particle_index<0) return;
            TWIST<TV>& wrench=F.rigid_V.array(rigid_body->particle_index);
            for(int t=0;t<rigid_body->simplicial_object->mesh.elements.m;t++){
                TV negative_V=-rigid_body->Pointwise_Object_Velocity(optimization(t).center);
                TV simplex_force=wind_viscosity*(negative_V-TV::Dot_Product(negative_V,optimization(t).inward_normal)*optimization(t).inward_normal);
                // total force
                simplex_force*=optimization(t).one_third_area; // one third of the triangle force is distriduted to each node
                wrench.linear+=simplex_force;wrench.angular+=TV::Cross_Product(optimization(t).center-rigid_body->Frame().t,simplex_force);}}}
    if(linear_normal_viscosity && triangulated_surface){
        for(int p:force_particles){
            F.V.array(p)-=particles.mass(p)*linear_normal_viscosity*TV::Dot_Product(V.V.array(p),vertex_normals(p))*vertex_normals(p);}}
}
//#####################################################################
namespace PhysBAM{
#define INSTANTIATION_HELPER(T) \
    template class WIND_DRAG_3D<T>; \
    template WIND_DRAG_3D<T>::WIND_DRAG_3D(TRIANGULATED_SURFACE<T>&,RIGID_BODY_COLLECTION<VECTOR<T,3> >&); \
    template WIND_DRAG_3D<T>::WIND_DRAG_3D(TETRAHEDRALIZED_VOLUME<T>&,RIGID_BODY_COLLECTION<VECTOR<T,3> >&);
INSTANTIATION_HELPER(float)
INSTANTIATION_HELPER(double)
}

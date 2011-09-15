//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FIRE_MELTING_EXAMPLE_S3D
//#####################################################################
#ifndef __FIRE_MELTING_EXAMPLE_S3D__
#define __FIRE_MELTING_EXAMPLE_S3D__

#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <Fracture/TRIANGLES_OF_MATERIAL_3D.h>
#include <Heat_Flows/HEAT_3D.h>
#include <Level_Sets/LEVELSET_TRIANGULATED_OBJECT.h>
#include <Solids_And_Fluids/MELTING_EXAMPLE_S3D.h>
namespace PhysBAM{

template<class T,class RW>
class FIRE_MELTING_EXAMPLE_S3D:public MELTING_EXAMPLE_S3D<T,RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::solids_parameters;using MELTING_EXAMPLE_S3D<T,RW>::melting_parameters;
    using MELTING_EXAMPLE_S3D<T,RW>::Update_Solids_Topology_For_Melting;

    bool initialized;
    T cloth_temperature_time_constant;
    T cloth_ignition_temperature;
    T reaction_start,reaction_peak_start,reaction_peak_stop,reaction_stop;
    T divergence_source_bandwidth;
    T divergence_source_rate;
    T phi_source_bandwidth;
    T phi_source_depth;
    int temperature_smoothing_steps;
    T negative_reaction_coefficient_multiplier;
    bool phi_source_smooth_start;
    LEVELSET_2D<T>* inflammable_levelset;

    FIRE_MELTING_EXAMPLE_S3D()
        :MELTING_EXAMPLE_S3D<T,RW>(fluids_parameters.FIRE),initialized(false),
        cloth_temperature_time_constant(2),cloth_ignition_temperature(300),reaction_start(0),reaction_peak_start(.5),reaction_peak_stop(1.5),reaction_stop(2),
        divergence_source_bandwidth(1.5),divergence_source_rate(1),phi_source_bandwidth(2.5),phi_source_depth(1.5),temperature_smoothing_steps(5),
        negative_reaction_coefficient_multiplier(0),phi_source_smooth_start(true),inflammable_levelset(0)
    {
        fluids_parameters.particle_levelset_evolution.particle_levelset.reincorporate_removed_particles_everywhere=true;
        fluids_parameters.use_non_zero_divergence=true;
    }

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    MELTING_EXAMPLE_S3D<T,RW>::Initialize_Bodies();
    if(!initialized){initialized=true;
        for(int object=1;object<=melting_parameters.levelsets.m;object++){
            melting_parameters.temperature(object)=new ARRAY<T>(melting_parameters.levelsets(object)->grid.number_of_nodes);
            ARRAY<T>::copy(fluids_parameters.temperature_container.ambient_temperature,*melting_parameters.temperature(object));
            melting_parameters.reaction(object)=new ARRAY<T>(melting_parameters.levelsets(object)->grid.number_of_nodes);
            ARRAY<T>::copy(0,*melting_parameters.reaction(object));}}
}   
//#####################################################################
// Function Melting_Substep
//#####################################################################
void Melting_Substep(const T dt,const T time)
{
    LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION<T,T> interpolation;
    GRID<TV>& grid=fluids_parameters.grid;
    ARRAY<T,VECTOR<int,3> > smoothed_temperature(fluids_parameters.temperature_container.array_3d);
    HEAT_3D<T>::Smooth(grid,smoothed_temperature,temperature_smoothing_steps);
    for(int object=1;object<=melting_parameters.levelsets.m;object++){
        LEVELSET_TRIANGULATED_OBJECT<T,VECTOR<T,3> >& levelset=*melting_parameters.levelsets(object);
        DEFORMABLE_OBJECT_3D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(melting_parameters.body_index(object));
        PARTICLES<T,VECTOR<T,3> >& particles=deformable_object.particles;
        ARRAY<T>& temperature=*melting_parameters.temperature(object);ARRAY<T> old_temperature(temperature);
        ARRAY<T>& reaction=*melting_parameters.reaction(object);
        for(int n=1;n<=temperature.m;n++){
            int p=levelset.node_to_particle_mapping(n);if(!p)continue;
            temperature(n)-=dt*cloth_temperature_time_constant*(temperature(n)-interpolation.Clamped_To_Array(grid,smoothed_temperature,particles.X(p)));}
        LOG::cout<<"maximum cloth temperature = "<<temperature.Max()<<std::endl;
        for(int n=1;n<=reaction.m;n++){
            if(reaction(n)<=0){
                if(temperature(n)>cloth_ignition_temperature) reaction(n)=dt*(temperature(n)-cloth_ignition_temperature)/(temperature(n)-old_temperature(n));
                else reaction(n)=negative_reaction_coefficient_multiplier*(temperature(n)-cloth_ignition_temperature);}
            else reaction(n)+=dt;}
        LOG::cout<<"maximum cloth reaction = "<<reaction.Max()<<std::endl;}
    Update_Solids_Topology_For_Melting(dt,time);
}
//#####################################################################
// Function Melting_Levelset_Substep
//#####################################################################
void Melting_Levelset_Substep(const int object,const T dt,const T time)
{
    ARRAY<T>& phi=melting_parameters.levelsets(object)->phi;
    ARRAY<T>& reaction=*melting_parameters.reaction(object);
    for(int n=1;n<=reaction.m;n++)phi(n)=reaction(n)-reaction_stop;
    if(inflammable_levelset){
        ARRAY<VECTOR_2D<T> >& node_locations=melting_parameters.levelsets(object)->grid.Node_Locations();
        for(int n=1;n<=node_locations.m;n++)phi(n)=min(phi(n),inflammable_levelset->Phi(node_locations(n)));}
}
//#####################################################################
// Function Get_Divergence
//#####################################################################
void Get_Divergence(ARRAY<T,VECTOR<int,3> >& divergence,const T dt,const T time)
{
    ARRAY<T,VECTOR<int,3> >::copy(0,divergence);
    GRID<TV>& p_grid=fluids_parameters.p_grid;
    T bandwidth=divergence_source_bandwidth*p_grid.max_dx_dy_dz;
    T divergence_max=dt*divergence_source_rate/p_grid.max_dx_dy_dz;
    for(int object=1;object<=melting_parameters.levelsets.m;object++){
        LEVELSET_TRIANGULATED_OBJECT<T,VECTOR<T,3> >& levelset=*melting_parameters.levelsets(object);
        DEFORMABLE_OBJECT_3D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(melting_parameters.body_index(object));
        PARTICLES<T,VECTOR<T,3> >& particles=deformable_object.particles;
        TRIANGULATED_SURFACE<T>& material_surface=deformable_object.triangles_of_material->material_surface;
        PARTICLES<T,VECTOR<T,3> >& material_particles=material_surface.particles;
        ARRAY<VECTOR_2D<T> >& node_locations=levelset.grid.Node_Locations();
        ARRAY<T>& reaction=*melting_parameters.reaction(object);
        // absurd hack: replace velocities with reactions and use Update_Particle_Velocities to interpolate them to the material surface
        ARRAY<VECTOR<T,3> > V_save(particles.array_collection->Size());ARRAY<VECTOR<T,3> >::Exchange_Arrays(V_save,particles.V.array);
        for(int n=1;n<=reaction.m;n++)if(levelset.node_to_particle_mapping(n)){
            particles.V(levelset.node_to_particle_mapping(n)).x=reaction(n);
            if(inflammable_levelset) particles.V(levelset.node_to_particle_mapping(n)).y=inflammable_levelset->Phi(node_locations(n));}
        deformable_object.triangles_of_material->Update_Particle_Velocities();
        // compute divergence contribution
        material_surface.Update_Triangle_List(); // inefficient if positions haven't changed
        if(!material_surface.triangle_hierarchy) material_surface.Initialize_Triangle_Hierarchy();else material_surface.triangle_hierarchy->Update_Boxes();
        ARRAY<bool,VECTOR<int,3> > occupied(p_grid);deformable_object.collisions.Compute_Occupied_Cells(p_grid,occupied,false,bandwidth,0);
        for(int i=1;i<=p_grid.m;i++)for(int j=1;j<=p_grid.n;j++)for(int ij=1;ij<=p_grid.mn;ij++)if(occupied(i,j,ij)){
            int triangle;T distance;VECTOR<T,3> X=material_surface.Surface(p_grid.X(i,j,ij),bandwidth,1e-6,&triangle,&distance);
            if(distance>bandwidth) continue;
            int ti,tj,tk;material_surface.triangle_mesh.triangles.Get(triangle,ti,tj,tk);
            VECTOR<T,3> weights=TRIANGLE_3D<T>::Clamped_Barycentric_Coordinates(X,material_particles.X(ti),material_particles.X(tj),material_particles.X(tk));
            if(inflammable_levelset){
                distance-=min((T)0,weights.x*material_particles.V(ti).y+weights.y*material_particles.V(tj).y+weights.z*material_particles.V(tk).y); // V.y is inflammable phi
                if(distance>bandwidth) continue;}
            T r=weights.x*material_particles.V(ti).x+weights.y*material_particles.V(tj).x+weights.z*material_particles.V(tk).x; // V.x is reaction
            if(r<reaction_start) continue;
            else if(r<reaction_peak_start) divergence(i,j,ij)+=divergence_max*(bandwidth-distance)*(r-reaction_start)/(reaction_peak_start-reaction_start);
            else if(r<reaction_peak_stop) divergence(i,j,ij)+=divergence_max*(bandwidth-distance);
            else if(r<reaction_stop) divergence(i,j,ij)+=divergence_max*(bandwidth-distance)*(reaction_stop-r)/(reaction_stop-reaction_peak_stop);}
        // undo absurd hack
        ARRAY<VECTOR<T,3> >::Exchange_Arrays(V_save,particles.V.array);deformable_object.triangles_of_material->Update_Particle_Velocities();} 
    LOG::cout<<"divergence sum = "<<ARRAY<T,VECTOR<int,3> >::sum(divergence)<<std::endl;
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
void Adjust_Phi_With_Sources(const T time)
{
    GRID<TV>& grid=fluids_parameters.grid;
    ARRAY<T,VECTOR<int,3> >& phi=fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi;
    T bandwidth=phi_source_bandwidth*grid.max_dx_dy_dz;
    T depth=phi_source_depth*grid.max_dx_dy_dz;
    int count=0;
    for(int object=1;object<=melting_parameters.levelsets.m;object++){
        LEVELSET_TRIANGULATED_OBJECT<T,VECTOR<T,3> >& levelset=*melting_parameters.levelsets(object);
        DEFORMABLE_OBJECT_3D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(melting_parameters.body_index(object));
        PARTICLES<T,VECTOR<T,3> >& particles=deformable_object.particles;
        TRIANGULATED_SURFACE<T>& material_surface=deformable_object.triangles_of_material->material_surface;
        PARTICLES<T,VECTOR<T,3> >& material_particles=material_surface.particles;
        ARRAY<VECTOR_2D<T> >& node_locations=levelset.grid.Node_Locations();
        ARRAY<T>& reaction=*melting_parameters.reaction(object);
        // absurd hack: replace velocities with reactions and use Update_Particle_Velocities to interpolate them to the material surface
        ARRAY<VECTOR<T,3> > V_save(particles.array_collection->Size());ARRAY<VECTOR<T,3> >::Exchange_Arrays(V_save,particles.V.array);
        for(int n=1;n<=reaction.m;n++)if(levelset.node_to_particle_mapping(n)){
            particles.V(levelset.node_to_particle_mapping(n)).x=reaction(n);
            if(inflammable_levelset) particles.V(levelset.node_to_particle_mapping(n)).y=inflammable_levelset->Phi(node_locations(n));}
        deformable_object.triangles_of_material->Update_Particle_Velocities();
        // compute phi contribution
        material_surface.Update_Triangle_List(); // inefficient if positions haven't changed
        if(!material_surface.triangle_hierarchy) material_surface.Initialize_Triangle_Hierarchy();else material_surface.triangle_hierarchy->Update_Boxes();
        ARRAY<bool,VECTOR<int,3> > occupied(grid);deformable_object.collisions.Compute_Occupied_Cells(grid,occupied,false,bandwidth,0);
        for(int i=1;i<=grid.m;i++)for(int j=1;j<=grid.n;j++)for(int ij=1;ij<=grid.mn;ij++)if(occupied(i,j,ij)){
            int triangle;T distance;VECTOR<T,3> X=material_surface.Surface(grid.X(i,j,ij),bandwidth,1e-6,&triangle,&distance);
            if(distance>bandwidth) continue;
            int ti,tj,tk;material_surface.triangle_mesh.triangles.Get(triangle,ti,tj,tk);
            VECTOR<T,3> weights=TRIANGLE_3D<T>::Clamped_Barycentric_Coordinates(X,material_particles.X(ti),material_particles.X(tj),material_particles.X(tk));
            if(inflammable_levelset){
                distance-=min((T)0,weights.x*material_particles.V(ti).y+weights.y*material_particles.V(tj).y+weights.z*material_particles.V(tk).y); // V.y is inflammable phi
                if(distance>bandwidth) continue;}
            T r=weights.x*material_particles.V(ti).x+weights.y*material_particles.V(tj).x+weights.z*material_particles.V(tk).x; // V.x is reaction
            if(r<=reaction_start) continue;
            else if(r<reaction_peak_start && phi_source_smooth_start){phi(i,j,ij)=min(phi(i,j,ij),distance-depth*(r-reaction_start)/(reaction_peak_start-reaction_start));count++;}
            else if(r<reaction_stop){phi(i,j,ij)=min(phi(i,j,ij),distance-depth);count++;}}
        // undo absurd hack
        ARRAY<VECTOR<T,3> >::Exchange_Arrays(V_save,particles.V.array);deformable_object.triangles_of_material->Update_Particle_Velocities();} 
    LOG::cout<<"cloth sourcing adjusted phi at "<<count<<" nodes"<<std::endl; 
}
//#####################################################################
// Function Modify_Particles_And_Phi_For_Melting
//#####################################################################
void Modify_Fluid_For_Melting(const T dt,const T time)
{}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time)
{
    if(fluids_parameters.use_solid_fluid_coupling) SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Set_Fluid_Boundary_Conditions();
    else{
        GRID<TV>& p_grid=fluids_parameters.p_grid,&u_grid=fluids_parameters.u_grid,&v_grid=fluids_parameters.v_grid,&w_grid=fluids_parameters.w_grid;
        for(int object=1;object<=melting_parameters.levelsets.m;object++){
            DEFORMABLE_OBJECT_3D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(melting_parameters.body_index(object));
            TRIANGULATED_SURFACE<T>& material_surface=*deformable_object.collisions.triangulated_surface;
            material_surface.Update_Triangle_List(); // inefficient if positions haven't changed
            if(!material_surface.triangle_hierarchy) material_surface.Initialize_Triangle_Hierarchy();else material_surface.triangle_hierarchy->Update_Boxes();
            ARRAY<bool,VECTOR<int,3> > occupied_cell(p_grid,1);deformable_object.collisions.Compute_Occupied_Cells(p_grid,occupied_cell,false,2*p_grid.max_dx_dy_dz,1);
            // calculate psi_N_u
            for(int i=1;i<=u_grid.m;i++) for(int j=1;j<=u_grid.n;j++) for(int ij=1;ij<=u_grid.mn;ij++) if(occupied_cell(i-1,j,ij) || occupied_cell(i,j,ij)){
                RAY<VECTOR<T,3> > ray(p_grid.X(i-1,j,ij),VECTOR<T,3>(1,0,0),true);ray.t_max=p_grid.dx;ray.semi_infinite=false;
                if(deformable_object.collisions.Triangulated_Surface_Intersection(ray)){
                    fluids_parameters.incompressible.projection.u(i,j,ij)=deformable_object.collisions.Pointwise_Object_Velocity(ray.aggregate_id,ray.Point(ray.t_max)).x;
                    fluids_parameters.incompressible.projection.u_fuel(i,j,ij)=fluids_parameters.incompressible.projection.u(i,j,ij);
                    fluids_parameters.incompressible.projection.elliptic_solver->psi_N_u(i,j,ij)=true;}}
            // calculate psi_N_v
            for(int i=1;i<=v_grid.m;i++) for(int j=1;j<=v_grid.n;j++) for(int ij=1;ij<=v_grid.mn;ij++) if(occupied_cell(i,j-1,ij) || occupied_cell(i,j,ij)){
                RAY<VECTOR<T,3> > ray(p_grid.X(i,j-1,ij),VECTOR<T,3>(0,1,0),true);ray.t_max=p_grid.dy;ray.semi_infinite=false;
                if(deformable_object.collisions.Triangulated_Surface_Intersection(ray)){
                    fluids_parameters.incompressible.projection.v(i,j,ij)=deformable_object.collisions.Pointwise_Object_Velocity(ray.aggregate_id,ray.Point(ray.t_max)).y;
                    fluids_parameters.incompressible.projection.v_fuel(i,j,ij)=fluids_parameters.incompressible.projection.v(i,j,ij);
                    fluids_parameters.incompressible.projection.elliptic_solver->psi_N_v(i,j,ij)=true;}}
            // calculate psi_N_w
            for(int i=1;i<=w_grid.m;i++) for(int j=1;j<=w_grid.n;j++) for(int ij=1;ij<=w_grid.mn;ij++) if(occupied_cell(i,j,ij-1) || occupied_cell(i,j,ij)){
                RAY<VECTOR<T,3> > ray(p_grid.X(i,j,ij-1),VECTOR<T,3>(0,0,1),true);ray.t_max=p_grid.dz;ray.semi_infinite=false;
                if(deformable_object.collisions.Triangulated_Surface_Intersection(ray)){
                    fluids_parameters.incompressible.projection.w(i,j,ij)=deformable_object.collisions.Pointwise_Object_Velocity(ray.aggregate_id,ray.Point(ray.t_max)).z;
                    fluids_parameters.incompressible.projection.w_fuel(i,j,ij)=fluids_parameters.incompressible.projection.w(i,j,ij);
                    fluids_parameters.incompressible.projection.elliptic_solver->psi_N_w(i,j,ij)=true;}}}}
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const
{
    MELTING_EXAMPLE_S3D<T,RW>::Write_Output_Files(frame);
}
//#####################################################################
};
}
#endif

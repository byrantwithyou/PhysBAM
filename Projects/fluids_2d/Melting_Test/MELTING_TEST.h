//#####################################################################
// Copyright 2004, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MELTING_TEST 
//##################################################################### 
#ifndef __MELTING_TEST__
#define __MELTING_TEST__

#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include "../WATER_MELTING_EXAMPLE_2D.h"
#include <Level_Sets/LEVELSET_TRIANGULATED_OBJECT.h>
#include <PhysBAM_Geometry/Red_Green/RED_GREEN_GRID_2D.h>
namespace PhysBAM{

template<class T,class RW=T>
class MELTING_TEST:public WATER_MELTING_EXAMPLE_2D<T,RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::solids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::use_collision_aware_velocity_extrapolation;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::thin_shells_semi_lagrangian_phi;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::thin_shells_semi_lagrangian_velocity;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::use_collision_aware_signed_distance;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::verbose_dt;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::write_time;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::data_directory;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::clamp_phi_with_collision_bodies;
    using MELTING_EXAMPLE_2D<T,RW>::melting_parameters;using MELTING_EXAMPLE_2D<T,RW>::initialized;

    // fluid parameters
    bool use_source;
    BOX_2D<T> source;
    MATRIX<T,3> world_to_source;
    VECTOR_2D<T> source_velocity;
    T initial_water_level;
    bool delete_positive_particles_crossing_bodies;

    // solid parameters
    ARRAY<CIRCLE<T> > circle;
    ARRAY<VECTOR_2D<T> >  initial_position;
    ARRAY<T> initial_orientation;
    ARRAY<VECTOR_2D<T> > initial_velocity;
    ARRAY<T> initial_angular_velocity;

    T side_length;
    int m,n;
    GRID<TV> base_grid;

    int example_number;
    int number_of_objects;

    MELTING_TEST()
        :side_length(3),m(1),n(1)
    {
        example_number=1;

        // Two way coupling
        fluids_parameters.use_solid_fluid_coupling=true;
        solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies=true;
        solids_parameters.deformable_body_parameters.write=true;
        use_collision_aware_velocity_extrapolation=true;
        use_collision_aware_signed_distance=true;
        clamp_phi_with_collision_bodies=true;

        if(fluids_parameters.use_solid_fluid_coupling){
            fluids_parameters.phi_advection=&thin_shells_semi_lagrangian_phi;
            fluids_parameters.vector_advection=&thin_shells_semi_lagrangian_velocity;
            fluids_parameters.levelset_substeps=0.95;
            delete_positive_particles_crossing_bodies=false;}

        output_directory="Melting_Test/output";
        first_frame=0;last_frame=2000;frame_rate=24;
        restart=false;restart_frame=0;

        // Fluids parameters
        fluids_parameters.simulate=true;
        fluids_parameters.grid.Initialize(51,51,3,10,0,7);
        fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=false;
        fluids_parameters.reseeding_frame_rate=10;
        fluids_parameters.bias_towards_negative_particles=true;fluids_parameters.number_particles_per_cell=16;
        fluids_parameters.use_removed_positive_particles=false;fluids_parameters.use_removed_negative_particles=true;
        
        initial_water_level=3;
        use_source=false;
        source=BOX_2D<T>(-0.6,0.6,-1.2,0);
        MATRIX<T,3> source_to_world=MATRIX<T,3>::Translation_Matrix(VECTOR_2D<T>(6.5,6));
        world_to_source=source_to_world.Inverse();
        source_velocity=VECTOR_2D<T>(0,-1);
        
        fluids_parameters.cfl=0.9;
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;fluids_parameters.variable_viscosity=false;fluids_parameters.second_order_pressure=false;
        fluids_parameters.incompressible_iterations=100;
                
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;
        fluids_parameters.write_particles=true;fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=true;fluids_parameters.write_debug_data=true;
        fluids_parameters.store_particle_ids=true;
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.enforce_divergence_free_extrapolation=false;

        // Solids parameters
        solids_parameters.gravity=4.8;
        solids_parameters.cfl=(T)1;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        verbose_dt=true;
        //solids_pre_roll=30;

        // Solids for melting parameters
        solids_parameters.collide_with_interior=false;

        // Melting parameters
        melting_parameters.maximum_depth=7;
        melting_parameters.refine_near_interface=true;
        melting_parameters.refine_for_high_deformation=true;

        if(example_number==1) number_of_objects=1;
        if(example_number==2) number_of_objects=2;        

        circle.Resize(number_of_objects);
        initial_position.Resize(number_of_objects);
        initial_orientation.Resize(number_of_objects);
        initial_velocity.Resize(number_of_objects);
        initial_angular_velocity.Resize(number_of_objects);
        if(example_number==1){
            initial_position(1)=VECTOR_2D<T>(0,0);
            initial_orientation(1)=0;
            initial_velocity(1)=VECTOR_2D<T>(0,0);
            initial_angular_velocity(1)=0;}
        else if(example_number==2){
            initial_position(1)=VECTOR_2D<T>(4.7,3);
            initial_orientation(1)=0;
            initial_velocity(1)=VECTOR_2D<T>(2,0);
            initial_angular_velocity(1)=0;            
            initial_position(2)=VECTOR_2D<T>(6,3);
            initial_orientation(2)=0;
            initial_velocity(2)=VECTOR_2D<T>(-1,0);
            initial_angular_velocity(2)=0;}

        fluids_parameters.Initialize_Domain_Boundary_Conditions(); // sets up the proper wall states
    }

    ~MELTING_TEST() 
    {}

//#####################################################################
// Function Initialize_Deformable_And_Rigid_Bodies
//#####################################################################
void Initialize_Deformable_And_Rigid_Bodies()
{
    for(int i=0;i<number_of_objects;i++){
        int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Embedded_Triangulated_Area();
        Add_Melting_Object(melting_parameters.DEFORMABLE,index);}

    int index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>("../../Public_Data/Rigid_Bodies_2D/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->is_static=true;

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Initialize_Grid
//#####################################################################
void Initialize_Grid(const int object,RED_GREEN_GRID_2D<T>& grid)
{
    grid.Initialize(GRID_2D<T>(4*m+1,4*n+1,0,side_length,0,side_length*n/m),melting_parameters.maximum_depth);
    circle(object).center=VECTOR_2D<T>(grid.uniform_grid.Domain().Center());
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi(const int object,ARRAY<T>& phi)
{
    RED_GREEN_GRID_2D<T>& grid=melting_parameters.levelsets(object)->grid;
    ARRAY<VECTOR_2D<T> >& node_locations=grid.Node_Locations();
    
    for(int p=0;p<phi.m;p++) phi(p)=circle(object).Signed_Distance(node_locations(p));    
}
//#####################################################################
// Function Initialize_Levelset_Velocity
//#####################################################################
void Initialize_Levelset_Velocity(const int object,ARRAY<VECTOR_2D<T> >& V)
{
    RED_GREEN_GRID_2D<T>& grid=melting_parameters.levelsets(1)->grid;
    //ARRAY<VECTOR_2D<T> >& node_locations=grid.Node_Locations();
    
    VECTOR_2D<T> circle_center(grid.uniform_grid.xmin+(T).5*side_length,0);
    for(int p=0;p<V.m;p++) V(p)=VECTOR_2D<T>(0,0);    
//    for(int p=0;p<V.m;p++) V(p)=VECTOR_2D<T>(node_locations(p).y-circle.center.y,0);
//    for(int p=0;p<V.m;p++) V(p)=-node_locations(p)+circle.center;
}
//#####################################################################
// Function Initialize_Particle_Positions_And_Velocities
//#####################################################################
void Initialize_Particle_Positions_And_Velocities(const int object)
{
    TRIANGULATED_AREA<T>& triangulated_area=solids_parameters.deformable_body_parameters.list(object).embedded_triangulated_area->triangulated_area;
    DEFORMABLE_PARTICLES<T,VECTOR_2D<T> >& particles=triangulated_area.particles;

    particles.Update_Velocity();
    triangulated_area.Update_Bounding_Box();
    VECTOR_2D<T> center(triangulated_area.bounding_box->Center());
    for(int i=0;i<particles.array_collection->Size();i++){
        particles.X(i)=center+MATRIX<T,2>::Rotation_Matrix(initial_orientation(object))*(particles.X(i)-center);
        VECTOR_2D<T> radial = particles.X(i)-center;
        particles.V(i)=initial_velocity(object)+initial_angular_velocity(object)*radial.Rotate_Counterclockwise_90();
        particles.X(i)+=initial_position(object);}
}
//#####################################################################
// Function Initialize_Forces
//#####################################################################
void Initialize_Forces()
{
    for(int i=0;i<number_of_objects;i++){
        DEFORMABLE_OBJECT_2D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(i);
        TRIANGULATED_AREA<T>& triangulated_area=deformable_object.embedded_triangulated_area->triangulated_area;
        
        triangulated_area.particles.Store_Mass();
        triangulated_area.Set_Density(.1);
        triangulated_area.Set_Mass_Of_Particles(false);
        
        deformable_object.Delete_Forces();
        solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_object.particles,solids_parameters.gravity,solids_parameters.gravity_direction));
        solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(triangulated_area,new DIAGONALIZED_SPLINE_MODEL_2D<T>(triangulated_area,(T)5e1,(T).3,(T).5,7,(T).01)));

        deformable_object.embedded_triangulated_area->embedded_curve.Set_Density(.1);
        
        //solids_parameters.collision_body_list.Add_Body(&deformable_object.collisions);
        Add_To_Fluid_Simulation(deformable_object,true,true,1);}
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    // Not so good to set up a heaviside function here because then the interface will
    // be exactly between the two nodes which can lead to roundoff issues when setting dirichlet cells, etc.
    GRID<TV>& grid=fluids_parameters.grid;
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++)
        fluids_parameters.particle_levelset_evolution.phi(i,j)=grid.y(j)-grid.ymin-initial_water_level;
}




//#####################################################################
// Function Extrapolate_Phi_Into_Objects
//#####################################################################
void Extrapolate_Phi_Into_Objects(const T time)
{
    if(!initialized) return;
    
    GRID<TV>& grid=fluids_parameters.grid;

    for(int object=0;object<melting_parameters.levelsets.m;object++){
        solids_parameters.deformable_body_parameters.list(object).triangles_of_material->material_area.Initialize_Triangle_Hierarchy();
        solids_parameters.deformable_body_parameters.list(object).triangles_of_material->material_area.Update_Bounding_Box();
        for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) if(solids_parameters.deformable_body_parameters.list(object).triangles_of_material->material_area.Inside(grid.X(i,j)))
            fluids_parameters.particle_levelset_evolution.phi(i,j)=grid.max_dx_dy;}
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
void Adjust_Phi_With_Sources(const T time)
{
    if(!use_source) return;
    GRID<TV>& grid=fluids_parameters.grid;
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++){
        VECTOR_2D<T> source_X=world_to_source*grid.X(i,j);
        if(source.Lazy_Inside(source_X)) fluids_parameters.particle_levelset_evolution.phi(i,j)=min(fluids_parameters.particle_levelset_evolution.phi(i,j),source.Signed_Distance(source_X));}
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,2> >*& cell_centered_mask,const T time)
{
    if(!use_source) return;
    GRID<TV>& grid=fluids_parameters.grid;
    if(cell_centered_mask) delete cell_centered_mask;cell_centered_mask=new ARRAY<bool,VECTOR<int,2> >(grid);
    T padding=3*grid.max_dx_dy;
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) if(!source.Outside(world_to_source*grid.X(i,j),padding)) (*cell_centered_mask)(i,j)=true;
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time)
{
    if(!use_source) return;
    GRID<TV> &u_grid=fluids_parameters.u_grid,&v_grid=fluids_parameters.v_grid;
    PROJECTION_2D<T>& projection=fluids_parameters.incompressible.projection;
    for(int i=0;i<u_grid.m;i++) for(int j=0;j<u_grid.n;j++) if(source.Lazy_Inside(world_to_source*u_grid.X(i,j))){projection.elliptic_solver->psi_N_u(i,j)=true;projection.u(i,j)=source_velocity.x;}
    for(int i=0;i<v_grid.m;i++) for(int j=0;j<v_grid.n;j++) if(source.Lazy_Inside(world_to_source*v_grid.X(i,j))){projection.elliptic_solver->psi_N_v(i,j)=true;projection.v(i,j)=source_velocity.y;}
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time)
{
    if(fluids_parameters.use_solid_fluid_coupling) SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Set_Fluid_Boundary_Conditions();
}
//#####################################################################
// Function Adjust_Particle_For_Objects
//#####################################################################
bool Adjust_Particle_For_Objects(PARTICLE_LEVELSET_PARTICLES<T,VECTOR_2D<T> >& particles,const int index,VECTOR_2D<T>& V,const typename PARTICLE_LEVELSET<T,VECTOR_2D<T> >::PARTICLE_TYPE particle_type,const T dt,const T time)
{
    if(fluids_parameters.use_solid_fluid_coupling)
        return fluids_parameters.Adjust_Particle_For_Objects(fluids_parameters.collision_bodies_affecting_fluid,fluids_parameters.particle_levelset_evolution.particle_levelset.Particle_Collision_Distance(particles.quantized_collision_distance(index)),particles,index,V,particle_type,dt,time);
    else
        return true;
}
//#####################################################################
};    
}
#endif  



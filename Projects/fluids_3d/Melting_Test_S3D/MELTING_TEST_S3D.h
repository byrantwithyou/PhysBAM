//#####################################################################
// Copyright 2005, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MELTING_TEST_S3D
//##################################################################### 
#ifndef __MELTING_TEST_S3D__
#define __MELTING_TEST_S3D__

#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include "../../solids_2d/Embedded_Circle/ANALYTIC_CIRCLE.h"
#include "../FIRE_MELTING_EXAMPLE_S3D.h"
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>

namespace PhysBAM{

template<class T,class RW=T>
class MELTING_TEST_S3D:public FIRE_MELTING_EXAMPLE_S3D<T,RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::solids_parameters;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::thin_shells_semi_lagrangian_density;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::thin_shells_semi_lagrangian_temperature;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::thin_shells_semi_lagrangian_velocity;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::thin_shells_semi_lagrangian_phi;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::first_frame;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::thin_shells_semi_lagrangian_velocity_fuel;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_output_files;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::verbose_dt;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::abort_when_dt_below;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::data_directory;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::use_collision_aware_signed_distance;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_time;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_frame_title;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_substeps;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::clamp_phi_with_collision_bodies;
    using MELTING_EXAMPLE_S3D<T,RW>::melting_parameters;using MELTING_EXAMPLE_S3D<T,RW>::Update_Solids_Topology_For_Melting;

    enum SOURCE_TYPE {BOX_SOURCE,CYLINDER_SOURCE};
    int example_number;
    T rho,rho_bottom,rho_top,buoyancy_constant;
    bool use_source;
    SOURCE_TYPE source_type;
    BOX_3D<T> box_source;
    CYLINDER<T> cylinder_source;
    MATRIX<T,4> world_to_source;
    VECTOR<T,3> source_velocity;
    ARRAY<ARRAY<int> > deformable_object_enslaved_nodes;
    PARAMETER_LIST parameter_list;

    // solid parameters
    ANALYTIC_CIRCLE<T> circle;
    VECTOR<T,3>  initial_position;
    QUATERNION<T> initial_orientation;
    VECTOR<T,3> initial_velocity;
    VECTOR<T,3> initial_angular_velocity;
    T side_length;
    int m,n,mn;
 
    MELTING_TEST_S3D()
        :use_source(false),source_type(BOX_SOURCE),world_to_source(MATRIX<T,4>::Identity_Matrix()),
         initial_position(.25,.25,.5),initial_orientation(QUATERNION<T>(-(T).5*pi,VECTOR<T,3>(1,0,0))),initial_velocity(),initial_angular_velocity(),side_length(.5),m(5),n(5)
    {
        // Two way coupling
        fluids_parameters.use_solid_fluid_coupling=true;
        use_collision_aware_signed_distance=true;
        clamp_phi_with_collision_bodies=true;
        thin_shells_semi_lagrangian_phi.phi_interpolation.compute_leaving=false;
        fluids_parameters.density_container.Set_Custom_Advection(thin_shells_semi_lagrangian_density);
        fluids_parameters.temperature_container.Set_Custom_Advection(thin_shells_semi_lagrangian_temperature);
        fluids_parameters.phi_advection=&thin_shells_semi_lagrangian_phi;
        fluids_parameters.vector_advection=&thin_shells_semi_lagrangian_velocity;
        fluids_parameters.vector_advection_fuel=&thin_shells_semi_lagrangian_velocity_fuel;
        fluids_parameters.levelset_substeps=0.95;

        // Common parameters
        output_directory="Melting_Test_S3D/output";
        first_frame=0;last_frame=500;frame_rate=24;
        restart=false;restart_frame=0;
        write_substeps=true;write_frame_title=true;

        // Fluids parameters 
        fluids_parameters.grid.Initialize(21,31,21,T(0),T(1),0,(T)1.5,T(0),T(1));
        
        fluids_parameters.gravity=0;fluids_parameters.cfl=0.9;
        fluids_parameters.domain_walls[0][0]=fluids_parameters.domain_walls[0][1]=fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[2][0]=fluids_parameters.domain_walls[2][1]=false;
        fluids_parameters.domain_walls[1][0]=true;
        fluids_parameters.kolmogorov=(T)0;
        fluids_parameters.implicit_viscosity=false;
        fluids_parameters.incompressible_iterations=100;

        fluids_parameters.normal_flame_speed=T(.1);
        fluids_parameters.buoyancy_constant=T(.001);
        fluids_parameters.use_vorticity_confinement_fuel=fluids_parameters.use_vorticity_confinement=true;
        fluids_parameters.confinement_parameter_fuel=fluids_parameters.confinement_parameter=(T).6;
    
        // density
        fluids_parameters.density_container.Set_Ambient_Density(0);
        fluids_parameters.density_fuel=1;
        fluids_parameters.density=T(0.1);

        // temperature
        fluids_parameters.temperature_container.Set_Ambient_Temperature((T)283.15);
        fluids_parameters.temperature_container.Set_Cooling_Constant((T)4000);
        fluids_parameters.temperature_products=3000;fluids_parameters.temperature_fuel=298;

        // Solids parameters
        solids_parameters.gravity=(T)9.8;
        solids_parameters.cfl=(T)0.9;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies=true;
        solids_parameters.deformable_body_parameters.write=true;
        solids_parameters.perform_self_collision=false;
        verbose_dt=false;

        // Melting parameters
        melting_parameters.maximum_depth=4;
        melting_parameters.refine_near_interface=true;
        melting_parameters.refine_for_high_deformation=true;

        fluids_parameters.domain_walls[0][0]=fluids_parameters.domain_walls[0][1]=fluids_parameters.domain_walls[2][0]=fluids_parameters.domain_walls[2][1]=fluids_parameters.domain_walls[1][0]=true;
        fluids_parameters.domain_walls[1][1]=false;
        solids_parameters.cfl=(T)3;
        use_source=true;source_type=CYLINDER_SOURCE;cylinder_source=CYLINDER<T>(VECTOR<T,3>((T).5,(T)-.25,(T).5),VECTOR<T,3>((T).5,(T).25,(T).5),(T).15);source_velocity=VECTOR<T,3>(0,1,0);

        // Debugging
        abort_when_dt_below=1e-7;
        write_time=true;
        write_frame_title=true;
        fluids_parameters.write_debug_data=true;
        fluids_parameters.write_thin_shells_advection_cache=false;

        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Parameters_From_Parameter_List(*this,parameter_list);
    }

    ~MELTING_TEST_S3D()
    {}

//#####################################################################
// Function Initialize_Deformable_And_Rigid_Bodies
//#####################################################################
void Initialize_Deformable_And_Rigid_Bodies()
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Embedded_Triangulated_Surface();
    Add_Melting_Object(melting_parameters.DEFORMABLE,index);

    index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>("../../Public_Data/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->is_static=true;

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Initialize_Particle_Positions_And_Velocities
//#####################################################################
void Initialize_Grid(const int object,RED_GREEN_GRID_2D<T>& grid)
{
    assert(object==1);
    grid.Initialize(GRID_2D<T>(4*m+1,4*n+1,0,side_length,0,side_length*n/m),melting_parameters.maximum_depth);
    circle.center=grid.uniform_grid.Domain().Center();circle.radius=.2;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi(const int object,ARRAY<T>& phi)
{
    assert(object==1);
    RED_GREEN_GRID_2D<T>& grid=melting_parameters.levelsets(1)->grid;
    ARRAY<VECTOR_2D<T> >& node_locations=grid.Node_Locations();
    
    for(int p=0;p<phi.m;p++) phi(p)=circle.Phi(node_locations(p));    
}
//#####################################################################
// Function Initialize_Levelset_Velocity
//#####################################################################
void Initialize_Levelset_Velocity(const int object,ARRAY<VECTOR_2D<T> >& V)
{
    assert(object==1);
    //RED_GREEN_GRID_2D<T>& grid=melting_parameters.levelsets(1)->grid;
    //ARRAY<VECTOR_2D<T> >& node_locations=grid.Node_Locations();
    
    for(int p=0;p<V.m;p++) V(p)=VECTOR_2D<T>(0,0);
}
//#####################################################################
// Function Initialize_Particle_Positions_And_Velocities
//#####################################################################
void Initialize_Particle_Positions_And_Velocities(const int object)
{
    TRIANGULATED_SURFACE<T>& triangulated_surface=solids_parameters.deformable_body_parameters.list(object).embedded_triangulated_surface->triangulated_surface;
    DEFORMABLE_PARTICLES<T,VECTOR<T,3> >& particles=solids_parameters.deformable_body_parameters.list(object).particles;

    particles.Update_Velocity();
    triangulated_surface.Update_Bounding_Box();
    VECTOR<T,3> center(triangulated_surface.bounding_box->Center());
    for(int i=0;i<particles.array_collection->Size();i++){
        particles.X(i)=center+initial_orientation.Rotate(particles.X(i)-center);
        particles.V(i)=initial_velocity+VECTOR<T,3>::Cross_Product(initial_angular_velocity,particles.X(i)-center);
        particles.X(i)+=initial_position;}
}
//#####################################################################
// Function Set_Parameter
//#####################################################################
template<class T2> void Set_Parameter(T2& parameter,const std::string& name)
{
    if(parameter_list.Is_Defined(name)){
        parameter=parameter_list.template Get_Parameter<T2>(name);
        std::cout << "[param] set " << name << " = " << parameter << std::endl;}
}
//#####################################################################
// Function Initialize_Forces
//#####################################################################
void Initialize_Forces()
{
    TRIANGULATED_SURFACE<T>& triangulated_surface=solids_parameters.deformable_body_parameters.list(1).embedded_triangulated_surface->triangulated_surface;

    solids_parameters.deformable_body_parameters.list(1).particles.Store_Mass();
    triangulated_surface.Set_Density(1);
    triangulated_surface.Set_Mass_Of_Particles(false);

    solids_parameters.deformable_body_parameters.list(1).Delete_Forces();
    solids_parameters.deformable_body_parameters.list(1).Add_Body_Forces(triangulated_surface,solids_parameters.gravity,solids_parameters.gravity_direction);
    solids_parameters.deformable_body_parameters.list(1).Add_Diagonalized_Linear_Finite_Volume(triangulated_surface,(T)500,(T).3,(T).01);
    solids_parameters.deformable_body_parameters.list(1).Add_Bending_Elements(triangulated_surface.triangle_mesh,(T)1,(T).01);
    Add_To_Fluid_Simulation(solids_parameters.deformable_body_parameters.list(1),true,true,1);
}

//#####################################################################
// Function Solids_Example_Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY<VECTOR<T,3> >& V,const T time)
{
//    Zero_Out_Enslaved_Velocity_Nodes(V,time,id_number);
}
//#####################################################################
// Function Solids_Example_Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR<T,3> >& V,const T time)
{
//    assert(id_number<=deformable_object_enslaved_nodes.m);
//    for(int i=0;i<deformable_object_enslaved_nodes(id_number).m;i++) V(deformable_object_enslaved_nodes(id_number)(i))=VECTOR<T,3>();
}
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
void Update_Time_Varying_Material_Properties(const T time)
{
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time)
{
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Update_Fluid_Parameters(dt,time);
    //if(example_number==4) use_source=time<5;
}
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
template<class SOURCE> void Adjust_Density_And_Temperature_With_Sources(const SOURCE& source,const T time)
{
    for(int i=0;i<fluids_parameters.grid.m;i++) for(int j=0;j<fluids_parameters.grid.n;j++) for(int ij=0;ij<fluids_parameters.grid.mn;ij++)
        if(source.Lazy_Inside(world_to_source*fluids_parameters.grid.X(i,j,ij))){
            if(fluids_parameters.use_density) fluids_parameters.density_container.density_3d(i,j,ij)=rho;
            if(fluids_parameters.use_temperature) fluids_parameters.temperature_container.temperature_3d(i,j,ij)=fluids_parameters.temperature_products;}
}
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Adjust_Density_And_Temperature_With_Sources(const T time)
{
    if(fluids_parameters.fire) return;
    if(use_source){
        if(source_type==BOX_SOURCE) Adjust_Density_And_Temperature_With_Sources(box_source,time);
        else if(source_type==CYLINDER_SOURCE) Adjust_Density_And_Temperature_With_Sources(cylinder_source,time);}
    // keep density >= 0 and T >=ambient temperature
    for(int i=0;i<fluids_parameters.grid.m;i++)for(int j=0;j<fluids_parameters.grid.n;j++)  for(int ij=0;ij<fluids_parameters.grid.mn;ij++){
        if(fluids_parameters.use_density) fluids_parameters.density_container.density_3d(i,j,ij)=max((T)0,fluids_parameters.density_container.density_3d(i,j,ij));
        if(fluids_parameters.use_temperature) fluids_parameters.temperature_container.temperature_3d(i,j,ij)=max((T)fluids_parameters.temperature_container.ambient_temperature,fluids_parameters.temperature_container.temperature_3d(i,j,ij));}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
template<class SOURCE> void Get_Source_Velocities(const SOURCE& source,const T time)
{
    GRID<TV> &u_grid=fluids_parameters.u_grid,&v_grid=fluids_parameters.v_grid,&w_grid=fluids_parameters.w_grid;
    PROJECTION_3D<T>& projection=fluids_parameters.incompressible.projection;
    for(int i=0;i<u_grid.m;i++) for(int j=0;j<u_grid.n;j++) for(int ij=0;ij<u_grid.mn;ij++) if(source.Lazy_Inside(world_to_source*u_grid.X(i,j,ij))){
        projection.elliptic_solver->psi_N_u(i,j,ij)=true;projection.u(i,j,ij)=source_velocity.x;if(fluids_parameters.fire)projection.u_fuel(i,j,ij)=source_velocity.x;}
    for(int i=0;i<v_grid.m;i++) for(int j=0;j<v_grid.n;j++) for(int ij=0;ij<v_grid.mn;ij++) if(source.Lazy_Inside(world_to_source*v_grid.X(i,j,ij))){
        projection.elliptic_solver->psi_N_v(i,j,ij)=true;projection.v(i,j,ij)=source_velocity.y;if(fluids_parameters.fire)projection.v_fuel(i,j,ij)=source_velocity.y;}
    for(int i=0;i<w_grid.m;i++) for(int j=0;j<w_grid.n;j++) for(int ij=0;ij<w_grid.mn;ij++) if(source.Lazy_Inside(world_to_source*w_grid.X(i,j,ij))){
        projection.elliptic_solver->psi_N_w(i,j,ij)=true;projection.w(i,j,ij)=source_velocity.z;if(fluids_parameters.fire)projection.w_fuel(i,j,ij)=source_velocity.z;}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time)
{
    if(!use_source) return;
    else if(source_type==BOX_SOURCE) Get_Source_Velocities(box_source,time);
    else if(source_type==CYLINDER_SOURCE) Get_Source_Velocities(cylinder_source,time);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time)
{
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Set_Fluid_Boundary_Conditions();
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
template<class SOURCE> void Initialize_Phi(const SOURCE& source)
{
    for(int i=0;i<fluids_parameters.grid.m;i++) for(int j=0;j<fluids_parameters.grid.n;j++) for(int ij=0;ij<fluids_parameters.grid.mn;ij++) 
        if(source.Lazy_Inside(world_to_source*fluids_parameters.grid.X(i,j,ij)))
            fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j,ij)=fluids_parameters.grid.dx;
        else
            fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j,ij)=-fluids_parameters.grid.dx;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    if(!use_source) return;
    else if(source_type==BOX_SOURCE) Initialize_Phi(box_source);
    else if(source_type==CYLINDER_SOURCE) Initialize_Phi(cylinder_source);
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
template<class SOURCE> void Adjust_Phi_With_Sources(const SOURCE& source,const T time)
{
    for(int i=0;i<fluids_parameters.grid.m;i++) for(int j=0;j<fluids_parameters.grid.n;j++) for(int ij=0;ij<fluids_parameters.grid.mn;ij++) 
//        if(source.Lazy_Inside(world_to_source*fluids_parameters.grid.X(i,j,ij)))
        if(source.Lazy_Inside(world_to_source*fluids_parameters.grid.X(i,j,ij))&&fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j,ij)<0)
            fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j,ij)=fluids_parameters.grid.dx;
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
void Adjust_Phi_With_Sources(const T time)
{
    if(!use_source) return;
    else if(source_type==BOX_SOURCE) Adjust_Phi_With_Sources(box_source,time);
    else if(source_type==CYLINDER_SOURCE) Adjust_Phi_With_Sources(cylinder_source,time);
}
//#####################################################################
};    
}
#endif

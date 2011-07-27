//#####################################################################
// Copyright 2005, Tamar Shinar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNDERWATER_DOME 
//##################################################################### 
#ifndef __UNDERWATER_DOME__
#define __UNDERWATER_DOME__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/FLOOD_FILL_3D.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/SIGNED_DISTANCE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>
namespace PhysBAM{

template<class T,class RW=T>
class UNDERWATER_DOME:public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::solids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::verbose_dt;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_time;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_frame_title;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::data_directory;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::abort_when_dt_below;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Use_Thin_Shells_Fluid_Coupling_Defaults;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::phi_collidable_interpolation;

    SPHERE<T> dome_sphere;
    T initial_water_level;
    T initial_trapped_air_pressure;
    int initial_trapped_air_cells; 
    bool initialized_air_cells;
    int last_number_of_air_regions;
    T last_trapped_air_pressure;
    ARRAY<T,VECTOR<int,3> > last_region_colors;

    UNDERWATER_DOME(T pressure_input=1.0)
        :SOLIDS_FLUIDS_EXAMPLE_3D<RW>(fluids_parameters.WATER),initial_water_level((T).5),initial_trapped_air_pressure(pressure_input),initialized_air_cells(false),last_number_of_air_regions(2),last_trapped_air_pressure(pressure_input)

    {
        fluids_parameters.use_solid_fluid_coupling=true;
        Use_Thin_Shells_Fluid_Coupling_Defaults();
        first_frame=0;last_frame=1000;frame_rate=60;
        restart=false;restart_frame=0;
        output_directory=STRING_UTILITIES::string_sprintf("output__pressure_%.2f",initial_trapped_air_pressure);

        // Fluids parameters
        fluids_parameters.grid.Initialize(51,51,51,0,1,0,1,0,1);
        fluids_parameters.number_particles_per_cell=8;
        fluids_parameters.reseeding_frame_rate=10;
        fluids_parameters.incompressible_iterations=200;
        fluids_parameters.cfl=0.9;
        fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][2]=fluids_parameters.domain_walls[2][1]=fluids_parameters.domain_walls[3][1]=fluids_parameters.domain_walls[3][2]=true;fluids_parameters.domain_walls[2][2]=false;
        fluids_parameters.bias_towards_negative_particles=true;
        fluids_parameters.use_removed_positive_particles=false;fluids_parameters.use_removed_negative_particles=true;
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;fluids_parameters.variable_viscosity=false;fluids_parameters.second_order_cut_cell_method=false;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;
        fluids_parameters.write_particles=true;fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=true;
        fluids_parameters.store_particle_ids=true;
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.enforce_divergence_free_extrapolation=false;
        fluids_parameters.use_density=true;

        // Solids parameters
        solids_parameters.gravity=9.8;
        solids_parameters.cfl=(T).9;
        solids_parameters.perform_self_collision=true;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies=true;
        solids_parameters.deformable_body_parameters.write=true;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=false;solids_parameters.rigid_body_parameters.spatial_partition_based_on_scene_size=true;

        // Debugging
        fluids_parameters.write_debug_data=true;
        abort_when_dt_below=1e-7;
        write_time=true;
        write_frame_title=true;
        fluids_parameters.simulate=true;
        fluids_parameters.write_ghost_values=false;
        verbose_dt=true;
    }

    ~UNDERWATER_DOME() 
    {}

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(true,data_directory+"/Rigid_Bodies/Thin_Shells/sphere_with_hole",(T).4,true,false,false,false);
    RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list(id);
    rigid_body.frame.r=QUATERNION<T>(2*pi/5,VECTOR<T,3>(0,0,1));
    rigid_body.frame.t=VECTOR<T,3>(.5,-.1,.5);
    solids_parameters.rigid_body_parameters.list(id)->is_static=true;solids_parameters.rigid_body_parameters.list(id)->is_kinematic=false;
    Add_To_Fluid_Simulation(*solids_parameters.rigid_body_parameters.list(id),true,false,100);
    dome_sphere=SPHERE<T>(rigid_body.frame.t,.4);

    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Initialize_Bodies();
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    GRID<TV>& grid=fluids_parameters.grid;
    ARRAY<T,VECTOR<int,3> >& phi=fluids_parameters.particle_levelset_evolution.phi;

    SIGNED_DISTANCE::Calculate(dome_sphere,grid,phi);phi*=-1;
    for(int i=1;i<=grid.m;i++) for(int j=1;j<=grid.n;j++) for(int ij=1;ij<=grid.mn;ij++) phi(i,j,ij)=max(phi(i,j,ij),grid.y(j)-initial_water_level);
}
//#####################################################################
// Function Adjust_Particle_For_Objects
//#####################################################################
bool Adjust_Particle_For_Objects(PARTICLE_LEVELSET_PARTICLES<T,VECTOR<T,3> >& particles,const int index,VECTOR<T,3>& V,const typename PARTICLE_LEVELSET<T,VECTOR<T,3> >::PARTICLE_TYPE particle_type,const T dt,const T time)
{
    if(fluids_parameters.use_solid_fluid_coupling) 
        return fluids_parameters.Adjust_Particle_For_Objects(fluids_parameters.collision_bodies_affecting_fluid,fluids_parameters.particle_levelset_evolution.particle_levelset.Particle_Collision_Distance(particles.quantized_collision_distance(index)),particles,index,V,particle_type,dt,time);
    else return true;
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time)
{
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Set_Fluid_Boundary_Conditions();
}
//#####################################################################
// Function Set_Dirichlet_Nodes
//#####################################################################
void Set_Dirichlet_Nodes()
{
    PROJECTION_3D<T>& projection=fluids_parameters.incompressible.projection;
    ARRAY<T,VECTOR<int,3> >& phi=fluids_parameters.particle_levelset_evolution.phi;

    // set the dirichlet nodes
    for(int i=1;i<=projection.p_grid.m;i++) for(int j=1;j<=projection.p_grid.n;j++) for(int ij=1;ij<=projection.p_grid.mn;ij++){
        if(phi_collidable_interpolation.From_Base_Node(fluids_parameters.grid,phi,projection.p_grid.X(i,j,ij),i,j,ij) > 0)
            projection.elliptic_solver->psi_D(i,j,ij)=true;}
}
//#####################################################################
// Function Find_Air_Regions
//#####################################################################
int Find_Air_Regions(ARRAY<int,VECTOR<int,3> >& filled_region_colors,int& trapped_air_cells)
{
    PROJECTION_3D<T>& projection=fluids_parameters.incompressible.projection;
    int number_of_regions;

    // flood fill the domain to find air regions
    filled_region_colors.Resize(projection.p_grid,1);
    FLOOD_FILL_3D flood_fill;
    for(int i=0;i<=projection.p_grid.m+1;i++)for(int j=0;j<=projection.p_grid.n+1;j++)for(int ij=0;ij<=projection.p_grid.mn+1;ij++){
        if(!projection.elliptic_solver->psi_D(i,j,ij)) filled_region_colors(i,j,ij)=-1;
        else filled_region_colors(i,j,ij)=0;}
    number_of_regions=flood_fill.Flood_Fill(filled_region_colors,projection.elliptic_solver->psi_N_u,projection.elliptic_solver->psi_N_v,projection.elliptic_solver->psi_N_w);

    // count trapped air cells
    int ambient_air_color=filled_region_colors(1,fluids_parameters.grid.n,1); // assume top front corner is ambient air
    trapped_air_cells=0;
    for(int i=1;i<=projection.p_grid.m;i++) for(int j=1;j<=projection.p_grid.n;j++) for(int ij=1;ij<=projection.p_grid.mn;ij++)
        if(projection.elliptic_solver->psi_D(i,j,ij)&&filled_region_colors(i,j,ij)!=ambient_air_color) trapped_air_cells+=1;
    
    return number_of_regions;
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
void Set_Dirichlet_Boundary_Conditions(const T time)
{
    PROJECTION_3D<T>& projection=fluids_parameters.incompressible.projection;

    Set_Dirichlet_Nodes();

    // find air regions
    int trapped_air_cells=0,number_of_regions=-1;
    ARRAY<int,VECTOR<int,3> > filled_region_colors;number_of_regions=Find_Air_Regions(filled_region_colors,trapped_air_cells);

    if(!initialized_air_cells){
        initial_trapped_air_cells=trapped_air_cells;
        last_trapped_air_pressure=initial_trapped_air_pressure;last_number_of_air_regions=number_of_regions;
        initialized_air_cells=true;}

    if(number_of_regions<last_number_of_air_regions){initial_trapped_air_cells=trapped_air_cells;initial_trapped_air_pressure=last_trapped_air_pressure;}

    // set the pressure at dirichlet nodes
    int ambient_air_color=filled_region_colors(1,fluids_parameters.grid.n,1); // assume top front corner is ambient air
    if(trapped_air_cells!=0) last_trapped_air_pressure=initial_trapped_air_pressure*initial_trapped_air_cells/trapped_air_cells;
    for(int i=1;i<=projection.p_grid.m;i++) for(int j=1;j<=projection.p_grid.n;j++) for(int ij=1;ij<=projection.p_grid.mn;ij++){
        if(projection.elliptic_solver->psi_D(i,j,ij)) 
            if(filled_region_colors(i,j,ij)!=ambient_air_color) projection.p(i,j,ij)=last_trapped_air_pressure;
            else projection.p(i,j,ij)=1;}

    last_number_of_air_regions=number_of_regions;
}
//#####################################################################
};    
}
#endif  

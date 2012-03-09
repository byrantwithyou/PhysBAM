//#####################################################################
// Copyright 2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BURNING_PAPER
//##################################################################### 
#ifndef __BURNING_PAPER__
#define __BURNING_PAPER__

#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include "../FIRE_MELTING_EXAMPLE_S3D.h"
#include <Geometry/EMBEDDED_TRIANGULATED_SURFACE.h>
namespace PhysBAM{

template<class T,class RW=T>
class BURNING_PAPER:public FIRE_MELTING_EXAMPLE_S3D<T,RW>
{
public:
    using FIRE_MELTING_EXAMPLE_S3D<T,RW>::fluids_parameters;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::solids_parameters;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::thin_shells_semi_lagrangian_density;
    using FIRE_MELTING_EXAMPLE_S3D<T,RW>::thin_shells_semi_lagrangian_temperature;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::thin_shells_semi_lagrangian_velocity;
    using FIRE_MELTING_EXAMPLE_S3D<T,RW>::thin_shells_semi_lagrangian_phi;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::first_frame;
    using FIRE_MELTING_EXAMPLE_S3D<T,RW>::thin_shells_semi_lagrangian_velocity_fuel;
    using FIRE_MELTING_EXAMPLE_S3D<T,RW>::last_frame;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::frame_rate;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::write_output_files;
    using FIRE_MELTING_EXAMPLE_S3D<T,RW>::output_directory;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::restart;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::restart_frame;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::verbose_dt;
    using FIRE_MELTING_EXAMPLE_S3D<T,RW>::abort_when_dt_below;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::data_directory;
    using FIRE_MELTING_EXAMPLE_S3D<T,RW>::use_collision_aware_signed_distance;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::write_time;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::write_frame_title;
    using FIRE_MELTING_EXAMPLE_S3D<T,RW>::write_substeps;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::clamp_phi_with_collision_bodies;
    using FIRE_MELTING_EXAMPLE_S3D<T,RW>::melting_parameters;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::Update_Solids_Topology_For_Melting;
    using FIRE_MELTING_EXAMPLE_S3D<T,RW>::cloth_temperature_time_constant;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::cloth_ignition_temperature;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::reaction_peak_start;
    using FIRE_MELTING_EXAMPLE_S3D<T,RW>::reaction_peak_stop;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::reaction_stop;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::divergence_source_bandwidth;
    using FIRE_MELTING_EXAMPLE_S3D<T,RW>::divergence_source_rate;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::phi_source_bandwidth;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::phi_source_depth;
    using FIRE_MELTING_EXAMPLE_S3D<T,RW>::inflammable_levelset;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::phi_source_smooth_start;

    enum SOURCE_TYPE {BOX_SOURCE,CYLINDER_SOURCE};
    int example_number;
    bool use_source;
    T source_start_time,source_stop_time;
    SOURCE_TYPE source_type;
    BOX_3D<T> box_source;
    CYLINDER<T> cylinder_source;
    MATRIX<T,4> world_to_source;
    VECTOR<T,3> source_velocity;
    PARAMETER_LIST parameter_list;

    // solid parameters
    VECTOR<T,3>  initial_position;
    QUATERNION<T> initial_orientation;
    VECTOR<T,3> initial_velocity;
    VECTOR<T,3> initial_angular_velocity;
    T side_length;
    int m,n,mn;
    ARRAY<int> constrained_nodes;

    // cutout
    bool cutout;
    GRID<TV> cutout_grid;
    ARRAY<T,VECTOR<int,2> > cutout_phi;
    LEVELSET_2D<T> cutout_levelset;

    T cloth_temperature_time_constant_after;

    BURNING_PAPER(const int example_number_input)
        :example_number(example_number_input),use_source(false),source_start_time(-1),source_type(BOX_SOURCE),world_to_source(MATRIX<T,4>::Identity_Matrix()),
        initial_position(0,.25,.5),initial_orientation(),initial_velocity(),initial_angular_velocity(),side_length(1),m(22),n(17),
        cutout(false),cutout_levelset(cutout_grid,cutout_phi),cloth_temperature_time_constant_after(0)
    {
        initial_orientation=QUATERNION<T>((T).3*pi,VECTOR<T,3>(0,0,1))*QUATERNION<T>(-(T).5*pi,VECTOR<T,3>(1,0,0));

        // Two way coupling
        fluids_parameters.use_solid_fluid_coupling=true;
        if(fluids_parameters.use_solid_fluid_coupling){
            use_collision_aware_signed_distance=true;
            clamp_phi_with_collision_bodies=true;
            fluids_parameters.density_container.Set_Custom_Advection(thin_shells_semi_lagrangian_density);
            fluids_parameters.temperature_container.Set_Custom_Advection(thin_shells_semi_lagrangian_temperature);
            fluids_parameters.phi_advection=&thin_shells_semi_lagrangian_phi;
            fluids_parameters.vector_advection=&thin_shells_semi_lagrangian_velocity;
            fluids_parameters.vector_advection_fuel=&thin_shells_semi_lagrangian_velocity_fuel;}
        fluids_parameters.scalar_substeps=0.95;

        int r=-10;
        //example_number=8;
        LOG::cout<<"example number = "<<example_number<<std::endl;
        if(example_number<=0) exit(1);

        // Common parameters
        output_directory=STRING_UTILITIES::string_sprintf("Burning_Paper/output%d",example_number);
        first_frame=0;last_frame=750;frame_rate=24*2;
        restart=false;restart_frame=40;
        //write_substeps=true;write_frame_title=true;

        fluids_parameters.gravity=0;fluids_parameters.cfl=0.9;
        fluids_parameters.domain_walls[0][0]=fluids_parameters.domain_walls[0][1]=fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[2][0]=fluids_parameters.domain_walls[2][1]=false;
        fluids_parameters.domain_walls[1][0]=true;
        fluids_parameters.kolmogorov=(T)0;
        fluids_parameters.implicit_viscosity=false;
        fluids_parameters.incompressible_iterations=100;
        fluids_parameters.normal_flame_speed=(T).1;
        fluids_parameters.buoyancy_constant=(T).001;
        fluids_parameters.use_vorticity_confinement_fuel=fluids_parameters.use_vorticity_confinement=true;
        fluids_parameters.confinement_parameter_fuel=fluids_parameters.confinement_parameter=(T).05;
        fluids_parameters.write_debug_data=false;
        source_velocity=VECTOR<T,3>(0,2,0);
        fluids_parameters.domain_walls[0][0]=fluids_parameters.domain_walls[0][1]=fluids_parameters.domain_walls[2][0]=fluids_parameters.domain_walls[2][1]=fluids_parameters.domain_walls[1][0]=true;
        fluids_parameters.domain_walls[1][1]=false;
        use_source=true;source_type=CYLINDER_SOURCE;
        cylinder_source=CYLINDER<T>(VECTOR<T,3>((T).5,(T)-.25,(T).5),VECTOR<T,3>((T).5,(T).25,(T).5),(T).15);

        solids_parameters.cfl=(T).5;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        solids_parameters.implicit_solve_parameters.cg_iterations=500;

        if(example_number==13){
            //restart=true;restart_frame=67;
            r=12;
            m*=3;n*=3;
            fluids_parameters.grid.Initialize(10*r+1,16*r+1,10*r+1,-.25,1,0,2,-.125,1.125);
            melting_parameters.maximum_depth=1;
            source_stop_time=1.5;
            cloth_ignition_temperature=291;
            cloth_temperature_time_constant=.015;
            cloth_temperature_time_constant_after=.0075;
            reaction_peak_start=.25;
            reaction_peak_stop=1.25;
            reaction_stop=1.5;
            divergence_source_bandwidth=2.5;
            divergence_source_rate=300;
            phi_source_bandwidth=3;
            phi_source_depth=2;
            PHYSBAM_FATAL_ERROR("Preroll is gone; this example is likely broken.");
            fluids_parameters.normal_flame_speed=(T).1;
            source_velocity=VECTOR<T,3>(0,1,0);
            cylinder_source=CYLINDER<T>(VECTOR<T,3>((T).5,(T)-.25,(T).5),VECTOR<T,3>((T).5,(T).25,(T).5),(T).1);
            fluids_parameters.confinement_parameter_fuel=fluids_parameters.confinement_parameter=(T).2;}
        else if(example_number==16){
            //restart=true;restart_frame=67;
            r=12;
            m*=5;n*=5;
            //r=4;
            //m*=2;n*=2;
            fluids_parameters.grid.Initialize(10*r+1,16*r+1,7*r+1,-.25,1,0,2,.0625,1-.0625);
            Initialize_Siggraph();
            melting_parameters.maximum_depth=1;
            source_start_time=.1;
            source_stop_time=1;
            cloth_ignition_temperature=288;
            cloth_temperature_time_constant=.015;
            cloth_temperature_time_constant_after=.05;
            reaction_peak_start=.25;
            reaction_peak_stop=1.25;
            reaction_stop=1.5;
            divergence_source_bandwidth=2.5;
            divergence_source_rate=300;
            phi_source_bandwidth=3;
            phi_source_depth=2;
            //solids_parameters.cfl=(T).125;
            //solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            PHYSBAM_FATAL_ERROR("Preroll is gone; this example is likely broken.");
            fluids_parameters.normal_flame_speed=(T).1;
            source_velocity=VECTOR<T,3>(0,1,0);
            cylinder_source=CYLINDER<T>(VECTOR<T,3>(0,(T)-.25,(T).5),VECTOR<T,3>(0,(T).05,(T).5),(T).075);
            fluids_parameters.confinement_parameter_fuel=fluids_parameters.confinement_parameter=(T).2;}
        else if(example_number==17){
            //restart=true;restart_frame=67;
            r=12;
            m*=5;n*=5;
            fluids_parameters.grid.Initialize(10*r+1,16*r+1,7*r+1,-.25,1,0,2,.0625,1-.0625);
            Initialize_SCA();
            melting_parameters.maximum_depth=1;
            source_start_time=.1;
            source_stop_time=1;
            cloth_ignition_temperature=288;
            cloth_temperature_time_constant=.015;
            cloth_temperature_time_constant_after=.05;
            reaction_peak_start=.25;
            reaction_peak_stop=1.25;
            reaction_stop=1.5;
            divergence_source_bandwidth=2.5;
            divergence_source_rate=300;
            phi_source_bandwidth=3;
            phi_source_depth=2;
            //solids_parameters.cfl=(T).125;
            //solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            PHYSBAM_FATAL_ERROR("Preroll is gone; this example is likely broken.");
            fluids_parameters.normal_flame_speed=(T).1;
            source_velocity=VECTOR<T,3>(0,1,0);
            cylinder_source=CYLINDER<T>(VECTOR<T,3>(0,(T)-.25,(T).5),VECTOR<T,3>(0,(T).05,(T).5),(T).075);
            fluids_parameters.confinement_parameter_fuel=fluids_parameters.confinement_parameter=(T).2;}
        else if(example_number==18){
            //restart=true;restart_frame=67;
            r=12;
            m*=5;n*=5;
            fluids_parameters.grid.Initialize(10*r+1,16*r+1,7*r+1,-.25,1,0,2,.0625,1-.0625);
            Initialize_SCA();
            melting_parameters.maximum_depth=1;
            source_start_time=.1;
            source_stop_time=1;
            cloth_ignition_temperature=288;
            cloth_temperature_time_constant=.015;
            cloth_temperature_time_constant_after=.05;
            reaction_peak_start=.25;
            reaction_peak_stop=1.25;
            reaction_stop=1.5;
            divergence_source_bandwidth=2.5;
            divergence_source_rate=300;
            phi_source_bandwidth=3;
            phi_source_depth=2;
            phi_source_smooth_start=false;
            //solids_parameters.cfl=(T).125;
            //solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            PHYSBAM_FATAL_ERROR("Preroll is gone; this example is likely broken.");
            fluids_parameters.normal_flame_speed=(T).1;
            source_velocity=VECTOR<T,3>(0,1,0);
            cylinder_source=CYLINDER<T>(VECTOR<T,3>(0,(T)-.25,(T).5),VECTOR<T,3>(0,(T).05,(T).5),(T).075);
            fluids_parameters.confinement_parameter_fuel=fluids_parameters.confinement_parameter=(T).2;}
        else if(example_number==19){
            //restart=true;restart_frame=67;
            r=12;
            m*=5;n*=5;
            fluids_parameters.grid.Initialize(10*r+1,16*r+1,7*r+1,-.25,1,0,2,.0625,1-.0625);
            Initialize_TVCG();
            melting_parameters.maximum_depth=1;
            source_start_time=.1;
            source_stop_time=1;
            cloth_ignition_temperature=288;
            cloth_temperature_time_constant=.015;
            cloth_temperature_time_constant_after=.05;
            reaction_peak_start=.25;
            reaction_peak_stop=1.25;
            reaction_stop=1.5;
            divergence_source_bandwidth=2.5;
            divergence_source_rate=300;
            phi_source_bandwidth=3;
            phi_source_depth=2;
            phi_source_smooth_start=false;
            //solids_parameters.cfl=(T).125;
            //solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            PHYSBAM_FATAL_ERROR("Preroll is gone; this example is likely broken.");
            fluids_parameters.normal_flame_speed=(T).1;
            source_velocity=VECTOR<T,3>(0,1,0);
            cylinder_source=CYLINDER<T>(VECTOR<T,3>(0,(T)-.25,(T).5),VECTOR<T,3>(0,(T).05,(T).5),(T).075);
            fluids_parameters.confinement_parameter_fuel=fluids_parameters.confinement_parameter=(T).2;}
        else if(example_number==20){
            //restart=true;restart_frame=67;
            r=12;
            m*=5;n*=5;
            fluids_parameters.grid.Initialize(10*r+1,16*r+1,7*r+1,-.25,1,0,2,.0625,1-.0625);
            Initialize_TVCG();
            melting_parameters.maximum_depth=1;
            source_start_time=.1;
            source_stop_time=1;
            cloth_ignition_temperature=288;
            cloth_temperature_time_constant=.015;
            cloth_temperature_time_constant_after=.073;
            reaction_peak_start=.25;
            reaction_peak_stop=1.25;
            reaction_stop=1.5;
            divergence_source_bandwidth=2.5;
            divergence_source_rate=300;
            phi_source_bandwidth=3;
            phi_source_depth=2;
            phi_source_smooth_start=false;
            //solids_parameters.cfl=(T).125;
            //solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            PHYSBAM_FATAL_ERROR("Preroll is gone; this example is likely broken.");
            fluids_parameters.normal_flame_speed=(T).1;
            source_velocity=VECTOR<T,3>(0,1,0);
            cylinder_source=CYLINDER<T>(VECTOR<T,3>(0,(T)-.25,(T).5),VECTOR<T,3>(0,(T).05,(T).5),(T).075);
            fluids_parameters.confinement_parameter_fuel=fluids_parameters.confinement_parameter=(T).2;}
        else{LOG::cout<<"Invalid example number "<<example_number<<std::endl;exit(1);}
    
        // density
        fluids_parameters.density_container.Set_Ambient_Density(0);
        fluids_parameters.density_fuel=1;
        fluids_parameters.density=(T).1;

        // temperature
        fluids_parameters.temperature_container.Set_Ambient_Temperature((T)283.15);
        fluids_parameters.temperature_container.Set_Cooling_Constant((T)4000);
        fluids_parameters.temperature_products=3000;fluids_parameters.temperature_fuel=298;

        // Solids parameters
        solids_parameters.gravity=(T)9.8;
        solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies=true;
        solids_parameters.deformable_body_parameters.write=true;
        solids_parameters.perform_self_collision=false;
        solids_parameters.collisions_disable_repulsions_based_on_proximity_factor=1.5;
        solids_parameters.check_initial_mesh_for_self_intersection=false;
        verbose_dt=false;

        // Melting parameters
        melting_parameters.refine_near_interface=true;
        melting_parameters.refine_for_high_deformation=true;

        // Debugging
        abort_when_dt_below=1e-7;
        write_time=true;
        write_frame_title=true;
        fluids_parameters.write_thin_shells_advection_cache=false;

        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Parameters_From_Parameter_List(*this,parameter_list);
    }

    ~BURNING_PAPER()
    {}

//#####################################################################
// Function Initialize_Deformable_And_Rigid_Bodies
//#####################################################################
void Initialize_Deformable_And_Rigid_Bodies()
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Embedded_Triangulated_Surface(15,1e-2);
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
    grid.Initialize(GRID_2D<T>(4*m+1,4*n+1,0,side_length,-side_length*n/m/2,side_length*n/m/2),melting_parameters.maximum_depth);
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi(const int object,ARRAY<T>& phi)
{
    assert(object==1);
    ARRAY<T>::copy(-10,phi);
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
    for(int i=0;i<particles.array_collection->Size();i++){
        particles.X(i)=initial_orientation.Rotate(particles.X(i));
        particles.V(i)=initial_velocity+VECTOR<T,3>::Cross_Product(initial_angular_velocity,particles.X(i));
        particles.X(i)+=initial_position;}
}
//#####################################################################
// Function Set_Parameter
//#####################################################################
template<class T2> void Set_Parameter(T2& parameter,const std::string& name)
{
    if(parameter_list.Is_Defined(name)){
        parameter=parameter_list.template Get_Parameter<T2>(name);
        LOG::cout << "[param] set " << name << " = " << parameter << std::endl;}
}
//#####################################################################
// Function Initialize_Forces
//#####################################################################
void Initialize_Forces()
{
    if(example_number==15 && restart && restart_frame==210){
        solids_parameters.cfl=(T).125;
        solids_parameters.implicit_solve_parameters.cg_tolerance=1e-3;}
    else if(example_number==16 && restart && restart_frame==212){
        solids_parameters.cfl=(T).125;
        solids_parameters.implicit_solve_parameters.cg_tolerance=1e-3;}
    else if(example_number==18 && restart && restart_frame==238){
        solids_parameters.cfl=(T).125;
        solids_parameters.implicit_solve_parameters.cg_tolerance=1e-3;}
    else if(example_number==20 && restart && restart_frame==247){
        solids_parameters.cfl=(T).125;
        solids_parameters.implicit_solve_parameters.cg_tolerance=1e-3;}

    TRIANGULATED_SURFACE<T>& triangulated_surface=solids_parameters.deformable_body_parameters.list(1).embedded_triangulated_surface->triangulated_surface;
    DEFORMABLE_PARTICLES<T,VECTOR<T,3> >& particles=triangulated_surface.particles;

    solids_parameters.deformable_body_parameters.list(1).particles.Store_Mass();
    triangulated_surface.Set_Density(1);
    triangulated_surface.Set_Mass_Of_Particles(false);

    solids_parameters.deformable_body_parameters.list(1).Delete_Forces();
    if(0 && example_number==5){
        solids_parameters.deformable_body_parameters.list(1).Add_Body_Forces(triangulated_surface,solids_parameters.gravity/10,solids_parameters.gravity_direction);
        solids_parameters.deformable_body_parameters.list(1).Add_Diagonalized_Linear_Finite_Volume(triangulated_surface,(T)50,(T).3,(T).1);
        solids_parameters.deformable_body_parameters.list(1).Add_Bending_Elements(triangulated_surface.triangle_mesh,(T).1,(T).01);}
    else{
        solids_parameters.deformable_body_parameters.list(1).Add_Body_Forces(triangulated_surface,solids_parameters.gravity,solids_parameters.gravity_direction);
        solids_parameters.deformable_body_parameters.list(1).Add_Diagonalized_Linear_Finite_Volume(triangulated_surface,(T)500,(T).3,(T).1);
        solids_parameters.deformable_body_parameters.list(1).Add_Bending_Elements(triangulated_surface.triangle_mesh,(T)1,(T).1);}

    if(fluids_parameters.use_solid_fluid_coupling){
        static int count=0;count++;
        if(count==1+restart) Add_To_Fluid_Simulation(solids_parameters.deformable_body_parameters.list(1),true,false,1);}

    BOX_2D<T> box=melting_parameters.levelsets(1)->grid.uniform_grid.Domain();
    constrained_nodes.Resize(0);ARRAY<VECTOR<T,3> > corners(4);
    corners(1)=VECTOR<T,3>(box.xmin,box.ymin,0);corners(2)=VECTOR<T,3>(box.xmin,box.ymax,0);corners(3)=VECTOR<T,3>(box.xmax,box.ymin,0);corners(4)=VECTOR<T,3>(box.xmax,box.ymax,0);
    T lower_edge=box.xmin;
    for(int p=0;p<particles.array_collection->Size();p++)for(int c=0;c<4;c++){
        VECTOR<T,3>& X=particles.X(p);
        if((corners(c)-X).Magnitude()<1e-5 || fabs(X.x-lower_edge)<1e-5){constrained_nodes.Append(p);break;}}
}
//#####################################################################
// Function Melting_Levelset_Substep
//#####################################################################
void Melting_Levelset_Substep(const int object,const T dt,const T time)
{
    FIRE_MELTING_EXAMPLE_S3D<T,RW>::Melting_Levelset_Substep(object,dt,time);

    if(0 && cutout){
        LEVELSET_TRIANGULATED_OBJECT<T,VECTOR<T,3> >& levelset=*melting_parameters.levelsets(1);
        ARRAY<VECTOR_2D<T> >& node_locations=levelset.grid.Node_Locations();
        for(int n=0;n<levelset.phi.m;n++) levelset.phi(n)=max(levelset.phi(n),cutout_levelset.Phi(node_locations(n)));}
}
//#####################################################################
// Function Solids_Example_Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY<VECTOR<T,3> >& V,const T time)
{
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Set_External_Velocities(V,time,id_number);
    for(int c=0;c<constrained_nodes.m;c++)V(constrained_nodes(c))=VECTOR<T,3>();
}
//#####################################################################
// Function Solids_Example_Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR<T,3> >& V,const T time)
{
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Zero_Out_Enslaved_Velocity_Nodes(V,time,id_number);
    for(int c=0;c<constrained_nodes.m;c++)V(constrained_nodes(c))=VECTOR<T,3>();
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
    if(cloth_temperature_time_constant_after && melting_parameters.reaction(1)->Max()>0){
        cloth_temperature_time_constant=cloth_temperature_time_constant_after;cloth_temperature_time_constant_after=0;}
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Update_Fluid_Parameters(dt,time);
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
    if(!use_source || time<source_start_time || time>source_stop_time) return;
    else if(source_type==BOX_SOURCE) Get_Source_Velocities(box_source,time);
    else if(source_type==CYLINDER_SOURCE) Get_Source_Velocities(cylinder_source,time);
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
template<class SOURCE> void Initialize_Phi(const SOURCE& source)
{
    for(int i=0;i<fluids_parameters.grid.m;i++) for(int j=0;j<fluids_parameters.grid.n;j++) for(int ij=0;ij<fluids_parameters.grid.mn;ij++) 
        if(source.Lazy_Inside(world_to_source*fluids_parameters.grid.X(i,j,ij)))
            fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j,ij)=-fluids_parameters.grid.dx;
        else
            fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j,ij)=fluids_parameters.grid.dx;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    ARRAY<T,VECTOR<int,3> >::copy((T)1,fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi);
    //if(!use_source) return;
    //else if(source_type==BOX_SOURCE) Initialize_Phi(box_source);
    //else if(source_type==CYLINDER_SOURCE) Initialize_Phi(cylinder_source);
}
//#####################################################################
// Function Initialize_Levelset_Velocity
//#####################################################################
void Initialize_Levelset_Velocity(const int object,ARRAY<VECTOR_2D<T> >& V)
{
    // do nothing, since we don't use levelset velocity in this example 
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
template<class SOURCE> void Adjust_Phi_With_Sources(const SOURCE& source,const T time)
{
    for(int i=0;i<fluids_parameters.grid.m;i++) for(int j=0;j<fluids_parameters.grid.n;j++) for(int ij=0;ij<fluids_parameters.grid.mn;ij++) 
//        if(source.Lazy_Inside(world_to_source*fluids_parameters.grid.X(i,j,ij)))
        if(source.Lazy_Inside(world_to_source*fluids_parameters.grid.X(i,j,ij))&&fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j,ij)>0)
            fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j,ij)=-fluids_parameters.grid.dx;
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
void Adjust_Phi_With_Sources(const T time)
{
    static bool first=true;
    if(first && restart && !restart_frame){
        ARRAY<T,VECTOR<int,3> >::copy((T)1,fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi);
        ARRAY<T,VECTOR<int,3> >::copy((T)0,fluids_parameters.density_container.array_3d);
        ARRAY<T,VECTOR<int,3> >::copy((T)283.15,fluids_parameters.temperature_container.array_3d);}
    first=false;
    FIRE_MELTING_EXAMPLE_S3D<T,RW>::Adjust_Phi_With_Sources(time);
    if(!use_source || time<source_start_time || time>source_stop_time) return;
    else if(source_type==BOX_SOURCE) Adjust_Phi_With_Sources(box_source,time);
    else if(source_type==CYLINDER_SOURCE) Adjust_Phi_With_Sources(cylinder_source,time);
}
//#####################################################################
// Function Initialize_Siggraph
//#####################################################################
void Initialize_Siggraph()
{
    cutout=true;
    ARRAY<VECTOR<T,3> ,VECTOR<int,2> > image;
    IMAGE<T>::Read(data_directory+"/Images/siggraph_print.png",image);
    T y=(T).3*image.n/image.n;
    cutout_grid.Initialize(image.m,image.n,0,1,-y/2,y/2);cutout_phi.Resize(cutout_grid);
    //for(int i=0;i<image.m;i++)for(int j=0;j<image.n;j++)LOG::cout<<image(i,j)<<std::endl;
    for(int i=0;i<image.m;i++)for(int j=0;j<image.n;j++)cutout_phi(i,j)=.5-image(i,image.n-j+1).x;
    cutout_levelset.Fast_Marching_Method();
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/cutout_levelset.phi",cutout_levelset);
    inflammable_levelset=&cutout_levelset;
}
//#####################################################################
// Function Initialize_SCA
//#####################################################################
void Initialize_SCA()
{
    cutout=true;
    ARRAY<VECTOR<T,3> ,VECTOR<int,2> > image;
    IMAGE<T>::Read(data_directory+"/Images/sca2005.png",image);
    image.Resize(1,1293+23,1,2*758-697-483);
    T x=(T)13/(1293+23);
    T y=(T).3;
    cutout_grid.Initialize(image.m,image.n,x,1-x,-y/2,y/2);cutout_phi.Resize(cutout_grid);
    //for(int i=0;i<image.m;i++)for(int j=0;j<image.n;j++)LOG::cout<<image(i,j)<<std::endl;
    for(int i=0;i<image.m;i++)for(int j=0;j<image.n;j++)cutout_phi(i,j)=.5-image(i,image.n-j+1).x;
    cutout_levelset.Fast_Marching_Method();
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/cutout_levelset.phi",cutout_levelset);
    inflammable_levelset=&cutout_levelset;
}
//#####################################################################
// Function Initialize_SCA
//#####################################################################
void Initialize_TVCG()
{
    cutout=true;
    ARRAY<VECTOR<T,3> ,VECTOR<int,2> > image;
    IMAGE<T>::Read(data_directory+"/Images/tvcg2005.png",image);
    T x=(T).02;
    T y=(T).22;
    cutout_grid.Initialize(image.m,image.n,x,1-x,-y/2,y/2);cutout_phi.Resize(cutout_grid);
    for(int i=0;i<image.m;i++)for(int j=0;j<image.n;j++)cutout_phi(i,j)=.5-image(i,image.n-j+1).x;
    cutout_levelset.Fast_Marching_Method();
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/cutout_levelset.phi",cutout_levelset);
    inflammable_levelset=&cutout_levelset;
}
//#####################################################################
// Function Refinement_Criteria
//#####################################################################
bool Refinement_Criteria(const int object,const RED_TRIANGLE<T>* triangle)
{
    if(!cutout || triangle->Depth()>=melting_parameters.maximum_depth) return false;
    for(int i=0;i<3;i++)if(fabs(cutout_levelset.Phi(triangle->Node_Location(0)))<.01) return true;
    return false;
}
//#####################################################################
};    
}
#endif

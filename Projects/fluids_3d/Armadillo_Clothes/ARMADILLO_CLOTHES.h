//#####################################################################
// Copyright 2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARMADILLO_CLOTHES
//##################################################################### 
#ifndef __ARMADILLO_CLOTHES__
#define __ARMADILLO_CLOTHES__

#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include "../FIRE_MELTING_EXAMPLE_S3D.h"
#include <Forces_And_Torques/BODY_FORCES_3D.h>
#include <Geometry/LEVELSET_IMPLICIT_SURFACE.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>
namespace PhysBAM{

template<class T,class RW=T>
class ARMADILLO_CLOTHES:public FIRE_MELTING_EXAMPLE_S3D<T,RW>
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
    using FIRE_MELTING_EXAMPLE_S3D<T,RW>::cloth_temperature_time_constant;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::cloth_ignition_temperature;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::reaction_peak_start;
    using FIRE_MELTING_EXAMPLE_S3D<T,RW>::reaction_peak_stop;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::reaction_stop;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::divergence_source_bandwidth;
    using FIRE_MELTING_EXAMPLE_S3D<T,RW>::divergence_source_rate;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::phi_source_bandwidth;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::phi_source_depth;
    using FIRE_MELTING_EXAMPLE_S3D<T,RW>::reaction_start;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::negative_reaction_coefficient_multiplier;using FIRE_MELTING_EXAMPLE_S3D<T,RW>::temperature_smoothing_steps;

    enum SOURCE_TYPE {BOX_SOURCE,CYLINDER_SOURCE};
    int example_number;
    T rho,rho_bottom,rho_top,buoyancy_constant;
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
    T armadillo_offset;

    T cloth_temperature_time_constant_after;

    int armadillo_index;
    ARRAY<T,VECTOR<int,3> > phi_object;
    bool armadillo_triangle_collisions;
 
    ARMADILLO_CLOTHES(const int example_number_input)
        :example_number(example_number_input),use_source(false),source_type(BOX_SOURCE),world_to_source(MATRIX<T,4>::Identity_Matrix()),
        initial_position(.5,-.18,.4),initial_orientation(-(T).9*pi,VECTOR<T,3>(1,0,0)),initial_velocity(),initial_angular_velocity(),side_length(1.5),m(20),n(25),armadillo_offset(0),
        cloth_temperature_time_constant_after(0),armadillo_triangle_collisions(false)
    {
        // Two way coupling
        fluids_parameters.simulate=true;
        fluids_parameters.use_solid_fluid_coupling=true;
        if(fluids_parameters.use_solid_fluid_coupling){
            use_collision_aware_signed_distance=true;
            clamp_phi_with_collision_bodies=true;
            fluids_parameters.density_container.Set_Custom_Advection(thin_shells_semi_lagrangian_density);
            fluids_parameters.temperature_container.Set_Custom_Advection(thin_shells_semi_lagrangian_temperature);
            fluids_parameters.phi_advection=&thin_shells_semi_lagrangian_phi;
            fluids_parameters.vector_advection=&thin_shells_semi_lagrangian_velocity;
            fluids_parameters.vector_advection_fuel=&thin_shells_semi_lagrangian_velocity_fuel;}
        fluids_parameters.levelset_substeps=.95;
        fluids_parameters.domain_walls[0][0]=fluids_parameters.domain_walls[0][1]=fluids_parameters.domain_walls[2][0]=fluids_parameters.domain_walls[2][1]=fluids_parameters.domain_walls[1][0]=true;
        fluids_parameters.domain_walls[1][1]=false;
        use_source=true;source_type=CYLINDER_SOURCE;cylinder_source=CYLINDER<T>(VECTOR<T,3>((T).5,(T)-.25,(T).5),VECTOR<T,3>((T).5,(T).25,(T).5),(T).15);source_velocity=VECTOR<T,3>(0,1,0);

        int r=-10;
        //example_number=2;
        std::cout<<"example number = "<<example_number<<std::endl;
        if(example_number<=0) exit(1);

        // Common parameters
        output_directory=STRING_UTILITIES::string_sprintf("Armadillo_Clothes/output%d",example_number);
        first_frame=0;last_frame=1000;frame_rate=24*2;
        restart=false;restart_frame=20;

        // Fluids parameters 
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
        solids_parameters.cfl=(T).5;

        if(example_number==13){
            //restart=false;restart_frame=39;
            r=11;
            m*=3;n*=3;
            //fluids_parameters.grid.Initialize(r*10+1,r*12+1,r*8+1,-.125,1.125,0,(T)1.5,0,1);
            fluids_parameters.grid.Initialize(r*12+1,r*11+1,r*7+1,-.25,1.25,0,(T)1.375,0,.875);
            melting_parameters.maximum_depth=1;
            solids_parameters.cfl=(T).25;
            source_start_time=.1;
            source_stop_time=1;
            cloth_ignition_temperature=297;
            cloth_temperature_time_constant=.03;
            cloth_temperature_time_constant_after=.03;
            //negative_reaction_coefficient_multiplier=.617-.590; // approximate time to heat up by 1 K
            reaction_start=.05;
            reaction_peak_start=.25;
            reaction_peak_stop=1.25;
            reaction_stop=1.5;
            divergence_source_bandwidth=2.5;
            divergence_source_rate=500;
            phi_source_bandwidth=5;
            phi_source_depth=2.5;
            PHYSBAM_FATAL_ERROR("Preroll is gone; this example is likely broken.");
            fluids_parameters.normal_flame_speed=(T).1;
            source_velocity=VECTOR<T,3>(0,1,0);
            cylinder_source=CYLINDER<T>(VECTOR<T,3>((T).5,(T)-.25,(T).34),VECTOR<T,3>((T).5,(T).063,(T).34),(T).075);
            fluids_parameters.confinement_parameter_fuel=fluids_parameters.confinement_parameter=(T).2;
            armadillo_triangle_collisions=true;
            solids_parameters.max_collision_loops=2;}
    
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
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies=true;
        solids_parameters.deformable_body_parameters.write=true;
        solids_parameters.perform_self_collision=true;
        //solids_parameters.collisions_disable_repulsions_based_on_proximity_factor=1.5;
        solids_parameters.check_initial_mesh_for_self_intersection=true;
        verbose_dt=false;

        // Melting parameters
        melting_parameters.refine_near_interface=true;
        melting_parameters.refine_for_high_deformation=true;
        //melting_parameters.use_constant_melting_speed=true;
        //melting_parameters.constant_melting_speed=0;

        // Debugging
        abort_when_dt_below=1e-7;
        write_time=true;
        write_frame_title=true;
        fluids_parameters.write_thin_shells_advection_cache=false;

        fluids_parameters.write_debug_data=true;
        //write_substeps=true;write_frame_title=true;

        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Parameters_From_Parameter_List(*this,parameter_list);
    }

    ~ARMADILLO_CLOTHES()
    {}

//#####################################################################
// Function Initialize_Armadillo
//#####################################################################
void Initialize_Armadillo()
{
    std::string tri_file=output_directory+"/armadillo.tri",phi_file=output_directory+"/armadillo.phi";
    FILE_UTILITIES::Create_Directory(output_directory);
    if(!FILE_UTILITIES::File_Exists(tri_file)){
        TETRAHEDRON_MESH tetrahedron_mesh;DEFORMABLE_PARTICLES<T,VECTOR<T,3> > particles;TETRAHEDRALIZED_VOLUME<T> tetrahedralized_volume(tetrahedron_mesh,particles);
        FILE_UTILITIES::Read_From_File<RW>(data_directory+"/Tetrahedralized_Volumes/armadillo_380K.tet",tetrahedralized_volume);
        tetrahedralized_volume.Initialize_Triangulated_Surface();tetrahedralized_volume.triangulated_surface->Discard_Valence_Zero_Particles_And_Renumber();
        FILE_UTILITIES::Write_To_File<RW>(tri_file,*tetrahedralized_volume.triangulated_surface);}
    if(!FILE_UTILITIES::File_Exists(phi_file)){
        TRIANGLE_MESH triangle_mesh;DEFORMABLE_PARTICLES<T,VECTOR<T,3> > particles;TRIANGULATED_SURFACE<T> triangulated_surface(triangle_mesh,particles);
        GRID<TV> grid;ARRAY<T,VECTOR<int,3> > phi;LEVELSET_3D<GRID<TV> > levelset(grid,phi);
        FILE_UTILITIES::Read_From_File<RW>(tri_file,triangulated_surface);
        triangulated_surface.Update_Bounding_Box();BOX_3D<T>& box=*triangulated_surface.bounding_box;
        grid=GRID_3D<T>::Create_Grid_Given_Cell_Size(box,.005*box.Size().Max(),false,5);
        phi.Resize(grid);
        LEVELSET_MAKER_UNIFORM<T> levelset_maker;
        levelset_maker.Verbose_Mode();
        levelset_maker.Set_Surface_Padding_For_Flood_Fill((T)1e-3);
        levelset_maker.Use_Fast_Marching_Method(true,20);
        levelset_maker.Compute_Level_Set(triangulated_surface,grid,phi);
        FILE_UTILITIES::Write_To_File<RW>(phi_file,levelset);}
}
//#####################################################################
// Function Initialize_Deformable_And_Rigid_Bodies
//#####################################################################
void Initialize_Deformable_And_Rigid_Bodies()
{
    {int index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/ground");
    RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list.rigid_bodies(index);
    rigid_body.Set_Coefficient_Of_Friction((T).3);rigid_body.is_static=true;}

    {Initialize_Armadillo();
    int index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(output_directory+"/armadillo",.0075,true,false,false,false);
    RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list.rigid_bodies(index);
    rigid_body.Set_Coefficient_Of_Friction(0);rigid_body.is_static=true;
    rigid_body.position=VECTOR<T,3>(.5,0,.5);
    rigid_body.Update_Bounding_Box();BOX_3D<T> box=rigid_body.Axis_Aligned_Bounding_Box();
    rigid_body.position.y-=box.ymin-5e-3;
    initial_position.y+=box.ymax-box.ymin;
    rigid_body.triangulated_surface->Update_Bounding_Box();std::cout<<*rigid_body.triangulated_surface->bounding_box<<std::endl;
    rigid_body.triangulated_surface->particles.X.array+=rigid_body.position;
    rigid_body.triangulated_surface->Update_Bounding_Box();std::cout<<*rigid_body.triangulated_surface->bounding_box<<std::endl;
    GRID<TV> armadillo_grid;ARRAY<T,VECTOR<int,3> > armadillo_phi;LEVELSET_IMPLICIT_SURFACE<T> armadillo_levelset(armadillo_grid,armadillo_phi);
    FILE_UTILITIES::Read_From_File<RW>(output_directory+"/armadillo.phi",armadillo_levelset);armadillo_levelset.Rescale(.0075);
    // initialize phi object
    armadillo_index=2;
    GRID<TV>& grid=fluids_parameters.grid;phi_object.Resize(grid);
    std::cout<<armadillo_grid<<std::endl;
    std::cout<<grid<<std::endl;
    //for(int i=-2;i<=grid.m+3;i++)for(int j=-2;j<=grid.n+3;j++)for(int ij=-2;ij<=grid.mn+3;ij++)
    for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)for(int ij=0;ij<grid.mn;ij++)
        phi_object(i,j,ij)=armadillo_levelset.Extended_Phi(grid.X(i,j,ij)-rigid_body.position);
    LEVELSET_3D<GRID<TV> > tmp_levelset(grid,phi_object);
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/tmp2.phi",tmp_levelset);
    rigid_body.position=VECTOR<T,3>();rigid_body.orientation=QUATERNION<T>();}

    {int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Embedded_Triangulated_Surface(15,1e-2);
    solids_parameters.deformable_body_parameters.list(index).use_nonembedded_self_collision=true;
    Add_Melting_Object(melting_parameters.DEFORMABLE,index);}

    if(armadillo_triangle_collisions){
        TRIANGULATED_SURFACE<T>& triangulated_surface=*solids_parameters.rigid_body_parameters.list(2)->triangulated_surface;
        triangulated_surface.particles.Store_Velocity();triangulated_surface.particles.Store_Mass();ARRAY<T>::copy(1e10,triangulated_surface.particles.mass.array);
        solids_parameters.extra_collision_surfaces.Append(&triangulated_surface);   
        solids_parameters.collision_body_list.Add_Body(solids_parameters.rigid_body_parameters.list(1));}
    else solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Initialize_Particle_Positions_And_Velocities
//#####################################################################
void Initialize_Grid(const int object,RED_GREEN_GRID_2D<T>& grid)
{
    assert(object==1);
    grid.Initialize(GRID_2D<T>(4*m+1,4*n+1,-side_length/2,side_length/2,-side_length*n/m*3/5,side_length*n/m*2/5),melting_parameters.maximum_depth);
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi(const int object,ARRAY<T>& phi)
{
    assert(object==1);
    RED_GREEN_GRID_2D<T>& grid=melting_parameters.levelsets(1)->grid;
    ARRAY<VECTOR_2D<T> >& node_locations=grid.Node_Locations();
    CIRCLE<T> head(VECTOR_2D<T>(0,.01),.14*1.05*.5);
    CIRCLE<T> arm1(VECTOR_2D<T>(.38,.07),.1*1.05*.5);
    CIRCLE<T> arm2(VECTOR_2D<T>(-.37,.1),.1*1.05*.5);
    for(int p=0;p<phi.m;p++){
        VECTOR_2D<T>& X=node_locations(p);
        phi(p)=max(-head.Signed_Distance(X),-arm1.Signed_Distance(X),-arm2.Signed_Distance(X));}
}
//#####################################################################
// Function Initialize_Particle_Positions_And_Velocities
//#####################################################################
void Initialize_Particle_Positions_And_Velocities(const int object)
{
    TRIANGULATED_SURFACE<T>& triangulated_surface=solids_parameters.deformable_body_parameters.list(object).embedded_triangulated_surface->triangulated_surface;
    DEFORMABLE_PARTICLES<T,VECTOR<T,3> >& particles=solids_parameters.deformable_body_parameters.list(object).particles;

    if(example_number>4){
        T influence_factor=3.5;T pull_factor=2;T normal_push=.2;T normal_influence_factor=2;
        CIRCLE<T> head(VECTOR_2D<T>(0,.01),.14*1.05*.5);
        CIRCLE<T> arm1(VECTOR_2D<T>(.38,.07),.1*1.05*.5);
        CIRCLE<T> arm2(VECTOR_2D<T>(-.37,.1),.1*1.05*.5);
        for(int p=0;p<particles.array_collection->Size();p++){
            VECTOR<T,3> d=particles.X(p)-VECTOR<T,3>(head.center);T r=d.Magnitude()/head.radius;
            if(r<influence_factor){
                particles.X(p)=VECTOR<T,3>(head.center)+(influence_factor+(r-influence_factor)/pull_factor)/r*d;
                if(r<normal_influence_factor) particles.X(p).z+=head.radius*normal_push*sqr(normal_influence_factor-r);}}
        for(int p=0;p<particles.array_collection->Size();p++){
            VECTOR<T,3> d=particles.X(p)-VECTOR<T,3>(arm1.center);T r=d.Magnitude()/arm1.radius;
            if(r<influence_factor){
                particles.X(p)=VECTOR<T,3>(arm1.center)+(influence_factor+(r-influence_factor)/pull_factor)/r*d;
                if(r<normal_influence_factor) particles.X(p).z+=arm1.radius*normal_push*sqr(normal_influence_factor-r);}}
        for(int p=0;p<particles.array_collection->Size();p++){
            VECTOR<T,3> d=particles.X(p)-VECTOR<T,3>(arm2.center);T r=d.Magnitude()/arm2.radius;
            if(r<influence_factor){
                particles.X(p)=VECTOR<T,3>(arm2.center)+(influence_factor+(r-influence_factor)/pull_factor)/r*d;
                if(r<normal_influence_factor) particles.X(p).z+=arm2.radius*normal_push*sqr(normal_influence_factor-r);}}}

    particles.Update_Velocity();
    triangulated_surface.Update_Bounding_Box();
    for(int p=0;p<particles.array_collection->Size();p++){
        particles.X(p)=initial_orientation.Rotate(particles.X(p));
        particles.V(p)=initial_velocity+VECTOR<T,3>::Cross_Product(initial_angular_velocity,particles.X(p));
        particles.X(p)+=initial_position;}

    VECTOR<T,3> pivot(0,initial_position.y+.2,10);QUATERNION<T> rotation((T).4*pi,VECTOR<T,3>(1,0,0));
    for(int p=0;p<particles.array_collection->Size();p++)if(particles.X(p).y>pivot.y) pivot.z=min(pivot.z,particles.X(p).z);
    for(int p=0;p<particles.array_collection->Size();p++)if(particles.X(p).y>pivot.y) particles.X(p)=pivot+rotation.Rotate(particles.X(p)-pivot);
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
    if(example_number==13 && restart && restart_frame==412){
        solids_parameters.cfl=(T).125;
        solids_parameters.implicit_solve_parameters.cg_tolerance=1e-3;}

    TRIANGULATED_SURFACE<T>& triangulated_surface=solids_parameters.deformable_body_parameters.list(1).embedded_triangulated_surface->triangulated_surface;
    solids_parameters.deformable_body_parameters.list(1).use_nonembedded_self_collision=true;

    solids_parameters.deformable_body_parameters.list(1).particles.Store_Mass();
    triangulated_surface.Set_Density(1);
    triangulated_surface.Set_Mass_Of_Particles(false);

    solids_parameters.deformable_body_parameters.list(1).Delete_Forces();
    solids_parameters.deformable_body_parameters.list(1).Add_Body_Forces(triangulated_surface,solids_parameters.gravity,solids_parameters.gravity_direction);
    solids_parameters.deformable_body_parameters.list(1).Add_Diagonalized_Linear_Finite_Volume(triangulated_surface,(T)200,(T).1,(T).05);
    solids_parameters.deformable_body_parameters.list(1).Add_Bending_Elements(triangulated_surface.triangle_mesh,(T).002,(T).005);

    std::cout<<"collision tolerance = "<<solids_parameters.deformable_body_parameters.list(1).collisions.collision_tolerance<<std::endl;

    if(fluids_parameters.use_solid_fluid_coupling){
        static int count=0;count++;
        if(count==1+restart){
            Add_To_Fluid_Simulation(solids_parameters.deformable_body_parameters.list(1),true,false,1);}}
            //Add_To_Fluid_Simulation(*solids_parameters.rigid_body_parameters.list(2),true,false,1);}}
}
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
void Update_Time_Varying_Material_Properties(const T time)
{
/*
    static T start_time=-1000;if(start_time<-100) start_time=time;
    T delta=-armadillo_offset;
    if(!restart && time<start_time+.1){
        T t=(start_time+.1-time)/.1;delta+=.1*t-.005;
        solids_parameters.deformable_body_parameters.list(1).body_forces(1)->gravity=0;}
    else{
        delta-=.005;
        solids_parameters.deformable_body_parameters.list(1).body_forces(1)->gravity=9.8;}
    std::cout<<"start_time "<<start_time<<", time "<<time<<", offset "<<armadillo_offset<<", delta "<<delta<<std::endl;
    if(fabs(delta)>1e-5){
        armadillo_offset+=delta;
        IMPLICIT_SURFACE<T>& levelset=*solids_parameters.rigid_body_parameters.list(2)->implicit_surface;
        levelset.Inflate(-delta);levelset.Update_Box();levelset.Compute_Cell_Minimum_And_Maximum(true);}
*/
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time)
{
    if(cloth_temperature_time_constant_after && melting_parameters.reaction(1).Max()>0){
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
            fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j,ij)=fluids_parameters.grid.dx;
        else
            fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j,ij)=-fluids_parameters.grid.dx;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    if(1 || !use_source) return;
    else if(source_type==BOX_SOURCE) Initialize_Phi(box_source);
    else if(source_type==CYLINDER_SOURCE) Initialize_Phi(cylinder_source);
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
        if(source.Lazy_Inside(world_to_source*fluids_parameters.grid.X(i,j,ij))&&fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j,ij)<0)
            fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j,ij)=fluids_parameters.grid.dx;
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
void Adjust_Phi_With_Sources(const T time)
{
    FIRE_MELTING_EXAMPLE_S3D<T,RW>::Adjust_Phi_With_Sources(time);
    if(!use_source || time<source_start_time || time>source_stop_time) return;
    else if(source_type==BOX_SOURCE) Adjust_Phi_With_Sources(box_source,time);
    else if(source_type==CYLINDER_SOURCE) Adjust_Phi_With_Sources(cylinder_source,time);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time)
{
    GRID<TV> &u_grid=fluids_parameters.u_grid,&v_grid=fluids_parameters.v_grid,&w_grid=fluids_parameters.w_grid;
    PROJECTION_3D<T>& projection=fluids_parameters.incompressible.projection;LAPLACE_3D<T>& elliptic_solver=*fluids_parameters.incompressible.projection.elliptic_solver;
    for(int i=0;i<u_grid.m;i++) for(int j=0;j<u_grid.n;j++) for(int ij=0;ij<u_grid.mn;ij++) if(phi_object(i,j,ij)+phi_object(i,j+1,ij)+phi_object(i,j,ij+1)+phi_object(i,j+1,ij+1)<0){
        elliptic_solver.psi_N_u(i,j,ij)=true;projection.u(i,j,ij)=0;projection.u_fuel(i,j,ij)=0;}
    for(int i=0;i<v_grid.m;i++) for(int j=0;j<v_grid.n;j++) for(int ij=0;ij<v_grid.mn;ij++) if(phi_object(i,j,ij)+phi_object(i+1,j,ij)+phi_object(i,j,ij+1)+phi_object(i+1,j,ij+1)<0){
        elliptic_solver.psi_N_v(i,j,ij)=true;projection.v(i,j,ij)=0;projection.v_fuel(i,j,ij)=0;}
    for(int i=0;i<w_grid.m;i++) for(int j=0;j<w_grid.n;j++) for(int ij=0;ij<w_grid.mn;ij++) if(phi_object(i,j,ij)+phi_object(i,j+1,ij)+phi_object(i+1,j,ij)+phi_object(i+1,j+1,ij)<0){
        elliptic_solver.psi_N_w(i,j,ij)=true;projection.w(i,j,ij)=0;projection.w_fuel(i,j,ij)=0;}
    //fluids_parameters.Extrapolate_Velocity_Into_Object(phi_object,*solids_parameters.rigid_body_parameters.list(armadillo_index),time);
    FIRE_MELTING_EXAMPLE_S3D<T,RW>::Get_Object_Velocities(dt,time);
}
//#####################################################################
};    
}
#endif

//#####################################################################
// Copyright 2003-2006, Doug Enright, Ronald Fedkiw, Nipun Kwatra, Frank Losasso, Duc Nguyen, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LIGHTHOUSE 
//##################################################################### 
#ifndef __LIGHTHOUSE__
#define __LIGHTHOUSE__

#include <Core/Math_Tools/cube.h>
#include <Core/Random_Numbers/NOISE.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Incompressible/Boundaries/BOUNDARY_OPEN_CALLBACKS.h>
#include <Dynamics/Incompressible_Flows/SPH_EVOLUTION_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input,class TV=VECTOR<T_input,3> >
class LIGHTHOUSE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>,public BOUNDARY_OPEN_CALLBACKS<TV>
{
    typedef T_input T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,TV_INT> T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::data_directory;
    using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;using BASE::stream_type;using BASE::solid_body_collection;using BASE::resolution;
    using BASE::test_number;
    using BASE::user_last_frame;
    
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    int lighthouse,cove;

    T ballistic_particles_as_percentage_of_target;
    T depth,flip_ratio;
    FLUID_COLLISION_BODY_INACCURATE_UNION<TV> inaccurate_union;

    bool use_two_way_coupling_for_sph,use_variable_density_for_sph,convert_sph_particles_to_fluid,use_analytic_divergence,use_analytic_divergence_for_expansion_only;
    int extra_cells_m_start,extra_cells_m_end,extra_cells_n_start,extra_cells_n_end;

    T epsilon;
    T scaled_depth;
    T omega;

    /***************
    example explanation:
    1. Plain ballistic particles
    2. SPH one-way coupled, variable density
    3. SPH two-way coupled, uniform density
    ***************/

    LIGHTHOUSE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>(stream_type_input,parse_args,1,fluids_parameters.WATER),
        rigid_body_collection(solid_body_collection.rigid_body_collection),inaccurate_union(*fluids_parameters.grid)
    {
        parse_args.Parse();
        first_frame=0;if(!user_last_frame) last_frame=1000;
        if(!this->user_frame_rate) frame_rate=36;
        restart=false;restart_frame=18;
        int cells=1*resolution;
        fluids_parameters.grid->Initialize(TV_INT(14*cells+1,3*cells+1,8*cells+1),RANGE<TV>(TV(-80,0,0),TV(60,30,80)));
        fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;
        fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[2][0]=true;fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.number_particles_per_cell=32;
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;
        fluids_parameters.use_removed_positive_particles=true;fluids_parameters.use_removed_negative_particles=true;
        fluids_parameters.write_removed_positive_particles=true;fluids_parameters.write_removed_negative_particles=true;
        fluids_parameters.write_debug_data=true;
        fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).save_removed_particle_times=true;
        if(!this->user_output_directory)
            output_directory=LOG::sprintf("Lighthouse/Test_%d_Lighthouse_Resolution_%d_%d_%d",test_number,(fluids_parameters.grid->counts.x-1),(fluids_parameters.grid->counts.y-1),
                (fluids_parameters.grid->counts.z-1));
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.enforce_divergence_free_extrapolation=false;
        fluids_parameters.store_particle_ids=true;
        fluids_parameters.second_order_cut_cell_method=true;
        fluids_parameters.cfl=(T)1.9;
        fluids_parameters.incompressible_iterations=40;

        // SPH parameters
        fluids_parameters.use_sph_for_removed_negative_particles=false;
        use_variable_density_for_sph=false;
        use_two_way_coupling_for_sph=false;
        convert_sph_particles_to_fluid=false;
        ballistic_particles_as_percentage_of_target=(T).25;
        flip_ratio=1;
        use_analytic_divergence=false;
        use_analytic_divergence_for_expansion_only=false;
        
        if(test_number==2){
            fluids_parameters.use_sph_for_removed_negative_particles=true;
            use_two_way_coupling_for_sph=true;
            use_variable_density_for_sph=false;
            flip_ratio=(T).5;
            use_analytic_divergence=true;
            use_analytic_divergence_for_expansion_only=true;}
        if(test_number==3){
            fluids_parameters.use_sph_for_removed_negative_particles=true;
            use_two_way_coupling_for_sph=true;}

        // MacCormack parameters
        fluids_parameters.use_maccormack_semi_lagrangian_advection=true;
        fluids_parameters.use_maccormack_compute_mask=true;
        fluids_parameters.use_maccormack_for_incompressible=true;
        fluids_parameters.bandwidth_without_maccormack_near_interface=1;

        // Wave parameters
        depth=(T)5.5;
        epsilon=(T).5;
        scaled_depth=depth/30;
        omega=(T)1.5;
    }

    ~LIGHTHOUSE() 
    {}

//#####################################################################
// Function Get_Wave_Height
//#####################################################################
T Get_Wave_Height(const TV& X,const TV_INT& cell_index,const T time=0)
{
    T x=-(X.x+10*time-50+X.z/10)/50;
    return 50*((T)1./(2*(T)pi)*(epsilon*cos(2*(T)pi*x)+(T).5*sqr(epsilon)*cos(4*(T)pi*x)+(T)3/8*cube(epsilon)*cos(6*(T)pi*x)));
}
//#####################################################################
// Function Get_Wave_Velocity
//#####################################################################
T Get_Wave_Velocity(const TV& X,const int axis,const T time=0)
{
    T x=-(X.x+10*time-50+X.z/10)/50,y=X.y/50;
    T scaled_omega=omega*(-X.z/140+1);
    T first_part=scaled_omega*exp(-2*(T)pi*(scaled_depth-y));
    if(axis==1) 
        return -((2*(T)pi)*first_part*(epsilon*cos(2*(T)pi*x)+(T).5*sqr(epsilon)*cos(4*(T)pi*x)+(T)3/8*cube(epsilon)*cos(6*(T)pi*x)))-3;
    else if(axis==2)
        return (2*(T)pi)*first_part*(epsilon*sin(2*(T)pi*x)+(T).5*sqr(epsilon)*sin(4*(T)pi*x)+(T)3/8*cube(epsilon)*sin(6*(T)pi*x));
    return 0;
}
//#####################################################################
// Function Get_Wave_Attenuation
//#####################################################################
T Get_Wave_Attenuation(const TV X)
{
    return clamp((X.x+30)/30,(T)0,(T)1);
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() override
{
    GRID<TV>& grid=*fluids_parameters.grid;
    ARRAY<T,VECTOR<int,3> >& phi=fluids_parameters.particle_levelset_evolution->phi;

    LOG::cout<<"setting phi values"<<std::endl;
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        TV X=iterator.Location();
        phi(iterator.Cell_Index())=X.y-depth-Get_Wave_Attenuation(X)*Get_Wave_Height(X,iterator.Cell_Index());}
    
    for(int i=0;i<grid.counts.x;i++) for(int j=0;j<grid.counts.y;j++) for(int ij=0;ij<grid.counts.z;ij++){
        TV X=grid.X(TV_INT(i,j,ij));
        fluids_parameters.particle_levelset_evolution->phi(i,j,ij)=max(fluids_parameters.particle_levelset_evolution->phi(i,j,ij),-rigid_body_collection.Rigid_Body(lighthouse).Implicit_Geometry_Extended_Value(X),-rigid_body_collection.Rigid_Body(cove).Implicit_Geometry_Extended_Value(X));}
}
//#####################################################################
// Function Get_Dirichlet_Boundary_Conditions
//#####################################################################
void Set_Dirichlet_Boundary_Conditions(const T time) override
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>::Set_Dirichlet_Boundary_Conditions(time);
    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters_uniform=dynamic_cast<FLUIDS_PARAMETERS_UNIFORM<TV>&>(fluids_parameters);
    if(fluids_parameters_uniform.mpi_grid && fluids_parameters_uniform.mpi_grid->Neighbor(1,2)) return;

    ARRAY<bool,FACE_INDEX<TV::m> >& psi_N=fluids_parameters.incompressible->projection.elliptic_solver->psi_N;
    ARRAY<T,FACE_INDEX<TV::m> >& face_velocities=fluid_collection.incompressible_fluid_collection.face_velocities;

    RANGE<VECTOR<int,3> > right_grid_cells=RANGE<VECTOR<int,3> >(TV_INT(fluids_parameters.grid->counts.x-2,1,1),fluids_parameters.grid->Numbers_Of_Cells());
    for(int axis=0;axis<3;axis++){
        RANGE<VECTOR<int,3> > right_grid_faces=right_grid_cells+RANGE<VECTOR<int,3> >(TV_INT(),TV_INT::Axis_Vector(axis));
        for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid,right_grid_faces,axis);iterator.Valid();iterator.Next()){TV_INT face=iterator.Face_Index();
            psi_N.Component(axis)(face)=true;face_velocities.Component(axis)(face)=Get_Wave_Velocity(iterator.Location(),axis,time);}}
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
bool Adjust_Phi_With_Sources(const T time) override
{
    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters_uniform=dynamic_cast<FLUIDS_PARAMETERS_UNIFORM<TV>&>(fluids_parameters);
    if(fluids_parameters_uniform.mpi_grid && fluids_parameters_uniform.mpi_grid->Neighbor(1,2)) return false;
    
    RANGE<VECTOR<int,3> > right_grid_cells=RANGE<VECTOR<int,3> >(TV_INT(fluids_parameters.grid->counts.x-2,1,1),fluids_parameters.grid->Numbers_Of_Cells());
    for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid,right_grid_cells);iterator.Valid();iterator.Next()){
        TV X=iterator.Location();
        fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index())=X.y-depth-Get_Wave_Height(X,iterator.Cell_Index(),time);}

    T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES& removed_negative_particles=fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).removed_negative_particles;
    RANGE<VECTOR<int,3> > right_grid_nodes=RANGE<VECTOR<int,3> >(TV_INT(fluids_parameters.grid->counts.x-9,-2,-2),TV_INT(fluids_parameters.grid->counts.x+3,fluids_parameters.grid->counts.y+3,fluids_parameters.grid->counts.z+3));
    for(NODE_ITERATOR<TV> iterator(*fluids_parameters.grid,right_grid_nodes);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();
        if(removed_negative_particles(block)){removed_negative_particles(block)->Delete_All_Elements();removed_negative_particles(block)=0;}}
    return false;
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) override
{
    Write_To_File(stream_type,LOG::sprintf("%s/removed_particle_times.%d",output_directory.c_str(),frame),
        fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).removed_particle_times);
    std::cout<<"Wrote "<<fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).removed_particle_times.m<<" number of particles"<<std::endl;
    fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).removed_particle_times.Clean_Memory();
}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() override
{
    fluids_parameters.Use_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() override
{
    ARRAY<T,FACE_INDEX<TV::m> >& face_velocities=fluid_collection.incompressible_fluid_collection.face_velocities;
    ARRAY<T,VECTOR<int,3> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
        TV X=iterator.Location();
        VECTOR<double,3> noise_v;VECTOR<double,3> Xd((T).09*X);NOISE<double>::Noise3(Xd,noise_v,3);noise_v-=VECTOR<double,3>((T).25,(T).25,(T).25);
        if((T).5*(phi(iterator.First_Cell_Index())+phi(iterator.Second_Cell_Index()))<1)
            face_velocities(iterator.Axis(),iterator.Face_Index())=(T)2.5*(T)noise_v[iterator.Axis()]+Get_Wave_Attenuation(X)*Get_Wave_Velocity(X,iterator.Axis());}
}
//#####################################################################
// Function Initial_Phi_Object
//#####################################################################
T Initial_Phi_Object(const VECTOR<T,3>& X) const
{
    return min(rigid_body_collection.Rigid_Body(lighthouse).Implicit_Geometry_Extended_Value(X),rigid_body_collection.Rigid_Body(cove).Implicit_Geometry_Extended_Value(X));
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{    
    lighthouse=rigid_body_collection.Add_Rigid_Body(data_directory+"/Rigid_Bodies/lighthouse",(T)2,true,true,false);
    rigid_body_collection.rigid_body_particles.frame(lighthouse).t=VECTOR<T,3>((T)25,(T)-2,(T)30);rigid_body_collection.Rigid_Body(lighthouse).Is_Kinematic()=true;
    cove=rigid_body_collection.Add_Rigid_Body(data_directory+"/Rigid_Bodies/cove",(T)1,true,true,false);
    rigid_body_collection.rigid_body_particles.frame(cove).t=VECTOR<T,3>((T)-80,(T)11,(T)40);rigid_body_collection.Rigid_Body(cove).Is_Kinematic()=true;
    inaccurate_union.collision_bodies.Add_Bodies(rigid_body_collection);
    fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(&inaccurate_union,0,false);
}
//#####################################################################
// Function Initialize_SPH_Particles
//#####################################################################
void Initialize_SPH_Particles() override
{
    GRID<TV>& grid=fluids_parameters.particle_levelset_evolution->grid;
    SPH_EVOLUTION_UNIFORM<TV>& sph_evolution=*fluids_parameters.sph_evolution;
    sph_evolution.use_two_way_coupling=use_two_way_coupling_for_sph;
    sph_evolution.use_variable_density_solve=use_variable_density_for_sph;
    sph_evolution.convert_particles_to_fluid=convert_sph_particles_to_fluid;
    sph_evolution.target_particles_per_unit_volume=(T)fluids_parameters.number_particles_per_cell/grid.Cell_Size();
    sph_evolution.ballistic_particles_per_unit_volume=ballistic_particles_as_percentage_of_target*sph_evolution.target_particles_per_unit_volume;
    sph_evolution.flip_ratio=flip_ratio;
    sph_evolution.use_analytic_divergence=use_analytic_divergence;
    sph_evolution.use_analytic_divergence_for_expansion_only=use_analytic_divergence_for_expansion_only;
}
//#####################################################################
// Function Get_Body_Force
//#####################################################################
void Get_Body_Force(ARRAY<T,FACE_INDEX<TV::m> >& force,const T dt,const T time) override
{
    BASE::Get_Body_Force(force,dt,time);
}
//#####################################################################
};
}
#endif

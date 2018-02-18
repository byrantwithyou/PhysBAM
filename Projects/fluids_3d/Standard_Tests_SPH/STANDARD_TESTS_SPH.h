//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS_SPH
//#####################################################################
#ifndef __STANDARD_TESTS_SPH__
#define __STANDARD_TESTS_SPH__

#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Dynamics/Incompressible_Flows/SPH_EVOLUTION_UNIFORM.h>
#include <Dynamics/Particles/SPH_PARTICLES.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Dynamics/Standard_Tests/WATER_STANDARD_TESTS_3D.h>
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS_SPH:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,3> >
{
    typedef T_input T;
public:
    typedef VECTOR<T,3> TV;typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::data_directory;using BASE::solid_body_collection;using BASE::test_number;
    using BASE::Adjust_Phi_With_Sources;using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;using BASE::Get_Object_Velocities; // silence -Woverloaded-virtual
    using BASE::resolution;using BASE::user_last_frame;

    WATER_STANDARD_TESTS_3D<TV> tests;
    int number_of_particles;
    T wall_damping;

    STANDARD_TESTS_SPH(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>(stream_type_input,parse_args,0,fluids_parameters.SPH),
        tests(*this,fluids_parameters,solid_body_collection.rigid_body_collection),wall_damping((T).5)
    {
        parse_args.Parse();
        tests.Initialize(test_number,resolution);
        *fluids_parameters.grid=tests.grid;
        fluids_parameters.number_particles_per_cell=16;
        fluids_parameters.write_ghost_values=true;
        fluids_parameters.store_particle_ids=true;
        if(!user_last_frame) last_frame=1000;

        if(tests.test_number==1) number_of_particles=100000;
        else if(tests.test_number==2) number_of_particles=15000;
        else{LOG::cerr<<"unrecognzed SPH test number"<<std::endl;exit(1);}

        if(!this->user_output_directory)
            output_directory=LOG::sprintf("Standard_Tests_SPH/Test_%d__Resolution_%d_%d",test_number,(tests.grid.counts.x-1),(tests.grid.counts.y-1));
        LOG::cout<<"Running SPH simulation to "<<output_directory<<std::endl;
    }

    virtual ~STANDARD_TESTS_SPH()
    {}

//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_SPH_Particles() override
{
    fluids_parameters.sph_evolution->flip_ratio=0;
    GRID<TV>& grid=*fluids_parameters.grid;
    RANDOM_NUMBERS<T> random;
    SPH_PARTICLES<TV>& sph_particles=fluids_parameters.sph_evolution->sph_particles;

    fluids_parameters.sph_evolution->use_variable_density_solve=true;

    //sph_particles.Preallocate(number_of_particles);
    if(tests.test_number==1){
        fluids_parameters.sph_evolution->target_particles_per_unit_volume=10000;
        for(int i=0;i<number_of_particles;i++){
            TV X=random.Get_Uniform_Vector(TV((T).7,0,(T).4),TV(1,1,(T).6));
            int id=sph_particles.Add_Element();
            sph_particles.X(id)=X;}}
    else if(tests.test_number==2){
        int left_particles_number=number_of_particles/3,right_particles_number=number_of_particles-left_particles_number;
        fluids_parameters.sph_evolution->target_particles_per_unit_volume=number_of_particles*(T).5;
        for(int i=0;i<left_particles_number;i++){
            TV X=random.Get_Uniform_Vector(grid.domain.min_corner,TV((T).5*(grid.domain.max_corner.x-grid.domain.min_corner.x),grid.domain.max_corner.y,grid.domain.max_corner.z));
            int id=sph_particles.Add_Element();
            sph_particles.X(id)=X;}
        for(int i=0;i<right_particles_number;i++){
            TV X=random.Get_Uniform_Vector(TV((T).5*(grid.domain.max_corner.x-grid.domain.min_corner.x),grid.domain.min_corner.y,grid.domain.min_corner.z),grid.domain.max_corner);
            int id=sph_particles.Add_Element();
            sph_particles.X(id)=X;}}

    this->limit_dt=[this](T& dt,T time){tests.Limit_Dt(dt,time);};
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time) override
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>::Update_Fluid_Parameters(dt,time);
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Object_Velocities(const bool set_densities_for_coupling,const T dt,const T time) // TODO: override Problem
{
    if(tests.test_number==1 && time<4)
        for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) if(iterator.Location().x<(T).7){
            fluids_parameters.incompressible->projection.elliptic_solver->psi_N(iterator.Axis(),iterator.Face_Index())=true;
            fluids_parameters.incompressible->projection.face_velocities(iterator.Axis(),iterator.Face_Index())=0;}
}
//#####################################################################
// Function Adjust_SPH_Particle_For_Domain_Boundaries
//#####################################################################
void Adjust_SPH_Particle_For_Domain_Boundaries(SPH_PARTICLES<TV>& particles,const int index,TV& V,const T dt,const T time)const override
{
    GRID<TV> grid=*fluids_parameters.grid;
    if(tests.test_number==1 && time<4) grid.domain.min_corner.x=(T).7;
    TV new_X=particles.X(index)+dt*V;
    if(fluids_parameters.domain_walls[0][0]&&new_X.x<grid.domain.min_corner.x){V.x=0;new_X.x=grid.domain.min_corner.x;particles.X(index)=new_X-dt*V;}
    if(fluids_parameters.domain_walls[0][1]&&new_X.x>grid.domain.max_corner.x){V.x=0;new_X.x=grid.domain.max_corner.x;particles.X(index)=new_X-dt*V;}
    if(fluids_parameters.domain_walls[1][0]&&new_X.y<grid.domain.min_corner.y){V.y=00;new_X.y=grid.domain.min_corner.y;particles.X(index)=new_X-dt*V;}
    if(fluids_parameters.domain_walls[1][1]&&new_X.y>grid.domain.max_corner.y){V.y=0;new_X.y=grid.domain.max_corner.y;particles.X(index)=new_X-dt*V;}
    if(fluids_parameters.domain_walls[2][0]&&new_X.z<grid.domain.min_corner.z){V.z=0;new_X.z=grid.domain.min_corner.z;particles.X(index)=new_X-dt*V;}
    if(fluids_parameters.domain_walls[2][1]&&new_X.z>grid.domain.max_corner.z){V.z=0;new_X.z=grid.domain.max_corner.z;particles.X(index)=new_X-dt*V;}
}
//#####################################################################
};
}
#endif

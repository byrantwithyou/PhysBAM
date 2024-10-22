//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS_SPH
//#####################################################################
#ifndef __STANDARD_TESTS_SPH__
#define __STANDARD_TESTS_SPH__
#include <Core/Matrices/FRAME.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Dynamics/Incompressible_Flows/SPH_EVOLUTION_UNIFORM.h>
#include <Dynamics/Particles/SPH_PARTICLES.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
// 1 breaking dam
// 2 variable density seed showing divergence control
// 3 falling drop
// 5 filling box with sphere obstacle
// 6 falling box into pool
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS_SPH:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,2> >
{
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::viewer_dir;using BASE::test_number;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::data_directory;using BASE::fluid_collection;using BASE::solid_body_collection;
    using BASE::Adjust_Phi_With_Sources;using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;using BASE::Get_Object_Velocities; // silence -Woverloaded-virtual
    using BASE::stream_type;using BASE::resolution;
    using BASE::user_last_frame;
    
    int number_of_particles;
    T wall_damping;
    GRID<TV> grid;
    RANDOM_NUMBERS<T> random;
    RIGID_BODY<TV>* rigid;
    RIGID_BODY<TV>* rigid2; // not used right now
    SPHERE<TV> sphere;

    STANDARD_TESTS_SPH(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>(stream_type_input,parse_args,0,fluids_parameters.SPH),
        wall_damping((T).5),sphere(TV((T).5,(T).4),(T).2)
    {
        parse_args.Parse();
        fluids_parameters.number_particles_per_cell=16;
        fluids_parameters.write_ghost_values=true;
        fluids_parameters.store_particle_ids=true;
        if(!user_last_frame) last_frame=1000;

        int cells=1*resolution;
        if(test_number==1){
            if(!user_last_frame) last_frame=200;
            grid.Initialize(TV_INT(15*cells+1,10*cells+1),RANGE<TV>(TV(0,0),TV((T)1.5,1)));}
        else if(test_number==2){
            if(!user_last_frame) last_frame=100;
            grid.Initialize(TV_INT(10*cells+1,15*cells+1),RANGE<TV>(TV(0,0),TV(1,(T)1.5)));}
        else if(test_number==3){
            if(!user_last_frame) last_frame=150;
            grid.Initialize(TV_INT(10*cells+1,15*cells+1),RANGE<TV>(TV(0,0),TV(1,(T)1.5)));}
        else if(test_number==5){
            grid.Initialize(TV_INT(10*cells+1,10*cells+1),RANGE<TV>(TV(0,0),TV(1,1)));}
        else if(test_number==6){   
            grid.Initialize(TV_INT(10*cells+1,10*cells+1),RANGE<TV>(TV(0,0),TV(1,1)),true);
            int rigid_id=solid_body_collection.rigid_body_collection.Add_Rigid_Body(data_directory+"/Rigid_Bodies_2D/square_refined",(T).1,true,true,false);
            rigid=&solid_body_collection.rigid_body_collection.Rigid_Body(rigid_id);
            //rigid->is_kinematic=true;
            rigid->Frame().t=TV((T).5,(T).4);
            rigid->coefficient_of_restitution=(T)0;
            rigid2=&solid_body_collection.rigid_body_collection.Rigid_Body(rigid_id);
            //rigid->is_kinematic=true;
            //rigid2->Frame().t=TV((T).1,(T).4);
            //rigid2->is_kinematic=true;
            int ground_id=solid_body_collection.rigid_body_collection.Add_Rigid_Body(data_directory+"/Rigid_Bodies_2D/ground",(T)1,true,true,false);
            solid_body_collection.rigid_body_collection.Rigid_Body(ground_id).is_static=true;
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        }
        else{LOG::cerr<<"unrecognzed SPH test number"<<std::endl;exit(1);}
        *fluids_parameters.grid=grid;

        if(!this->user_output_directory)
            viewer_dir.output_directory=LOG::sprintf("Standard_Tests_SPH/Test_%d__Resolution_%d_%d",test_number,(grid.counts.x-1),(grid.counts.y-1));
        LOG::cout<<"Running SPH simulation to "<<viewer_dir.output_directory<<std::endl;
    }

    virtual ~STANDARD_TESTS_SPH()
    {}
    void Initialize_Bodies() override
    {
        solid_body_collection.Add_Force(new GRAVITY<TV>(solid_body_collection.deformable_body_collection.particles,
                solid_body_collection.rigid_body_collection,true,true));
    }

//#####################################################################
// Function Initialize_SPH_Particles
//#####################################################################
void Initialize_SPH_Particles() override
{
    fluids_parameters.sph_evolution->flip_ratio=(T).9;
    fluids_parameters.sph_evolution->divergence_expansion_multiplier=(T)1;
    GRID<TV>& grid=*fluids_parameters.grid;
    RANDOM_NUMBERS<T> random;
    SPH_PARTICLES<TV>& sph_particles=fluids_parameters.sph_evolution->sph_particles;
    fluids_parameters.sph_evolution->use_variable_density_solve=false;

    fluids_parameters.sph_evolution->target_particles_per_unit_volume=20000;
    if(test_number==1){
        int particles_added=0;
        RANGE<TV> seed_box(TV((T)0.7,0),TV(1,1));
        number_of_particles=(int)(seed_box.Size()*(T)fluids_parameters.sph_evolution->target_particles_per_unit_volume);
        while(particles_added<number_of_particles){
            TV X;
            random.Fill_Uniform(X,grid.domain);
            T phi=seed_box.Signed_Distance(X);
            if(phi<=0){
                particles_added++;
                int id=sph_particles.Add_Element();
                sph_particles.X(id)=X;}}}
    else if(test_number==2){
        number_of_particles=10000;
        int left_particles_number=number_of_particles/3,right_particles_number=number_of_particles-left_particles_number;
        fluids_parameters.sph_evolution->target_particles_per_unit_volume=number_of_particles*(T).5;
        for(int i=0;i<left_particles_number;i++){
            TV X;
            random.Fill_Uniform(X,grid.domain.min_corner,TV((T).5*(grid.domain.max_corner.x-grid.domain.min_corner.x),grid.domain.max_corner.y));
            int id=sph_particles.Add_Element();
            sph_particles.X(id)=X;}
        for(int i=0;i<right_particles_number;i++){
            TV X;
            random.Fill_Uniform(X,TV((T).5*(grid.domain.max_corner.x-grid.domain.min_corner.x),grid.domain.min_corner.y),grid.domain.max_corner);
            int id=sph_particles.Add_Element();
            sph_particles.X(id)=X;}}
    else if(test_number==3){
        int particles_added=0;
        SPHERE<TV> seed_sphere(TV((T).5,(T).75),(T).2);
        RANGE<TV> seed_box(TV(),TV(1,(T).5));
        // assume they don't overlap to compute total area
        number_of_particles=(int)((seed_box.Size()+seed_sphere.Size())*(T)fluids_parameters.sph_evolution->target_particles_per_unit_volume);
        while(particles_added<number_of_particles){
            TV X;
            random.Fill_Uniform(X,grid.domain);
            T phi=min(seed_sphere.Signed_Distance(X),seed_box.Signed_Distance(X));
            if(phi<=0){
                particles_added++;
                int id=sph_particles.Add_Element();
                sph_particles.X(id)=X;}}}
    else if(test_number==5);
    else if(test_number==6){
        RANGE<TV> initial_seed_box(TV(),TV(1,(T).25));
        number_of_particles=(int)(initial_seed_box.Size()*(T)fluids_parameters.sph_evolution->target_particles_per_unit_volume);
        int particles_added=0;
        while(particles_added<number_of_particles){
            TV X;
            random.Fill_Uniform(X,grid.domain);
            T phi=initial_seed_box.Signed_Distance(X);
            if(phi<=0){
                particles_added++;
                int id=sph_particles.Add_Element();
                sph_particles.X(id)=X;}}}
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
//######################### ############################################
void Get_Object_Velocities(LAPLACE_UNIFORM<TV>* elliptic_solver,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time) override
{
    if(test_number==1 && time<4)
        for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) if(iterator.Location().x<(T).7){
            fluids_parameters.incompressible->projection.elliptic_solver->psi_N(iterator.Axis(),iterator.Face_Index())=true;
            fluid_collection.incompressible_fluid_collection.face_velocities(iterator.Axis(),iterator.Face_Index())=0;}
    if(test_number==5){
        for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) if(sphere.Lazy_Inside(iterator.Location())){
            elliptic_solver->psi_N(iterator.Axis(),iterator.Face_Index())=true;
            fluid_collection.incompressible_fluid_collection.face_velocities(iterator.Axis(),iterator.Face_Index())=0;}}
    if(test_number==6){
        for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) if(rigid->Implicit_Geometry_Lazy_Inside(iterator.Location())){
            elliptic_solver->psi_N(iterator.Axis(),iterator.Face_Index())=true;
            fluid_collection.incompressible_fluid_collection.face_velocities(iterator.Axis(),iterator.Face_Index())=rigid->Pointwise_Object_Velocity(iterator.Location())[iterator.Axis()];}}   
}
//#####################################################################
// Function Adjust_SPH_Particle_For_Domain_Boundaries
//#####################################################################
void Adjust_SPH_Particle_For_Domain_Boundaries(SPH_PARTICLES<TV>& particles,const int index,TV& V,const T dt,const T time)const override
{
    GRID<TV> grid=*fluids_parameters.grid;
    if(test_number==1 && time<4) grid.domain.min_corner.x=(T).7;
    TV new_X=particles.X(index)+dt*V;
    if(test_number==5 && sphere.Lazy_Inside(new_X)){
        new_X=sphere.Surface(new_X);
        V=(new_X-particles.X(index))/dt;}
    T phi;
    if(test_number==6 && rigid->Implicit_Geometry_Lazy_Inside_And_Value(new_X,phi)){
        TV normal=rigid->Implicit_Geometry_Normal(new_X);
        new_X-=normal*phi;
        V=(new_X-particles.X(index))/dt;}
    if(fluids_parameters.domain_walls[0][0]&&new_X.x<grid.domain.min_corner.x){V.x=0;new_X.x=grid.domain.min_corner.x;particles.X(index)=new_X-dt*V;}
    if(fluids_parameters.domain_walls[0][1]&&new_X.x>grid.domain.max_corner.x){V.x=0;new_X.x=grid.domain.max_corner.x;particles.X(index)=new_X-dt*V;}
    if(fluids_parameters.domain_walls[1][0]&&new_X.y<grid.domain.min_corner.y){V.y=00;new_X.y=grid.domain.min_corner.y;particles.X(index)=new_X-dt*V;}
    if(fluids_parameters.domain_walls[1][1]&&new_X.y>grid.domain.max_corner.y){V.y=0;new_X.y=grid.domain.max_corner.y;particles.X(index)=new_X-dt*V;}
        
}
//#####################################################################
// Function Add_SPH_Particles_For_Sources
//#####################################################################
void Add_SPH_Particles_For_Sources(const T dt,const T time) override
{
    SPH_PARTICLES<TV>& sph_particles=fluids_parameters.sph_evolution->sph_particles;
    if(test_number==5){
        RANGE<TV> source_box(TV((T).2,(T).8),TV((T).3,(T).9));
        int particles_per_second=2000;
        int particles_to_add=max(1,(int)((T)particles_per_second*dt));
        for(int i=0;i<particles_to_add;i++)
            random.Fill_Uniform(sph_particles.X(sph_particles.Add_Element()),source_box);}
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) override
{
    return false;
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) override
{
    //if(rigid2->particle_index==id) frame.t=time<2?TV((T).1,.1):TV((T).1,.1+time);
}
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2005-2007, Geoffrey Irving, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLOW_PAST_EFTYCHIS
//#####################################################################
#ifndef __FLOW_PAST_EFTYCHIS__
#define __FLOW_PAST_EFTYCHIS__

#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>

namespace PhysBAM{

template<class T_input>
class FLOW_PAST_EFTYCHIS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,3> >
{
    typedef T_input T;
public:
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::viewer_dir;using BASE::Adjust_Phi_With_Sources;
    using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;using BASE::data_directory;
    using BASE::solid_body_collection;using BASE::stream_type;using BASE::test_number;
    using BASE::user_last_frame;
    
    T flow_speed;
    int body;
    SPHERE<TV> source_sphere;
    TV source_vector;


    FLOW_PAST_EFTYCHIS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>(stream_type_input,parse_args,0,fluids_parameters.SMOKE)
    {
        parse_args.Parse();

        GRID<TV> full_grid(TV_INT(75,75,75),RANGE<TV>(TV((T)-.3,(T)0.35,(T)-.3),TV((T).3,(T).95,(T).3)));
        *fluids_parameters.grid=full_grid;
        fluids_parameters.use_vorticity_confinement=true;fluids_parameters.confinement_parameter=(T).5;
        fluids_parameters.write_debug_data=true;

        if(!this->user_output_directory)
            viewer_dir.output_directory=LOG::sprintf("Flow_Past_Eftychis/Test_%d",test_number);

        if(test_number==1){
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;
            fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[2][0]=true;fluids_parameters.domain_walls[2][1]=true;
            fluids_parameters.confinement_parameter=(T).5;
        }
        else if(test_number==2){
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=false;fluids_parameters.domain_walls[1][0]=true;
            fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[2][0]=false;fluids_parameters.domain_walls[2][1]=true;
            fluids_parameters.gravity=TV();
            source_sphere.radius=(T).1;
            source_sphere.center=TV((T)-.24,(T).62,(T).25);
            source_vector=TV((T).66666,0,(T)-.66666);
            fluids_parameters.confinement_parameter=(T).15;
        }
        else if(test_number==3){
            fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=false;fluids_parameters.domain_walls[1][0]=true;
            fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[2][0]=false;fluids_parameters.domain_walls[2][1]=true;
            //fluids_parameters.gravity=TV();
            source_sphere.radius=(T).1;
            source_sphere.center=TV((T)-.24,(T).62,(T).25);
            source_vector=TV((T).66666,0,(T)-.66666);
            fluids_parameters.gravity=source_vector*fluids_parameters.gravity.Magnitude();
            fluids_parameters.confinement_parameter=(T).15;
        }
 
        if(!user_last_frame) last_frame=600;
        if(!this->user_frame_rate) frame_rate=120;
        
        flow_speed=(T)1;
    }

    ~FLOW_PAST_EFTYCHIS()
    {}

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
    if(test_number==1){
        for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) 
            fluid_collection.incompressible_fluid_collection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=(iterator.Axis()==3?-flow_speed:(T)0);}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    body=rigid_body_collection.Add_Rigid_Body("/solver/vol3/hair1/data/body_recede/Rigid_Bodies/body_recede",(T)1,true,true,false);
    rigid_body_collection.Rigid_Body(body).Is_Kinematic()=true;
    fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection);
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time) override
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>::Update_Fluid_Parameters(dt,time);
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
bool Adjust_Phi_With_Sources(const T time) override
{
    return false;
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,3> >*& cell_centered_mask,const T time) override
{
}
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Adjust_Density_And_Temperature_With_Sources(const T time) override
{
    if(test_number==1){
        if(!fluids_parameters.mpi_grid || !fluids_parameters.mpi_grid->Neighbor(3,2))
            for(int i=0;i<fluids_parameters.grid->counts.x;i++) for(int j=0;j<fluids_parameters.grid->counts.y;j++){
                TV X=fluids_parameters.grid->X(VECTOR<int,3>(i,j,0));
                if(VECTOR<T,2>(X.x,X.y-(T).6).Magnitude()<(T).1)
                    fluids_parameters.density_container.density(i,j,fluids_parameters.grid->counts.z-1)=(T)1;}}
    else if(test_number==2 || test_number==3){
        for(int i=0;i<fluids_parameters.grid->counts.x;i++) for(int j=0;j<fluids_parameters.grid->counts.y;j++) for(int ij=0;ij<fluids_parameters.grid->counts.z;ij++){
            if(source_sphere.Lazy_Inside(fluids_parameters.grid->X(TV_INT(i,j,ij)))){
                fluids_parameters.density_container.density(i,j,ij)=(T)1;}}}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,const T time) override
{
    GRID<TV> u_grid=fluids_parameters.grid->Get_Face_Grid(0),v_grid=fluids_parameters.grid->Get_Face_Grid(1),w_grid=fluids_parameters.grid->Get_Face_Grid(2);
    LAPLACE_UNIFORM<TV>& elliptic_solver=*fluids_parameters.incompressible->projection.elliptic_solver;

    if(test_number==1){
        if(!fluids_parameters.mpi_grid || !fluids_parameters.mpi_grid->Neighbor(3,1))
            for(int i=0;i<w_grid.counts.x;i++) for(int j=0;j<w_grid.counts.y;j++){
                fluid_collection.incompressible_fluid_collection.face_velocities.Component(2)(i,j,2)=-flow_speed;
                elliptic_solver.psi_N.Component(2)(i,j,2)=true;}
        if(!fluids_parameters.mpi_grid || !fluids_parameters.mpi_grid->Neighbor(3,2))
            for(int i=0;i<w_grid.counts.x;i++) for(int j=0;j<w_grid.counts.y;j++){
                fluid_collection.incompressible_fluid_collection.face_velocities.Component(2)(i,j,w_grid.counts.z-1)=-flow_speed;
                elliptic_solver.psi_N.Component(2)(i,j,w_grid.counts.z-1)=true;}
    }
    else if(test_number==2 || test_number==3){
        for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) if(source_sphere.Lazy_Inside(iterator.Location())) {
            fluid_collection.incompressible_fluid_collection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=source_vector[iterator.Axis()];
            elliptic_solver.psi_N.Component(iterator.Axis())(iterator.Face_Index())=true;}}
}
//#####################################################################
// Function Get_Variable_Viscosity
//#####################################################################
void Get_Variable_Viscosity(ARRAY<T,VECTOR<int,3> >& viscosity,const T time) override
{
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) override
{
    frame=FRAME<TV>();
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) override
{
    twist=TWIST<TV>();
    return false;
}
//#####################################################################
// Function Initialize_SPH_Particles
//#####################################################################
void Initialize_SPH_Particles() override
{
}
//#####################################################################
// Function Get_Analytic_Velocities
//#####################################################################
void Get_Analytic_Velocities(const T time) const override
{
}
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2007 Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOD_ST_DROP
//#####################################################################
#ifndef __SOD_ST_DROP__
#define __SOD_ST_DROP__

#include <fstream>
#include <iostream>

#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Vectors/VECTOR_UTILITIES.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Fluids/PhysBAM_Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP.h>
#include <Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_RF.h>
#include <Fluids/PhysBAM_Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class SOD_ST_DROP:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,1> > >
{
public:
    typedef T_input T;typedef VECTOR<T,1> TV;typedef VECTOR<int,1> TV_INT;typedef GRID<TV> T_GRID;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    typedef VECTOR<T,2*T_GRID::dimension> T_FACE_VECTOR;typedef VECTOR<TV,2*T_GRID::dimension> TV_FACE_VECTOR;
    using BASE::initial_time;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;
    using BASE::parse_args;using BASE::test_number;using BASE::resolution;

    /* 1: Shock impinging on rho=1000kg/m^3 drop
       2: Shock impinging on rho=10kg/m^3 drop
       2: rho=1000kg/m^3 drop moving to right
       3: rho=10kg/m^3 drop moving to the right */

    int eno_scheme;
    int eno_order;
    int rk_order;
    T cfl_number;
    bool timesplit;
    bool implicit_rk;
    bool exact;

    SOD_ST_DROP(const STREAM_TYPE stream_type)
        :BASE(stream_type,1,fluids_parameters.COMPRESSIBLE),eno_scheme(1),eno_order(2),rk_order(3),cfl_number((T).5),timesplit(false),
        implicit_rk(false),exact(false)
    {
    }
    
    virtual ~SOD_ST_DROP() {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add("-eno_scheme",&eno_scheme,"eno_scheme","eno scheme");
    parse_args->Add("-eno_order",&eno_order,"eno_order","eno order");
    parse_args->Add("-rk_order",&rk_order,"rk_order","runge kutta order");
    parse_args->Add("-cfl",&cfl_number,"CFL","cfl number");
    parse_args->Add("-timesplit",&timesplit,"split time stepping into an explicit advection part, and an implicit non-advection part");
    parse_args->Add("-implicit_rk",&implicit_rk,"perform runge kutta on the implicit part");
    parse_args->Add("-exact",&exact,"output a fully-explicit sim to (output_dir)_exact");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    implicit_rk=implicit_rk && !exact;

    fluids_parameters.number_particles_per_cell=16;
    fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;
    fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;fluids_parameters.write_debug_data=true;
    fluids_parameters.write_ghost_values=true;fluids_parameters.write_removed_positive_particles=true;fluids_parameters.write_removed_negative_particles=true;
    fluids_parameters.delete_fluid_inside_objects=true;
    fluids_parameters.incompressible_iterations=40;
    fluids_parameters.use_removed_positive_particles=true;fluids_parameters.use_removed_negative_particles=true;
    fluids_parameters.second_order_cut_cell_method=false;
    fluids_parameters.store_particle_ids=true;
    fluids_parameters.use_vorticity_confinement=false;
    fluids_parameters.use_vorticity_confinement_fuel=false;

    //grid
    fluids_parameters.grid->Initialize(TV_INT()+20*resolution+1,RANGE<TV>::Centered_Box()*(T).5);
    *fluids_parameters.grid=fluids_parameters.grid->Get_MAC_Grid_At_Regular_Positions();
    fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;
    //time
    initial_time=0.;last_frame=500;frame_rate=133333;
    fluids_parameters.cfl=cfl_number;;
    //custom stuff . . . 
    fluids_parameters.compressible_eos = new EOS_GAMMA<T>;
    if(eno_scheme==1) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,false,false);
    else if(eno_scheme==2) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,true,false);
    else fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,true,true);
    fluids_parameters.compressible_conservation_method->Set_Order(eno_order);
    fluids_parameters.compressible_conservation_method->Save_Fluxes();
    fluids_parameters.compressible_conservation_method->Scale_Outgoing_Fluxes_To_Clamp_Variable(true,0,(T)1e-5);
    fluids_parameters.compressible_rungekutta_order=rk_order;
    //fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_RF<T,T_GRID::dimension+2>;
    fluids_parameters.compressible_timesplit=timesplit;
    fluids_parameters.compressible_perform_rungekutta_for_implicit_part=implicit_rk;

    solids_parameters.triangle_collision_parameters.perform_self_collision=false;

    if(test_number==1||test_number==3) fluids_parameters.density=(T)1e3;
    else if(test_number==2||test_number==4) fluids_parameters.density=(T)10;
        
    if(timesplit) output_directory=STRING_UTILITIES::string_sprintf("Sod_ST_Drop/Test_%d__Resolution_%d_semiimplicit",test_number,(fluids_parameters.grid->counts.x));
    else output_directory=STRING_UTILITIES::string_sprintf("Sod_ST_Drop/Test_%d__Resolution_%d_explicit",test_number,(fluids_parameters.grid->counts.x));
    if(eno_scheme==2) output_directory+="_density_weighted";
    else if(eno_scheme==3) output_directory+="_velocity_weighted";
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    //set custom boundary
    fluids_parameters.compressible_boundary= new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<TV>(fluids_parameters.euler,T_FACE_VECTOR((T).125,(T).125),T_FACE_VECTOR((T).1,(T).1),
        TV_FACE_VECTOR(),(T).5,VECTOR_UTILITIES::Complement(fluids_parameters.domain_walls));
    fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() PHYSBAM_OVERRIDE
{
    if(test_number==1||test_number==2){
        for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) 
            fluid_collection.incompressible_fluid_collection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=(T)0;
        for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next())
            fluids_parameters.incompressible->projection.p(iterator.Cell_Index())=(T)98066.5;}
    else if(test_number==3||test_number==4){
        for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next())
            fluid_collection.incompressible_fluid_collection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=(T)100;
        for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next())
            fluids_parameters.incompressible->projection.p(iterator.Cell_Index())=(T)1e5;}
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=*fluids_parameters.grid;
    ARRAY<T,VECTOR<int,1> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    ARRAY<bool,VECTOR<int,1> > psi_cut_out(grid.Domain_Indices(),true);
    T incompressible_left_boundary=(T)-.1,incompressible_right_boundary=(T).1;
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        VECTOR<int,1> cell_index=iterator.Cell_Index();VECTOR<T,1> X=iterator.Location();
        if(X.x<=((incompressible_left_boundary+incompressible_right_boundary)*(T).5)) phi(cell_index)=incompressible_left_boundary-X.x;
        else phi(cell_index)=X.x-incompressible_right_boundary;
        if(phi(cell_index)>0) psi_cut_out(cell_index)=true;}
    fluids_parameters.euler->Set_Up_Cut_Out_Grid(psi_cut_out);
}
//#####################################################################
// Function Intialize_Euler_State
//#####################################################################
void Initialize_Euler_State() PHYSBAM_OVERRIDE
{
    T_GRID& grid=fluids_parameters.euler->grid;
    ARRAY<VECTOR<T,3> ,VECTOR<int,1> >& U=fluids_parameters.euler->U;
    EOS_GAMMA<T> *tmp_eos = dynamic_cast<EOS_GAMMA<T>*>(fluids_parameters.euler->eos);

    //initialize grid variables
    //1 == density, 2 == momentum, 3 == total energy
    if(test_number==1||test_number==2){
        for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            VECTOR<int,1> cell_index=iterator.Cell_Index();
            T rho=(T)0.,u=(T)0.,p=(T)0.;
            if(iterator.Location().x<(T)-.4){rho=(T)2.124;p=(T)148407.3;} else{rho=(T)1.58317;p=(T)98066.5;}
            U(cell_index)(1)=rho;U(cell_index)(2)=rho*u;U(cell_index)(3)=rho*(tmp_eos->e_From_p_And_rho(p,rho)+sqr(u)/(T)2.);}}
    else if(test_number==3||test_number==4){
        for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            VECTOR<int,1> cell_index=iterator.Cell_Index();
            T rho=(T)1.226,u=(T)0.,p=(T)1e5;
            U(cell_index)(1)=rho;U(cell_index)(2)=rho*u;U(cell_index)(3)=rho*(tmp_eos->e_From_p_And_rho(p,rho)+sqr(u)/(T)2.);}}
}
//#####################################################################
};
}
#endif

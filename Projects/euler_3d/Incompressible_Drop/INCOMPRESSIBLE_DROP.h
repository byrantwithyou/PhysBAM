//#####################################################################
// Copyright 2007 Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE_DROP
//#####################################################################
#ifndef __INCOMPRESSIBLE_DROP__
#define __INCOMPRESSIBLE_DROP__

#include <fstream>
#include <iostream>
#include "math.h"

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>

namespace PhysBAM{

template<class T_input>
class INCOMPRESSIBLE_DROP:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
public:
    typedef T_input T;typedef VECTOR<T,3> TV;typedef GRID<TV> T_GRID;typedef VECTOR<int,3> TV_INT;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
    typedef VECTOR<T,2*T_GRID::dimension> T_FACE_VECTOR;typedef VECTOR<TV,2*T_GRID::dimension> TV_FACE_VECTOR;

    using BASE::initial_time;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::stream_type;
    using BASE::parse_args;using BASE::test_number;using BASE::resolution;

    int eno_scheme;
    int eno_order;
    int rk_order;
    T cfl_number;
    bool timesplit;
    bool exact;

    INCOMPRESSIBLE_DROP(const STREAM_TYPE stream_type)
        :BASE(stream_type,1,fluids_parameters.COMPRESSIBLE),eno_scheme(1),eno_order(2),rk_order(3),cfl_number((T).6),timesplit(false),exact(false)
    {
    }
    
    virtual ~INCOMPRESSIBLE_DROP() {}

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
    parse_args->Add("-exact",&exact,"output a fully-explicit sim to (output_dir)_exact");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    timesplit=timesplit && !exact;

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
    fluids_parameters.gravity=false;

    //grid
    int cells=resolution;
    fluids_parameters.grid->Initialize(4*cells,cells,cells,RANGE<VECTOR<T,3> >((T)-2.,(T)2.,(T)-.5,(T).5,(T)-.5,(T).5));
    *fluids_parameters.grid=fluids_parameters.grid->Get_MAC_Grid_At_Regular_Positions();
    fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;fluids_parameters.domain_walls[1][0]=true;
    fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[2][0]=true;fluids_parameters.domain_walls[2][1]=true;
    //time
    initial_time=(T)0.;last_frame=400;frame_rate=(T)80.;
    fluids_parameters.cfl=cfl_number;
    //custom stuff . . . 
    fluids_parameters.compressible_eos = new EOS_GAMMA<T>;
    if(eno_scheme==1) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,false,false);
    else if(eno_scheme==2) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,true,false);
    else fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,true,true);
    fluids_parameters.compressible_conservation_method->Set_Order(eno_order);
    fluids_parameters.compressible_conservation_method->Save_Fluxes();
    fluids_parameters.compressible_conservation_method->Scale_Outgoing_Fluxes_To_Clamp_Variable(true,0,(T)1e-5);
    fluids_parameters.compressible_rungekutta_order=rk_order;
    fluids_parameters.compressible_conservation_method->Set_Order(eno_order);
    fluids_parameters.compressible_timesplit=timesplit;

    fluids_parameters.density=(T)10;

    if(timesplit) output_directory=STRING_UTILITIES::string_sprintf("Incompressible_Drop/Test_%d__Resolution_%d_%d_%d_semiimplicit",test_number,(fluids_parameters.grid->counts.x),
        (fluids_parameters.grid->counts.y),(fluids_parameters.grid->counts.z));
    else output_directory=STRING_UTILITIES::string_sprintf("Incompressible_Drop/Test_%d__Resolution_%d_%d_%d_explicit",test_number,(fluids_parameters.grid->counts.x),(fluids_parameters.grid->counts.y),
        (fluids_parameters.grid->counts.z));
    if(eno_scheme==2) output_directory+="_density_weighted";
    else if(eno_scheme==3) output_directory+="_velocity_weighted";
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Intialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    //set custom boundary
    VECTOR<VECTOR<bool,2>,T_GRID::dimension> valid_wall;
    for(int axis=0;axis<T_GRID::dimension;axis++) for(int axis_side=0;axis_side<2;axis_side++)
        valid_wall[axis][axis_side]=(fluids_parameters.mpi_grid?!fluids_parameters.mpi_grid->Neighbor(axis,axis_side):true) && !fluids_parameters.domain_walls[axis][axis_side];
    fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<T_GRID>(fluids_parameters.euler,
        T_FACE_VECTOR((T)1.58317,(T)1.58317,(T)1.58317,(T)1.58317,(T)1.58317,(T)1.58317),T_FACE_VECTOR((T).5,(T).5,(T).5,(T).5,(T).5,(T).5),TV_FACE_VECTOR(),(T).5,valid_wall);
    fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() PHYSBAM_OVERRIDE
{
    for(UNIFORM_GRID_ITERATOR_FACE<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next())
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=(T)0;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next())
        fluids_parameters.incompressible->projection.p(iterator.Cell_Index())=(T).5;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    T_GRID& grid=*fluids_parameters.grid;
    ARRAY<T,VECTOR<int,3> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    ARRAY<bool,VECTOR<int,3> > psi_cut_out(grid.Domain_Indices(),true);
    TV drop_center((T)-1,(T)0,(T)0);T drop_radius=(T).2;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){
        VECTOR<int,3> cell_index=iterator.Cell_Index();VECTOR<T,3> X=iterator.Location();
        phi(cell_index)=VECTOR<T,3>(X-drop_center).Magnitude()-drop_radius;
        if(phi(cell_index)>0) psi_cut_out(cell_index)=true;}
    fluids_parameters.euler->Set_Up_Cut_Out_Grid(psi_cut_out);
}
//#####################################################################
// Function Intialize_Euler_State
//#####################################################################
void Initialize_Euler_State() PHYSBAM_OVERRIDE
{
    T_GRID& grid=fluids_parameters.euler->grid;
    ARRAY<VECTOR<T,5> ,VECTOR<int,3> >& U=fluids_parameters.euler->U;
    EOS_GAMMA<T> *tmp_eos = dynamic_cast<EOS_GAMMA<T>*>(fluids_parameters.euler->eos);

    //initialize grid variables
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(fluids_parameters.euler->grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        T rho=(T)0.,u_vel=(T)0.,v_vel=(T)0.,w_vel=(T)0.,p=(T)0.;
        if(grid.X(cell_index).x<=-1.5){rho=(T)2.124;p=(T)1.5;}
        else{rho=(T)1.58317;p=(T).5;}
        U(cell_index)(0)=rho;U(cell_index)(0)=rho*u_vel;U(cell_index)(1)=rho*v_vel;U(cell_index)(2)=rho*w_vel;
        U(cell_index)(3)=rho*(tmp_eos->e_From_p_And_rho(p,rho)+(sqr(u_vel)+sqr(v_vel)+sqr(w_vel))/(T)2.);}
}
//#####################################################################
};
}
#endif

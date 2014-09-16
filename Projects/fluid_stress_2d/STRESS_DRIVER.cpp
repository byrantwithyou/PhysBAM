//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Extrapolation/EXTRAPOLATION_HIGHER_ORDER_POLY.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Grids_Uniform_Advection/ADVECTION_HAMILTON_JACOBI_ENO.h>
#include <Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <Tools/Grids_Uniform_Interpolation/QUADRATIC_INTERPOLATION_UNIFORM.h>
#include <Tools/Krylov_Solvers/GMRES.h>
#include <Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Log/DEBUG_SUBSTEPS.h>
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Tools/Vectors/VECTOR_UTILITIES.h>
#include <Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <Geometry/Finite_Elements/CELL_MANAGER_COLOR.h>
#include <Geometry/Finite_Elements/INTERFACE_STOKES_MULTIGRID.h>
#include <Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_COLOR.h>
#include <Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_VECTOR_COLOR.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <Geometry/Grids_Uniform_Computations/REINITIALIZATION.h>
#include <Geometry/Level_Sets/EXTRAPOLATION_HIGHER_ORDER.h>
#include "STRESS_DRIVER.h"
#include "STRESS_EXAMPLE.h"
#include <boost/function.hpp>
using namespace PhysBAM;
namespace{
    template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
    {
        ((STRESS_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
    }
};
//#####################################################################
// Constructor
//#####################################################################
template<class TV> STRESS_DRIVER<TV>::
STRESS_DRIVER(STRESS_EXAMPLE<TV>& example)
    :example(example)
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> STRESS_DRIVER<TV>::
~STRESS_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
}
//#####################################################################
// Execute_Main_Program
//#####################################################################
template<class TV> void STRESS_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void STRESS_DRIVER<TV>::
Initialize()
{
    LOG::cout<<std::setprecision(16)<<std::endl;
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.substeps_delay_frame<0?example.write_substeps_level:-1);

    // setup time
    current_frame=example.restart;
    output_number=current_frame;
    example.time=time=example.time_steps_per_frame*current_frame*example.dt;

    example.levelset.phi.Resize(example.grid.Domain_Indices(example.number_of_ghost_cells));
    example.face_velocities.Resize(example.grid,example.number_of_ghost_cells);
    example.prev_face_velocities.Resize(example.grid,example.number_of_ghost_cells);
    example.bc_phi.Resize(example.grid.Domain_Indices(example.number_of_ghost_cells));
    example.polymer_stress.Resize(example.grid.Domain_Indices(example.number_of_ghost_cells));
    example.prev_polymer_stress.Resize(example.grid.Domain_Indices(example.number_of_ghost_cells));

    example.Initialize();
    if(example.restart){
        example.Read_Output_Files(example.restart);}
    else{
        example.Get_Velocities(time);
        example.Get_Initial_Polymer_Stresses();}

    if(!example.restart) Write_Output_Files(0);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after init",0,1);
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class TV> void STRESS_DRIVER<TV>::
Advance_One_Time_Step(bool first_step)
{
    T dt=example.dt;
    example.Begin_Time_Step(time);
    example.Get_Velocities(time+dt);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before stress evolution",0,1);
    Update_Polymer_Stress(dt);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after stress evolution",0,1);
    example.time=time+=dt;
    example.End_Time_Step(time);
}
//#####################################################################
// Function Advection_And_BDF
//#####################################################################
template<class TV> void STRESS_DRIVER<TV>::
Advection_And_BDF(T dt,bool first_step)
{
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T,AVERAGING_UNIFORM<TV>,QUADRATIC_INTERPOLATION_UNIFORM<TV,T> > quadratic_advection;
    BOUNDARY_MAC_GRID_PERIODIC<TV,T> boundary;
    FACE_LOOKUP_UNIFORM<TV> lookup_face_velocities(example.face_velocities),lookup_prev_face_velocities(example.prev_face_velocities);
    PHYSBAM_DEBUG_ONLY(Assert_Advection_CFL(example.face_velocities,dt));
    if(!first_step){
        PHYSBAM_DEBUG_ONLY(Assert_Advection_CFL(example.prev_face_velocities,dt));
        ARRAY<T,FACE_INDEX<TV::dimension> > temp(example.grid,example.number_of_ghost_cells);
        quadratic_advection.Update_Advection_Equation_Face_Lookup(example.grid,temp,lookup_prev_face_velocities,lookup_prev_face_velocities,boundary,2*dt,time+dt);
        quadratic_advection.Update_Advection_Equation_Face_Lookup(example.grid,example.prev_face_velocities,lookup_face_velocities,lookup_face_velocities,boundary,dt,time+dt);
        example.prev_face_velocities.Copy((T)2/(T)1.5,example.prev_face_velocities,-(T).5/(T)1.5,temp);}
    else quadratic_advection.Update_Advection_Equation_Face_Lookup(example.grid,example.prev_face_velocities,lookup_face_velocities,lookup_face_velocities,boundary,dt,time+dt);
    example.prev_face_velocities.Exchange(example.face_velocities);
}
//#####################################################################
// Function Assert_Advection_CFL
//#####################################################################
template<class TV> void STRESS_DRIVER<TV>::
Assert_Advection_CFL(const ARRAY<T,FACE_INDEX<TV::m> >& u,T dt) const
{
    for(FACE_ITERATOR<TV> it(example.grid);it.Valid();it.Next()){
        PHYSBAM_ASSERT(u(it.Full_Index())*dt<example.grid.dX(it.Axis())*example.number_of_ghost_cells);}
}
//#####################################################################
// Function Update_Polymer_Stress
//#####################################################################
template<class TV> void STRESS_DRIVER<TV>::
Update_Polymer_Stress(T dt)
{
}
//#####################################################################
// Function Extrapolate_Velocity
//#####################################################################
template<class TV> void STRESS_DRIVER<TV>::
Extrapolate_Stress(const ARRAY<T,FACE_INDEX<TV::m> >& u,ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT>& S)
{
    for(CELL_ITERATOR<TV> it(example.grid,example.number_of_ghost_cells,GRID<TV>::GHOST_REGION);it.Valid();it.Next())
        S(it.index)=SYMMETRIC_MATRIX<T,TV::m>()+1e20;

    EXTRAPOLATION_HIGHER_ORDER<TV,SYMMETRIC_MATRIX<T,TV::m> > eho(example.grid,example.levelset,example.number_of_ghost_cells*10,3,example.number_of_ghost_cells);
    eho.periodic=true;

    eho.Extrapolate_Cell([&](const TV_INT& index){return example.levelset.phi(index)<=0;},S);

    example.boundary_stress.Apply_Boundary_Condition(example.grid,S,time+example.dt);
    example.boundary_stress.Fill_Ghost_Cells(example.grid,S,S,0,example.number_of_ghost_cells);
}
//#####################################################################
// Simulate_To_Frame
//#####################################################################
template<class TV> void STRESS_DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    for(;current_frame<frame;current_frame++){
        LOG::SCOPE scope("FRAME","frame %d",current_frame+1);
        if(example.substeps_delay_frame==current_frame)
            DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);
        for(int substep=0;substep<example.time_steps_per_frame;substep++){
            LOG::SCOPE scope("SUBSTEP","substep %d",substep+1);
            Advance_One_Time_Step(current_frame==0 && substep==0);}
        Write_Output_Files(++output_number);}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void STRESS_DRIVER<TV>::
Write_Substep(const std::string& title,const int substep,const int level)
{
    if(level<=example.write_substeps_level){
        example.frame_title=title;
        LOG::cout<<"Writing substep ["<<example.frame_title<<"]: output_number="<<output_number+1<<", time="<<time<<", frame="<<current_frame<<", substep="<<substep<<std::endl;
        Write_Output_Files(++output_number);
        example.frame_title="";}
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void STRESS_DRIVER<TV>::
Write_Output_Files(const int frame)
{
    FILE_UTILITIES::Create_Directory(example.output_directory);
    FILE_UTILITIES::Create_Directory(example.output_directory+STRING_UTILITIES::string_sprintf("/%d",frame));
    FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+STRING_UTILITIES::string_sprintf("/%d/frame_title",frame),example.frame_title);
    if(frame==0)
        FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/first_frame",frame,"\n");
    example.Write_Output_Files(frame);
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/last_frame",frame,"\n");
}
//#####################################################################
namespace PhysBAM{
template class STRESS_DRIVER<VECTOR<float,2> >;
template class STRESS_DRIVER<VECTOR<float,3> >;
template class STRESS_DRIVER<VECTOR<double,2> >;
template class STRESS_DRIVER<VECTOR<double,3> >;
}

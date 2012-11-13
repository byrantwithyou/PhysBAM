//#####################################################################
// Copyright 2009-2010, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include "PROJECTION_DRIVER.h"
#include "PROJECTION_EXAMPLE.h"

using namespace PhysBAM;
namespace{
template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
{
    ((PROJECTION_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
}
};
//#####################################################################
// Initialize
//#####################################################################
template<class TV> PROJECTION_DRIVER<TV>::
PROJECTION_DRIVER(PROJECTION_EXAMPLE<TV>& example)
    :example(example)
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> PROJECTION_DRIVER<TV>::
~PROJECTION_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void PROJECTION_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void PROJECTION_DRIVER<TV>::
Initialize()
{
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);

    // setup time
    if(example.restart) current_frame=example.restart;else current_frame=example.first_frame;    
    time=example.Time_At_Frame(current_frame);

    // mpi
    if(example.mpi_grid) example.mpi_grid->Initialize(example.domain_boundary);
    example.projection.elliptic_solver->mpi_grid=example.mpi_grid;
    if(example.mpi_grid) example.boundary=new BOUNDARY_MPI<GRID<TV>,T>(example.mpi_grid,example.boundary_scalar);
    else example.boundary=&example.boundary_scalar;

    //threading
    example.projection.elliptic_solver->thread_queue=example.thread_queue;

    // setup grids and velocities
    example.projection.Initialize_Grid(example.mac_grid);
    example.face_velocities.Resize(example.mac_grid);
    example.Initialize_Fields();

    // setup laplace
    example.projection.elliptic_solver->Set_Relative_Tolerance(1e-11);
    example.projection.elliptic_solver->pcg.Set_Maximum_Iterations(200);
    example.projection.elliptic_solver->pcg.evolution_solver_type=krylov_solver_cg;
    example.projection.elliptic_solver->pcg.cg_restart_iterations=40;

    if(example.restart) example.Read_Output_Files(example.restart);
    
    // setup domain boundaries
    VECTOR<VECTOR<bool,2>,TV::dimension> constant_extrapolation;constant_extrapolation.Fill(VECTOR<bool,2>::Constant_Vector(true));
    example.boundary->Set_Constant_Extrapolation(constant_extrapolation);
    example.Set_Boundary_Conditions(time); // get so CFL is correct

    if(!example.restart) Write_Output_Files(example.first_frame);
    output_number=example.first_frame;
}
//#####################################################################
// Advance_One_Time_Step_Implicit_Part
//#####################################################################
template<class TV> void PROJECTION_DRIVER<TV>::
Advance_One_Time_Step_Implicit_Part(ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,const T dt,const T time)
{
    example.boundary->Apply_Boundary_Condition_Face(example.mac_grid,face_velocities,time+dt);
    example.projection.Make_Divergence_Free(face_velocities,dt,time);
}
//#####################################################################
// Project
//#####################################################################
template<class TV> void PROJECTION_DRIVER<TV>::
Project(const T dt,const T time)
{
    LOG::Time("Start");
    example.Set_Boundary_Conditions(time+dt);
    example.projection.p*=dt; // rescale pressure for guess
    LOG::Time("Project");
    Advance_One_Time_Step_Implicit_Part(example.face_velocities,dt,time);
    LOG::Time("End");
    example.projection.p*=(1/dt); // unscale pressure
}
//#####################################################################
// Advance_To_Target_Time
//#####################################################################
template<class TV> void PROJECTION_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    T dt=target_time-time;
    Project(dt,time);
    time+=dt;
}
//#####################################################################
// Simulate_To_Frame
//#####################################################################
template<class TV> void PROJECTION_DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    while(current_frame<frame){
        LOG::SCOPE scope("FRAME","Frame %d",current_frame+1);
        Advance_To_Target_Time(example.Time_At_Frame(current_frame+1));
        Write_Output_Files(++output_number);
        current_frame++;}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void PROJECTION_DRIVER<TV>::
Write_Substep(const std::string& title,const int substep,const int level)
{
    if(level<=example.write_substeps_level){
        example.frame_title=title;
        LOG::cout<<"Writing substep ["<<example.frame_title<<"]: output_number="<<output_number+1<<", time="<<time<<", frame="<<current_frame<<", substep="<<substep<<std::endl;
        Write_Output_Files(++output_number);example.frame_title="";}
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void PROJECTION_DRIVER<TV>::
Write_Output_Files(const int frame)
{
    FILE_UTILITIES::Create_Directory(example.output_directory);
    FILE_UTILITIES::Create_Directory(example.output_directory+STRING_UTILITIES::string_sprintf("/%d",frame));
    FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+STRING_UTILITIES::string_sprintf("/%d/frame_title",frame),example.frame_title);
    if(frame==example.first_frame) 
        FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/first_frame",frame,"\n");
    example.Write_Output_Files(frame);
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/last_frame",frame,"\n");
}
//#####################################################################
template class PROJECTION_DRIVER<VECTOR<float,1> >;
template class PROJECTION_DRIVER<VECTOR<float,2> >;
template class PROJECTION_DRIVER<VECTOR<float,3> >;
template class PROJECTION_DRIVER<VECTOR<double,1> >;
template class PROJECTION_DRIVER<VECTOR<double,2> >;
template class PROJECTION_DRIVER<VECTOR<double,3> >;

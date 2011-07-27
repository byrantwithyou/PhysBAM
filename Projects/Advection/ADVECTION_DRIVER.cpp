//#####################################################################
// Copyright 2009-2010, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_THREADED.h>
#include "ADVECTION_DRIVER.h"
#include "ADVECTION_EXAMPLE.h"

using namespace PhysBAM;
namespace{
template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
{
    ((ADVECTION_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
}
};
//#####################################################################
// Initialize
//#####################################################################
template<class TV> ADVECTION_DRIVER<TV>::
ADVECTION_DRIVER(ADVECTION_EXAMPLE<TV>& example)
    :example(example)
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> ADVECTION_DRIVER<TV>::
~ADVECTION_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void ADVECTION_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void ADVECTION_DRIVER<TV>::
Initialize()
{
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);

    // setup time
    if(example.restart) current_frame=example.restart;else current_frame=example.first_frame;    
    time=example.Time_At_Frame(current_frame);

    // mpi
    if(example.mpi_grid) example.boundary=new BOUNDARY_MPI<GRID<TV> >(example.mpi_grid,example.boundary_scalar);
    else if(example.thread_queue) example.boundary=new BOUNDARY_THREADED<GRID<TV> >(*example.thread_queue,example.boundary_scalar);
    else example.boundary=&example.boundary_scalar;

    // setup grids and velocities
    example.face_velocities.Resize(example.mac_grid);
    example.density.Resize(example.mac_grid.Domain_Indices(3));
    example.Initialize_Fields();

    if(example.restart) example.Read_Output_Files(example.restart);
    
    // setup domain boundaries
    VECTOR<VECTOR<bool,2>,TV::dimension> constant_extrapolation;constant_extrapolation.Fill(VECTOR<bool,2>::Constant_Vector(true));
    example.boundary->Set_Constant_Extrapolation(constant_extrapolation);

    if(!example.restart) Write_Output_Files(example.first_frame);
    output_number=example.first_frame;
}
//#####################################################################
// Scalar_Advance
//#####################################################################
template<class TV> void ADVECTION_DRIVER<TV>::
Scalar_Advance(const T dt,const T time)
{
    ARRAY<T,TV_INT> density_ghost(example.mac_grid.Domain_Indices(3),false);
    example.boundary->Fill_Ghost_Cells(example.mac_grid,example.density,density_ghost,dt,time,3);
    example.advection_scalar.Update_Advection_Equation_Cell(example.mac_grid,example.density,density_ghost,example.face_velocities,*example.boundary,dt,time);
}
//#####################################################################
// Advance_To_Target_Time
//#####################################################################
template<class TV> void ADVECTION_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;for(int substep=1;!done;substep++){
        LOG::Time("Substep");
        T dt=example.cfl*example.CFL(example.face_velocities);
        if(example.mpi_grid) example.mpi_grid->Synchronize_Dt(dt);
        if(time+dt>=target_time){dt=target_time-time;done=true;}
        else if(time+2*dt>=target_time){dt=.5*(target_time-time);}
        Scalar_Advance(dt,time);
        time+=dt;}
}
//#####################################################################
// Simulate_To_Frame
//#####################################################################
template<class TV> void ADVECTION_DRIVER<TV>::
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
template<class TV> void ADVECTION_DRIVER<TV>::
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
template<class TV> void ADVECTION_DRIVER<TV>::
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
template class ADVECTION_DRIVER<VECTOR<float,1> >;
template class ADVECTION_DRIVER<VECTOR<float,2> >;
template class ADVECTION_DRIVER<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ADVECTION_DRIVER<VECTOR<double,1> >;
template class ADVECTION_DRIVER<VECTOR<double,2> >;
template class ADVECTION_DRIVER<VECTOR<double,3> >;
#endif

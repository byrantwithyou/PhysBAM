//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Tools/Ordinary_Differential_Equations/DRIVER.h>
#include <Tools/Ordinary_Differential_Equations/EXAMPLE.h>
using namespace PhysBAM;
//#####################################################################
// DRIVER
//#####################################################################
template<class TV> DRIVER<TV>::
DRIVER(EXAMPLE<TV>& example)
    :time(0),example(example),current_frame(0),output_number(0)
{
    DEBUG_SUBSTEPS::write_substeps_level=example.write_substeps_level;
    DEBUG_SUBSTEPS::writer=[=](const std::string& title){Write_Substep(title);};
}
//#####################################################################
// ~DRIVER
//#####################################################################
template<class TV> DRIVER<TV>::
~DRIVER()
{
    DEBUG_SUBSTEPS::writer=0;
}
//#####################################################################
// Execute_Main_Program
//#####################################################################
template<class TV> void DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void DRIVER<TV>::
Initialize()
{
    // setup time
    example.Setup_Log();
    if(example.auto_restart){Read_Last_Frame();example.restart=true;}
    if(example.restart){current_frame=example.restart_frame;Read_Time(current_frame);}
    else current_frame=0;
    output_number=current_frame;
    time=example.Time_At_Frame(current_frame);
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void DRIVER<TV>::
Write_Substep(const std::string& title)
{
    example.frame_title=title;
    LOG::cout<<"Writing substep ["<<example.frame_title<<"]: output_number="<<output_number+1<<", time="<<time<<", frame="<<current_frame<<std::endl;
    Write_Output_Files(++output_number);
    example.frame_title="";
}
//#####################################################################
// Simulate_To_Frame
//#####################################################################
template<class TV> void DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    while(current_frame<frame){
        LOG::SCOPE scope("FRAME","Frame %d",current_frame+1);
        Advance_To_Target_Time(example.Time_At_Frame(current_frame+1));
        Write_Output_Files(++output_number);
        current_frame++;}
}
//#####################################################################
// Function Read_Time
//#####################################################################
template<class TV> void DRIVER<TV>::
Read_Time(const int frame)
{
    time=example.Time_At_Frame(frame);
    std::string filename=LOG::sprintf("%s/%d/time",example.output_directory.c_str(),frame);
    if(File_Exists(filename)){
        T corrected_time;
        Read_From_File(filename,corrected_time);
        if(abs(time-corrected_time)>(T)1e-4*abs(time)){ // only adjust time if significantly different from default in order to get deterministic restarts
            time=corrected_time;
            // adjust initial time so that Simulate_To_Frame() returns correct time (essential when writing substeps)
        }}
}
//#####################################################################
// Function Read_Last_Frame
//#####################################################################
template<class TV> void DRIVER<TV>::
Read_Last_Frame()
{
    std::string filename=LOG::sprintf("%s/common/last_frame",example.output_directory.c_str());
    if(File_Exists(filename))
        Read_From_Text_File(filename, example.restart_frame);
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void DRIVER<TV>::
Write_Output_Files(const int frame)
{
    Create_Directory(example.output_directory);
    Create_Directory(example.output_directory+LOG::sprintf("/%d",frame));
    Create_Directory(example.output_directory+"/common");
    example.Write_Output_Files(frame);
    Write_Time(frame);
    Write_Last_Frame(frame);
    Write_To_Text_File(example.output_directory+LOG::sprintf("/%d/frame_title",frame),example.frame_title);
}
//#####################################################################
namespace PhysBAM{
template class DRIVER<VECTOR<float,1> >;
template class DRIVER<VECTOR<float,2> >;
template class DRIVER<VECTOR<float,3> >;
template class DRIVER<VECTOR<double,1> >;
template class DRIVER<VECTOR<double,2> >;
template class DRIVER<VECTOR<double,3> >;
}

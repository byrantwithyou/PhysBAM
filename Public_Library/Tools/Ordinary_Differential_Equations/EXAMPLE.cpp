//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Log/SCOPE.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Ordinary_Differential_Equations/EXAMPLE.h>
#include <Tools/Parallel_Computation/MPI_WORLD.h>
#include <Tools/Parsing/PARSE_ARGS.h>
using namespace PhysBAM;
//#####################################################################
// EXAMPLE
//#####################################################################
template<class TV> EXAMPLE<TV>::
EXAMPLE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :stream_type(stream_type_input),last_frame(120),frame_rate(24),frame_title(""),
    write_substeps_level(-1),
    viewer_dir("output"),data_directory("../../Public_Data"),auto_restart(false),restart(0),
    write_output_files(true),write_frame_title(true),abort_when_dt_below(0),
    mpi_world(0),need_finish_logging(true),test_number(0),
    substeps_delay_frame(-1),substeps_delay_level(-1),use_test_output(false),test_output_prefix(""),
    opt_all_verbose(false),opt_query_output(false),opt_nolog(false),opt_verbosity(1<<30)
{
    stored_args=parse_args.Print_Arguments();
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    PROCESS_UTILITIES::Set_Backtrace(true);

    if(const char* d=getenv("PHYSBAM_DATA_DIRECTORY")) data_directory=d;
    parse_args.Extra_Optional(&test_number,"example number","example number to run");
    parse_args.Add("-dt",&fixed_dt,"size","fix the time step size to this value.");
    parse_args.Add("-framerate",&frame_rate,&user_frame_rate,"rate","frame rate");
    parse_args.Add("-max_dt",&max_dt,"size","fix the time step size to be no larger than this value.");
    parse_args.Add("-delay_substeps",&substeps_delay_frame,"frame","delay substeps until later frame");
    parse_args.Add("-last_frame",&last_frame,&user_last_frame,"frame","last frame");
    parse_args.Add("-restart",&restart,"frame","restart frame");
    parse_args.Add("-substeps",&substeps_delay_level,"level","substep output level");
    parse_args.Add("-v",&opt_verbosity,"level","verbosity level");
    parse_args.Add("-auto_restart",&auto_restart,"restart from last_frame");
    parse_args.Add("-nolog",&opt_nolog,"disable log.txt");
    parse_args.Add("-query_output",&opt_query_output,"print the output directory and exit");
    parse_args.Add("-d",&data_directory,"dir","data directory");
    parse_args.Add("-o",&viewer_dir.output_directory,&user_output_directory,"dir","output directory");
    parse_args.Add("-test_output_prefix",&test_output_prefix,&use_test_output,"","prefix to use for test output");
    parse_args.Add("-m",&m,"scale","length unit");
    parse_args.Add("-s",&s,"scale","time unit");
    parse_args.Add("-kg",&kg,"scale","mass unit");
    parse_args.Add("-all_verbose",&opt_all_verbose,"all mpi processes write to stdout (not just the first)");
    parse_args.Parse(true);

    if(mpi_world && !opt_all_verbose && mpi_world->initialized && mpi_world->rank) opt_verbosity=0;
    LOG::Initialize_Logging(opt_verbosity<10,false,opt_verbosity,!opt_nolog);
    fixed_dt*=s;
    frame_rate/=s;
    max_dt*=s;

    if(substeps_delay_frame==-1) Set_Write_Substeps_Level(substeps_delay_level);
    else Set_Write_Substeps_Level(-1);
}
//#####################################################################
// ~EXAMPLE
//#####################################################################
template<class TV> EXAMPLE<TV>::
~EXAMPLE()
{
    delete mpi_world;
    if(need_finish_logging) LOG::Finish_Logging();
}
template<class TV> void EXAMPLE<TV>::
Setup_Log() const
{
    if(opt_query_output){LOG::cout<<viewer_dir.output_directory;exit(0);}
    LOG::cout<<stored_args<<std::endl;
}
template<class TV> typename TV::SCALAR EXAMPLE<TV>::
Time_At_Frame(const int frame) const
{
    return frame/frame_rate;
}
template<class TV> bool EXAMPLE<TV>::
Clamp_Time_Step_With_Target_Time(const T time,const T target_time,T& dt)
{
    if(limit_dt) limit_dt(dt,time);
    if(max_dt && dt>max_dt) dt=max_dt;
    if(min_dt && dt<min_dt) dt=min_dt;
    if(fixed_dt) dt=fixed_dt;
    if(time+dt>=target_time){dt=target_time-time;return true;}
    else if(time+2*dt>=target_time) dt=min(dt,(T).51*(target_time-time));
    return false;
}
template<class TV> void EXAMPLE<TV>::
Set_Write_Substeps_Level(const int level)
{
    write_substeps_level=level;
    DEBUG_SUBSTEPS::write_substeps_level=level;
}
template<class TV> void EXAMPLE<TV>::
Write_Frame_Title(const int frame) const
{
    if(write_frame_title) Write_To_Text_File(viewer_dir.current_directory+"/frame_title",frame_title);
}
//#####################################################################
// Function Log_Parameters
//#####################################################################
template<class TV> void EXAMPLE<TV>::
Log_Parameters() const
{
    LOG::SCOPE scope("EXAMPLE parameters");
    LOG::cout<<"last_frame="<<last_frame<<std::endl;
    LOG::cout<<"frame_rate="<<frame_rate<<std::endl;
    LOG::cout<<"auto_restart="<<auto_restart<<std::endl;
    LOG::cout<<"restart="<<restart<<std::endl;
    LOG::cout<<"write_output_files="<<write_output_files<<std::endl;
    LOG::cout<<"write_frame_title="<<write_frame_title<<std::endl;
    LOG::cout<<"write_substeps_level="<<write_substeps_level<<std::endl;
    LOG::cout<<"output_directory="<<viewer_dir.output_directory<<std::endl;
    LOG::cout<<"data_directory="<<data_directory<<std::endl;
    LOG::cout<<"frame_title="<<frame_title<<std::endl;
    LOG::cout<<"abort_when_dt_below="<<abort_when_dt_below<<std::endl;
}
//#####################################################################
namespace PhysBAM{
template class EXAMPLE<VECTOR<float,1> >;
template class EXAMPLE<VECTOR<float,2> >;
template class EXAMPLE<VECTOR<float,3> >;
template class EXAMPLE<VECTOR<double,1> >;
template class EXAMPLE<VECTOR<double,2> >;
template class EXAMPLE<VECTOR<double,3> >;
}

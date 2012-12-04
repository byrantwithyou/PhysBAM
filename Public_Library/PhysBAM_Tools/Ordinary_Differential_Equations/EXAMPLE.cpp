//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Ordinary_Differential_Equations/EXAMPLE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// EXAMPLE
//#####################################################################
template<class TV> EXAMPLE<TV>::
EXAMPLE(const STREAM_TYPE stream_type_input)
    :stream_type(stream_type_input),initial_time(0),first_frame(0),last_frame(120),frame_rate(24),frame_title(""),
    write_substeps_level(-1),write_first_frame(true),write_last_frame(true),write_time(true),
    output_directory("output"),data_directory("../../Public_Data"),auto_restart(false),restart(false),
    restart_frame(0),write_output_files(true),write_frame_title(true),abort_when_dt_below(0),parse_args(0),
    mpi_world(0),want_mpi_world(false),need_finish_logging(false),test_number(0),fixed_dt(0),max_dt(0),
    substeps_delay_frame(-1),substeps_delay_level(-1),use_test_output(false),test_output_prefix(""),
    stored_data_directory(data_directory),stored_output_directory(output_directory),opt_all_verbose(false),user_dt(false),
    user_frame_rate(false),user_max_dt(false),user_first_frame(false),user_last_frame(false),user_data_directory(false),
    user_output_directory(false),opt_query_output(false),opt_nolog(false),opt_verbosity(1<<30),stored_first_frame(0),
    stored_last_frame(0),stored_dt(0),stored_frame_rate(24),stored_max_dt(0)
{
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
template<class TV> typename TV::SCALAR EXAMPLE<TV>::
Time_At_Frame(const int frame) const
{
    return initial_time+(frame-first_frame)/frame_rate;
}
template<class TV> void EXAMPLE<TV>::
Clamp_Time_Step_With_Target_Time(const T time,const T target_time,T& dt,bool& done,const T min_dt,bool* min_dt_failed)
{
    if(dt<min_dt){dt=min_dt;if(min_dt_failed) *min_dt_failed=true;}
    if(time+dt>=target_time){dt=target_time-time;done=true;}
    else if(time+2*dt>=target_time) dt=min(dt,(T).51*(target_time-time));
}
template<class TV> void EXAMPLE<TV>::
Set_Write_Substeps_Level(const int level)
{
    write_substeps_level=level;
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(level);
}
template<class TV> void EXAMPLE<TV>::
Write_Frame_Title(const int frame) const
{
    if(write_frame_title) FILE_UTILITIES::Write_To_Text_File(STRING_UTILITIES::string_sprintf("%s/%d/frame_title",output_directory.c_str(),frame),frame_title);
}
template<class TV> void EXAMPLE<TV>::
Limit_Dt(T& dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Log_Parameters
//#####################################################################
template<class TV> void EXAMPLE<TV>::
Log_Parameters() const
{
    LOG::SCOPE scope("EXAMPLE parameters");
    LOG::cout<<"initial_time="<<initial_time<<std::endl;
    LOG::cout<<"first_frame="<<first_frame<<std::endl;
    LOG::cout<<"last_frame="<<last_frame<<std::endl;
    LOG::cout<<"frame_rate="<<frame_rate<<std::endl;
    LOG::cout<<"auto_restart="<<auto_restart<<std::endl;
    LOG::cout<<"restart="<<restart<<std::endl;
    LOG::cout<<"restart_frame="<<restart_frame<<std::endl;
    LOG::cout<<"write_output_files="<<write_output_files<<std::endl;
    LOG::cout<<"write_first_frame="<<write_first_frame<<std::endl;
    LOG::cout<<"write_last_frame="<<write_last_frame<<std::endl;
    LOG::cout<<"write_time="<<write_time<<std::endl;
    LOG::cout<<"write_frame_title="<<write_frame_title<<std::endl;
    LOG::cout<<"write_substeps_level="<<write_substeps_level<<std::endl;
    LOG::cout<<"output_directory="<<output_directory<<std::endl;
    LOG::cout<<"data_directory="<<data_directory<<std::endl;
    LOG::cout<<"frame_title="<<frame_title<<std::endl;
    LOG::cout<<"abort_when_dt_below="<<abort_when_dt_below<<std::endl;
}
//#####################################################################
// Function Register_Options
//#####################################################################
template<class TV> void EXAMPLE<TV>::
Register_Options()
{
    if(!parse_args) return;
    parse_args->Extra_Optional(&test_number,"example number","example number to run");
    parse_args->Add("-dt",&stored_dt,&user_dt,"size","fix the time step size to this value.");
    parse_args->Add("-framerate",&stored_frame_rate,&user_frame_rate,"rate","frame rate");
    parse_args->Add("-max_dt",&stored_max_dt,&user_max_dt,"size","fix the time step size to be no larger than this value.");
    parse_args->Add("-delay_substeps",&substeps_delay_frame,"frame","delay substeps until later frame");
    parse_args->Add("-first_frame",&stored_first_frame,&user_first_frame,"frame","first frame");
    parse_args->Add("-last_frame",&stored_last_frame,&user_last_frame,"frame","last frame");
    parse_args->Add("-restart",&restart_frame,&restart,"frame","restart frame");
    parse_args->Add("-substeps",&substeps_delay_level,"level","substep output level");
    parse_args->Add("-v",&opt_verbosity,"level","verbosity level");
    parse_args->Add("-auto_restart",&auto_restart,"restart from last_frame");
    parse_args->Add("-nolog",&opt_nolog,"disable log.txt");
    parse_args->Add("-query_output",&opt_query_output,"print the output directory and exit");
    parse_args->Add("-d",&stored_data_directory,&user_data_directory,"dir","data directory");
    parse_args->Add("-o",&stored_output_directory,&user_output_directory,"dir","output directory");
    parse_args->Add("-test_output_prefix",&test_output_prefix,&use_test_output,"","prefix to use for test output");
    if(mpi_world) parse_args->Add("-all_verbose",&opt_all_verbose,"all mpi processes write to stdout (not just the first)");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
template<class TV> void EXAMPLE<TV>::
Parse_Options()
{
    if(!parse_args) return;
    if(mpi_world && !opt_all_verbose && mpi_world->initialized && mpi_world->rank) opt_verbosity=0;
    need_finish_logging=true;
    LOG::Initialize_Logging(opt_verbosity<10,false,opt_verbosity,!opt_nolog);
}
//#####################################################################
// Function Override_Options
//#####################################################################
template<class TV> void EXAMPLE<TV>::
Override_Options()
{
    if(user_output_directory) output_directory=stored_output_directory;
    if(user_data_directory) data_directory=stored_data_directory;
    else if(const char* d=getenv("PHYSBAM_DATA_DIRECTORY")) data_directory=d;
    if(substeps_delay_frame==-1 && substeps_delay_level!=-1) Set_Write_Substeps_Level(substeps_delay_level);
    if(user_first_frame) first_frame=stored_first_frame;
    if(user_frame_rate) frame_rate=stored_frame_rate;
    if(user_dt) fixed_dt=stored_dt;
    if(user_max_dt) max_dt=stored_max_dt;

    if(opt_query_output){LOG::cout<<output_directory;exit(0);}
    if(!opt_nolog){
        if(!restart && !auto_restart) FILE_UTILITIES::Create_Directory(output_directory);
        FILE_UTILITIES::Create_Directory(output_directory+"/common");
        LOG::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",restart);}
}
//#####################################################################
// Function Parse_Late_Options
//#####################################################################
template<class TV> void EXAMPLE<TV>::
Parse_Late_Options()
{
    if(!parse_args) return;
    if(user_last_frame) last_frame=stored_last_frame;
}
//#####################################################################
// Function Parse
//#####################################################################
template<class TV> void EXAMPLE<TV>::
Parse(PARSE_ARGS& parse_args_input)
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    PROCESS_UTILITIES::Set_Backtrace(true);

    parse_args=&parse_args_input;
    if(want_mpi_world) mpi_world=new MPI_WORLD(*parse_args);
    Register_Options();

    std::string print_args=parse_args->Print_Arguments();
    parse_args->Parse();

    Parse_Options();
    LOG::cout<<print_args<<std::endl;
    Override_Options();
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

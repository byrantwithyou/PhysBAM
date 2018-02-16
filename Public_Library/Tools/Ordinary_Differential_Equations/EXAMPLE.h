//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __EXAMPLE__
#define __EXAMPLE__
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Vectors/VECTOR.h>
#include <Grid_PDE/Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
namespace PhysBAM{

class PARSE_ARGS;
class MPI_WORLD;
template<class TV>
class EXAMPLE
{
    typedef typename TV::SCALAR T;
    enum workaround1{d=TV::m};

public:
    const STREAM_TYPE stream_type;
    T initial_time;
    int first_frame,last_frame;
    T frame_rate;
    std::string frame_title;
    int write_substeps_level;
    bool write_first_frame,write_last_frame,write_time;
    std::string output_directory,data_directory;

    bool auto_restart;
    bool restart;
    int restart_frame;
    bool write_output_files,write_frame_title;

    T abort_when_dt_below;
    MPI_WORLD* mpi_world;
    bool need_finish_logging;
    int test_number;
    bool use_default_test;
    T fixed_dt=0;
    T min_dt=0;
    T max_dt=0;
    int substeps_delay_frame;
    int substeps_delay_level;
    bool use_test_output;
    std::string test_output_prefix;
    bool opt_all_verbose,user_last_frame=false,user_output_directory=false;
    bool user_frame_rate=false;
    bool opt_query_output,opt_nolog;
    int opt_verbosity;
    T m=1,s=1,kg=1;
    std::string stored_args;
    
    EXAMPLE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    virtual ~EXAMPLE();

    std::function<void(T& dt,T time)> limit_dt;

//#####################################################################
    virtual T Time_At_Frame(const int frame) const;
    bool Clamp_Time_Step_With_Target_Time(const T time,const T target_time,T& dt);
    void Set_Write_Substeps_Level(const int level);
    void Write_Frame_Title(const int frame) const;
    void Setup_Log() const;
    virtual void Write_Output_Files(const int frame) const=0;
    virtual void Log_Parameters() const;
//#####################################################################
};
}
#endif

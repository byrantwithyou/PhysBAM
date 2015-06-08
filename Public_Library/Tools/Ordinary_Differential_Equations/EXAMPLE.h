//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __EXAMPLE__
#define __EXAMPLE__
#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Log/DEBUG_SUBSTEPS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Vectors/VECTOR.h>
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
    T fixed_dt;
    T max_dt;
    int substeps_delay_frame;
    int substeps_delay_level;
    bool use_test_output;
    std::string test_output_prefix,stored_data_directory,stored_output_directory;
    bool opt_all_verbose,user_dt,user_frame_rate,user_max_dt,user_first_frame,user_last_frame,user_data_directory,user_output_directory;
    bool opt_query_output,opt_nolog;
    int opt_verbosity,stored_first_frame,stored_last_frame;
    T stored_dt,stored_frame_rate,stored_max_dt;
    T m,s,kg;
    std::string stored_args;
    
    EXAMPLE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    virtual ~EXAMPLE();
    
//#####################################################################
    virtual T Time_At_Frame(const int frame) const;
    static void Clamp_Time_Step_With_Target_Time(const T time,const T target_time,T& dt,bool& done,const T min_dt=0,bool* min_dt_failed=0);
    void Set_Write_Substeps_Level(const int level);
    void Write_Frame_Title(const int frame) const;
    virtual void Limit_Dt(T& dt,const T time);
    virtual void Write_Output_Files(const int frame) const=0;
    virtual void Log_Parameters() const;
    virtual void After_Construction(); // Call parent first
    virtual void After_Initialization(); // Call parent last
//#####################################################################
};
}
#endif

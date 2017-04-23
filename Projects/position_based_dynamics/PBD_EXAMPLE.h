//#####################################################################
// Copyright 2015
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PBD_EXAMPLE__
#define __PBD_EXAMPLE__
#include <Core/Vectors/VECTOR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
namespace PhysBAM{

template<class TV> class DEBUG_PARTICLES;
template<class TV> class PBD_CONSTRAINTS_BASE;

template<class TV>
class PBD_EXAMPLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    STREAM_TYPE stream_type;
    DEBUG_PARTICLES<TV>& debug_particles;

    ARRAY<T> w;
    ARRAY<TV> X;
    ARRAY<TV> P;
    ARRAY<TV> V;
    ARRAY<PBD_CONSTRAINTS_BASE<TV>*> constraints;

    T initial_time;
    int last_frame;
    std::string frame_title;
    int write_substeps_level;
    int substeps_delay_frame;
    std::string output_directory;
    int restart;
    T dt,time,frame_dt,min_dt,max_dt;
    bool print_stats;
    int solver_iterations;

    bool test_diff;
    int threads;

    PBD_EXAMPLE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    PBD_EXAMPLE(const PBD_EXAMPLE&) = delete;
    virtual ~PBD_EXAMPLE();
    void operator=(const PBD_EXAMPLE&) = delete;
    
    virtual void Write_Output_Files(const int frame);
    virtual void Read_Output_Files(const int frame);
    virtual void Initialize()=0;
    virtual void Begin_Frame(const int frame)=0;
    virtual void End_Frame(const int frame)=0;
    virtual void Begin_Time_Step(const T time)=0;
    virtual void End_Time_Step(const T time)=0;

    int Add_Constraints(PBD_CONSTRAINTS_BASE<TV>& constraints);
//#####################################################################
};
}
#endif

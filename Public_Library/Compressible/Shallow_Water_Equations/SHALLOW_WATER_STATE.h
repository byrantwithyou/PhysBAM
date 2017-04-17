//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SHALLOW_WATER_STATE__
#define __SHALLOW_WATER_STATE__
#include <Core/Utilities/NONCOPYABLE.h>
#include <Grid_Tools/Grids/GRID.h>
#include <functional>
namespace PhysBAM{

template<class T> class SHALLOW_WATER;
template<class TV> class DEBUG_PARTICLES;

template<class TV>
class SHALLOW_WATER_STATE:public NONCOPYABLE
{
    STATIC_ASSERT(TV::m<3);
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef VECTOR<T,TV::m+1> T_VEC;
public:
    GRID<TV> grid;
    STREAM_TYPE stream_type;
    int ghost;

    // initialization & output
    T initial_time;
    int last_frame;
    std::string frame_title;
    int write_substeps_level;
    int substeps_delay_frame;
    std::string output_directory;
    std::string data_directory;
    std::string test_output_prefix;
    bool use_test_output;
    int restart;
    T dt,time,frame_dt,min_dt,max_dt;
    T cfl;

    ARRAY<T_VEC,TV_INT> U;
    SHALLOW_WATER<TV>& shallow_water;
    DEBUG_PARTICLES<TV>& debug_particles;

    std::function<void()> initialize;
    std::function<void(int frame)> begin_frame;
    std::function<void(int frame)> end_frame;
    std::function<void(T time)> begin_time_step;
    std::function<void(T time)> end_time_step;
    std::function<void(int frame)> write_output_files;
    std::function<void(int frame)> read_output_files;

    SHALLOW_WATER_STATE(const STREAM_TYPE stream_type);
    ~SHALLOW_WATER_STATE();

    void Write_Output_Files(const int frame);
    void Read_Output_Files(const int frame);
//#####################################################################
};
}
#endif

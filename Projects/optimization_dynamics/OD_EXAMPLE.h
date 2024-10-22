//#####################################################################
// Copyright 2015
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __OD_EXAMPLE__
#define __OD_EXAMPLE__
#include <Core/Utilities/Find_Type.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
namespace PhysBAM{

template<class TV> class DEBUG_PARTICLES;
template<class TV> class DEFORMABLE_PARTICLES;
template<class TV> class STRUCTURE;
template<class T,int d> class CONSTITUTIVE_MODEL;

template<class TV>
class OD_EXAMPLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    STREAM_TYPE stream_type;
    DEBUG_PARTICLES<TV>& debug_particles;
    DEFORMABLE_PARTICLES<TV>& particles;
    ARRAY<CONSTITUTIVE_MODEL<T,TV::m>*> constitutive_models;

    ARRAY<STRUCTURE<TV>*> structures;

    int last_frame;
    std::string frame_title;
    int write_substeps_level;
    int substeps_delay_frame;
    VIEWER_DIR viewer_dir;
    int restart;
    T dt,time,frame_dt,min_dt,max_dt;
    bool print_stats;

    bool test_diff;
    int threads;

    OD_EXAMPLE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    OD_EXAMPLE(const OD_EXAMPLE&) = delete;
    virtual ~OD_EXAMPLE();
    void operator=(const OD_EXAMPLE&) = delete;
    
    virtual void Write_Output_Files();
    virtual void Read_Output_Files();

    int Add_Structure(STRUCTURE<TV>* structure)
    {return structures.Append(structure);}

    template<class T_STRUCTURE> T_STRUCTURE
    Find_Structure(const int index=0)
    {return Find_Type<T_STRUCTURE>(structures,index);}

    template<class T_STRUCTURE> const T_STRUCTURE
    Find_Structure(const int index=0) const
    {return Find_Type<T_STRUCTURE>(structures,index);}

    virtual void Initialize()=0;
    virtual void Begin_Frame(const int frame)=0;
    virtual void End_Frame(const int frame)=0;
    virtual void Begin_Time_Step(const T time)=0;
    virtual void End_Time_Step(const T time)=0;
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MOTION_SEQUENCE__
#define __MOTION_SEQUENCE__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>

namespace PhysBAM{

template<class T,class T2>
class MOTION_SEQUENCE
{
    typedef VECTOR<T,1> TV;
public:
    GRID<TV> time_grid;
    ARRAY<ARRAY<T2,VECTOR<int,1> > > trajectories;
    ARRAY<ARRAY<bool,VECTOR<int,1> > > valid;
    ARRAY<std::string> names;

    //INTERPOLATION_UNIFORM<GRID<TV>,T2>* interpolation;
    //LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T2> default_interpolation;
    HASHTABLE<std::string,int> name_to_track_index;

    MOTION_SEQUENCE()
        //:interpolation(&default_interpolation)
    {}

    MOTION_SEQUENCE(const MOTION_SEQUENCE<T,T2>& motion_input)
    {
        //if(motion_input.interpolation==&motion_input.default_interpolation) interpolation=&default_interpolation;
        //else interpolation=motion_input.interpolation;
        //interpolation=&default_interpolation;
        *this=motion_input;
    }

    MOTION_SEQUENCE<T,T2>& operator=(const MOTION_SEQUENCE<T,T2>& motion_input)
    {time_grid=motion_input.time_grid;
    trajectories=motion_input.trajectories;
    valid=motion_input.valid;
    names=motion_input.names;
    name_to_track_index=motion_input.name_to_track_index;
    //if(motion_input.interpolation==&motion_input.default_interpolation) interpolation=&default_interpolation;
    //else interpolation=motion_input.interpolation;
    return *this;}

    void Initialize(const int number_of_tracks,const GRID<TV>& time_grid_input)
    {time_grid=time_grid_input;trajectories.Resize(number_of_tracks);valid.Resize(number_of_tracks);names.Resize(number_of_tracks);
    for(int i=1;i<=number_of_tracks;i++){trajectories(i).Resize(1,time_grid.m);valid(i).Resize(1,time_grid.m);}}

    T2 X(const int marker,const T t) const
    {return Frame_X(marker,(int)((t-time_grid.domain.min_corner.x)/(time_grid.domain.max_corner.x-time_grid.domain.min_corner.x)*(time_grid.counts.x-1)+1));}
    //{if(time_grid.m==1) return trajectories(marker)(1);}
    //else return interpolation->Clamped_To_Array(time_grid,trajectories(marker),VECTOR<T,1>(t));}

    const T2& Frame_X(const int marker,const int frame) const
    {return trajectories(marker)(frame);}

    T2& Frame_X(const int marker,const int frame)
    {return trajectories(marker)(frame);}

    //void Set_Custom_Interpolation(const INTERPOLATION_UNIFORM<T,T2,GRID<VECTOR<T,1> > >* interpolation_input)
    //{interpolation=interpolation_input;}

    int Track_Index(const std::string& name)
    {int index;
    if(name_to_track_index.Get(name,index)) return index;
    else return -1;} 

    void Update_Name_Lookup()
    {name_to_track_index.Clean_Memory();for(int i=1;i<=names.m;i++) name_to_track_index.Insert(names(i),i);}

    void Set_Frame_Rate(const T frame_rate,const T initial_time=0)
    {time_grid.Initialize(time_grid.m,initial_time,initial_time+frame_rate*(time_grid.m-1));}

    template<class RW> 
    void Read(std::istream& input_stream)
    {Read_Binary<RW>(input_stream,time_grid,trajectories,valid,names);Update_Name_Lookup();}

    template<class RW> 
    void Write(std::ostream& output_stream) const
    {Write_Binary<RW>(output_stream,time_grid,trajectories,valid,names);}

//#####################################################################
};
}
#include <PhysBAM_Dynamics/Read_Write/Motion/READ_WRITE_MOTION_SEQUENCE.h>
#endif

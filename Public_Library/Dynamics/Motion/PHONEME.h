//#####################################################################
// Copyright 2005, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PHONEME__
#define __PHONEME__

#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_PDE/Interpolation/INTERPOLATION_UNIFORM.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_STATE.h>
namespace PhysBAM{

template<class T>
class PHONEME
{
    typedef VECTOR<T,3> TV;
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    std::string name,previous_name,next_name;
    int frame_length;
    T time_length;
    ARRAY<ARRAY<T> ,VECTOR<int,1> > controls;
    const INTERPOLATION_UNIFORM<VECTOR<T,1>,ARRAY<T> >* interpolation;

    PHONEME()
        :interpolation(0)
    {}

    ~PHONEME()
    {}

    void Set_Names(const std::string& name_input,const std::string& previous_name_input,const std::string& next_name_input)
    {name=name_input;previous_name=previous_name_input;next_name=next_name_input;}

    void Set_Frame_Length(const int frame_length_input)
    {frame_length=frame_length_input;}

    void Set_Time_Length(const T time_length_input)
    {time_length=time_length_input;}

    void Resize_Controls(const int leading_context_frames=0,const int trailing_context_frames=0)
    {controls.Resize(1-leading_context_frames,frame_length+trailing_context_frames);}

    template <class RW>
    void Interpolate_Rigid_Body_State(RIGID_BODY_STATE<VECTOR<T,3> >& rigid_body_state,const std::string rigid_transform_data_directory,const std::string rigid_transform_filename_prefix,const int initial_frame,const T time)
    {T fractional_frame=((T)(frame_length-1))*time/time_length+(T)initial_frame;
    int first_frame=(int)fractional_frame,second_frame=first_frame+1;T interpolation_fraction=fractional_frame-(T)first_frame;
    if(frame_length<=2){rigid_body_state=RIGID_BODY_STATE<VECTOR<T,3> >();return;} // HACK
    FRAME<TV> first_frame_rigid_transform,second_frame_rigid_transform;
    Read_From_File(rigid_transform_data_directory+"/"+rigid_transform_filename_prefix+"."+Number_To_String(first_frame),first_frame_rigid_transform);
    Read_From_File(rigid_transform_data_directory+"/"+rigid_transform_filename_prefix+"."+Number_To_String(second_frame),second_frame_rigid_transform);
    rigid_body_state.time=time;
    rigid_body_state.frame=FRAME<TV>::Interpolation(first_frame_rigid_transform,second_frame_rigid_transform,interpolation_fraction);
    rigid_body_state.velocity=(second_frame_rigid_transform.t-first_frame_rigid_transform.t)*((T)frame_length)/time_length;
    rigid_body_state.angular_velocity=(second_frame_rigid_transform.r*first_frame_rigid_transform.r.Inverse()).Rotation_Vector()*((T)frame_length)/time_length;}

    template <class RW>
    void Interpolate_Positions(ARRAY<VECTOR<T,3> >& positions,const std::string position_data_directory,const std::string position_filename_prefix,const int initial_frame,const T time)
    {T fractional_frame=((T)(frame_length-1))*time/time_length+(T)initial_frame;
    int first_frame=(int)fractional_frame,second_frame=first_frame+1;T interpolation_fraction=fractional_frame-(T)first_frame;
    if(frame_length<=2) first_frame=second_frame=0; // HACK
    ARRAY<VECTOR<T,3> > first_frame_positions,second_frame_positions;
    Read_From_File(position_data_directory+"/"+position_filename_prefix+"."+Number_To_String(first_frame),first_frame_positions);
    Read_From_File(position_data_directory+"/"+position_filename_prefix+"."+Number_To_String(second_frame),second_frame_positions);
    if(!positions.m) positions.Resize(first_frame_positions.m);
    ARRAY<VECTOR<T,3> >::copy((T)1-interpolation_fraction,first_frame_positions,interpolation_fraction,second_frame_positions,positions);}

    ARRAY<T> Controls(const T time) const
    {assert(interpolation);GRID<VECTOR<T,1> > time_grid(frame_length,0,time_length);return interpolation->Clamped_To_Array(time_grid,controls,time);}

    ARRAY<T> Controls(const T time,const GRID<VECTOR<T,1> >& time_grid) const
    {assert(interpolation);return interpolation->Clamped_To_Array(time_grid,controls,VECTOR<T,1>(time));}

    const ARRAY<T>& Frame_Controls(const int frame) const
    {return controls(frame);}

    ARRAY<T>& Frame_Controls(const int frame)
    {return controls(frame);}

    RANGE<VECTOR<T,1> > Time_Range() const
    {return RANGE<VECTOR<T,1> >(time_length*((T)(controls.domain.min_corner.x-1))/((T)(frame_length-1)),time_length*((T)(controls.domain.max_corner.x-1))/((T)(frame_length-1)));}

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,name,previous_name,next_name,frame_length,time_length,controls);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,name,previous_name,next_name,frame_length,time_length,controls);}

//#####################################################################
};
}
#endif

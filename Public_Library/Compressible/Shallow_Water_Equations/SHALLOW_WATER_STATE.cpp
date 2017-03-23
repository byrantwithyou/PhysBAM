//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Compressible/Shallow_Water_Equations/SHALLOW_WATER.h>
#include <Compressible/Shallow_Water_Equations/SHALLOW_WATER_STATE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SHALLOW_WATER_STATE<TV>::
SHALLOW_WATER_STATE(const STREAM_TYPE stream_type)
    :stream_type(stream_type),ghost(3),initial_time(0),
    last_frame(100),write_substeps_level(-1),substeps_delay_frame(-1),
    output_directory("output"),data_directory("../../Public_Data"),use_test_output(false),
    restart(0),dt(0),time(0),frame_dt((T)1/24),min_dt(0),max_dt(frame_dt),cfl(1),
    shallow_water(*new SHALLOW_WATER<TV>(grid,U))
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SHALLOW_WATER_STATE<TV>::
~SHALLOW_WATER_STATE()
{
    delete &shallow_water;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void SHALLOW_WATER_STATE<TV>::
Write_Output_Files(const int frame)
{
    std::string f=LOG::sprintf("%d",frame);
    if(this->use_test_output){
        std::string file=LOG::sprintf("%s/%s-%03d.txt",output_directory.c_str(),test_output_prefix.c_str(),frame);
        OCTAVE_OUTPUT<T> oo(file.c_str());}

    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",grid);
    FILE_UTILITIES::Write_To_File(stream_type,LOG::sprintf("%s/%d/restart_data",output_directory.c_str(),frame),time);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV> void SHALLOW_WATER_STATE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=LOG::sprintf("%d",frame);
    FILE_UTILITIES::Read_From_File(stream_type,LOG::sprintf("%s/%d/restart_data",output_directory.c_str(),frame),time);
}
//#####################################################################
namespace PhysBAM{
template class SHALLOW_WATER_STATE<VECTOR<float,1> >;
template class SHALLOW_WATER_STATE<VECTOR<double,1> >;
template class SHALLOW_WATER_STATE<VECTOR<float,2> >;
template class SHALLOW_WATER_STATE<VECTOR<double,2> >;
}

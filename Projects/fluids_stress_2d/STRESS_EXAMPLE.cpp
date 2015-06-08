//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include "STRESS_EXAMPLE.h"
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> STRESS_EXAMPLE<TV>::
STRESS_EXAMPLE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :stream_type(stream_type),initial_time(0),last_frame(100),write_substeps_level(-1),substeps_delay_frame(-1),
    write_output_files(true),output_directory("output"),restart(0),number_of_ghost_cells(5),dt(1),time(0),
    time_steps_per_frame(1),use_advection(true),use_reduced_advection(true),inv_Wi(0),use_du_terms(true),
    grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    levelset(grid,bc_phi,number_of_ghost_cells),debug_particles(*new DEBUG_PARTICLES<TV>)
{
    debug_particles.debug_particles.template Add_Array<TV>(ATTRIBUTE_ID_V);
    debug_particles.debug_particles.template Add_Array<T>(ATTRIBUTE_ID_DISPLAY_SIZE);
    parse_args.Parse();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> STRESS_EXAMPLE<TV>::
~STRESS_EXAMPLE()
{
    delete &debug_particles;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void STRESS_EXAMPLE<TV>::
Write_Output_Files(const int frame)
{
    if(!write_output_files) return;
    std::string f=LOG::sprintf("%d",frame);

    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",grid);
    FILE_UTILITIES::Write_To_File(stream_type,LOG::sprintf("%s/%d/levelset",output_directory.c_str(),frame),levelset);
    debug_particles.Write_Debug_Particles(stream_type,output_directory,frame);
    FILE_UTILITIES::Write_To_File(stream_type,LOG::sprintf("%s/%d/restart_data",output_directory.c_str(),frame),
        time,face_velocities,prev_face_velocities);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV> void STRESS_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=LOG::sprintf("%d",frame);
    FILE_UTILITIES::Read_From_File(stream_type,LOG::sprintf("%s/%d/levelset",output_directory.c_str(),frame),levelset);
    FILE_UTILITIES::Read_From_File(stream_type,LOG::sprintf("%s/%d/restart_data",output_directory.c_str(),frame),
        time,face_velocities,prev_face_velocities);
}
//#####################################################################
namespace PhysBAM{
template class STRESS_EXAMPLE<VECTOR<float,2> >;
template class STRESS_EXAMPLE<VECTOR<float,3> >;
template class STRESS_EXAMPLE<VECTOR<double,2> >;
template class STRESS_EXAMPLE<VECTOR<double,3> >;
}

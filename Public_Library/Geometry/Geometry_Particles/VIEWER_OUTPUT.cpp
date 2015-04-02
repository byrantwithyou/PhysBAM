//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Log/LOG.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> VIEWER_OUTPUT<TV>::
VIEWER_OUTPUT(STREAM_TYPE stream_type,const GRID<TV>& grid,const std::string& output_directory)
    :frame(0),output_directory(output_directory),stream_type(stream_type),grid(grid),debug_particles(*new DEBUG_PARTICLES<TV>)
{
    Singleton(this);
    debug_particles.debug_particles.template Add_Array<T>(ATTRIBUTE_ID_DISPLAY_SIZE);
    debug_particles.debug_particles.template Add_Array<TV>(ATTRIBUTE_ID_V);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> VIEWER_OUTPUT<TV>::
~VIEWER_OUTPUT()
{
    Singleton(0);
    delete &debug_particles;
}
//#####################################################################
// Function Singleton
//#####################################################################
template<class TV> VIEWER_OUTPUT<TV>* VIEWER_OUTPUT<TV>::
Singleton(VIEWER_OUTPUT* vo)
{
    static VIEWER_OUTPUT* viewer_output=0;
    VIEWER_OUTPUT* tmp=viewer_output;
    if(vo) viewer_output=vo;
    return tmp;
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void VIEWER_OUTPUT<TV>::
Flush_Frame(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const char* title)
{
    if(frame==0){
        FILE_UTILITIES::Create_Directory(output_directory);
        FILE_UTILITIES::Create_Directory(output_directory+"/common");
        LOG::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",false);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid.gz",grid);
        FILE_UTILITIES::Write_To_Text_File(output_directory+"/common/first_frame",frame,"\n");}

    debug_particles.Write_Debug_Particles(stream_type,output_directory,frame);

    std::string frame_directory=LOG::sprintf("%s/%i",output_directory.c_str(),frame);
    FILE_UTILITIES::Write_To_File(stream_type,frame_directory+"/mac_velocities.gz",face_velocities);
    if(title) FILE_UTILITIES::Write_To_Text_File(frame_directory+"/frame_title",title);

    FILE_UTILITIES::Write_To_Text_File(output_directory+"/common/last_frame",frame,"\n");
    frame++;
}
//#####################################################################
// Function Flush_Frame
//#####################################################################
template<class TV> void VIEWER_OUTPUT<TV>::
Flush_Frame(const char* title)
{
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities(grid);
    Flush_Frame(face_velocities,title);
}
namespace PhysBAM{
template class VIEWER_OUTPUT<VECTOR<float,1> >;
template class VIEWER_OUTPUT<VECTOR<float,2> >;
template class VIEWER_OUTPUT<VECTOR<float,3> >;
template class VIEWER_OUTPUT<VECTOR<double,1> >;
template class VIEWER_OUTPUT<VECTOR<double,2> >;
template class VIEWER_OUTPUT<VECTOR<double,3> >;
}

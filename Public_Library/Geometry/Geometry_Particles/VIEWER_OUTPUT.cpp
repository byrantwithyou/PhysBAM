//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
VIEWER_OUTPUT::
VIEWER_OUTPUT(STREAM_TYPE stream_type,VIEWER_DIR& viewer_dir)
    :viewer_dir(viewer_dir),stream_type(stream_type)
{
    Singleton(this);
}
//#####################################################################
// Destructor
//#####################################################################
VIEWER_OUTPUT::
~VIEWER_OUTPUT()
{
    Singleton(0);
}
//#####################################################################
// Function Singleton
//#####################################################################
VIEWER_OUTPUT* VIEWER_OUTPUT::
Singleton(VIEWER_OUTPUT* vo)
{
    static VIEWER_OUTPUT* viewer_output=0;
    VIEWER_OUTPUT* tmp=viewer_output;
    if(vo) viewer_output=vo;
    return tmp;
}
//#####################################################################
// Function Flush_Frame
//#####################################################################
void VIEWER_OUTPUT::
Flush_Frame(const char* title)
{
    viewer_dir.Start_Directory(0,title);
    if(viewer_dir.First_Frame())
        for(const auto& f:common_entries) f();
    for(const auto& f:entries) f();
    viewer_dir.Finish_Directory();
}
//#####################################################################
// Function Use_Debug_Particles
//#####################################################################
template<class TV> void
Use_Debug_Particles()
{
    VIEWER_OUTPUT::Singleton()->Use_Debug_Particles<TV>();
}
//#####################################################################
// Function Use_Debug_Particles
//#####################################################################
template<class TV> void VIEWER_OUTPUT::
Use_Debug_Particles()
{
    auto& p=Get_Debug_Particles<TV>();
    p.debug_particles.template Add_Array<typename TV::SCALAR>("display_size");
    p.debug_particles.template Add_Array<TV>("V");
    entries.Append([&p,this]()
        {
            p.Write_Debug_Particles(stream_type,viewer_dir);
        });
}
template void Use_Debug_Particles<VECTOR<float,1> >();
template void Use_Debug_Particles<VECTOR<float,2> >();
template void Use_Debug_Particles<VECTOR<float,3> >();
template void Use_Debug_Particles<VECTOR<double,1> >();
template void Use_Debug_Particles<VECTOR<double,2> >();
template void Use_Debug_Particles<VECTOR<double,3> >();
}

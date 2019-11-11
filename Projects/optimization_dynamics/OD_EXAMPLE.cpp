//#####################################################################
// Copyright 2015
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Read_Write/TYPED_STREAM.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Topology_Based_Geometry/STRUCTURE.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include "OD_EXAMPLE.h"
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> OD_EXAMPLE<TV>::
OD_EXAMPLE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :stream_type(stream_type_input),debug_particles(*new DEBUG_PARTICLES<TV>),particles(*new DEFORMABLE_PARTICLES<TV>),
    last_frame(100),
    write_substeps_level(-1),substeps_delay_frame(-1),viewer_dir("output"),
    restart(0),dt(0),time(0),frame_dt((T)1/24),min_dt(1e-8),max_dt(frame_dt),
    print_stats(false),test_diff(false),threads(1)
{
    parse_args.Parse();
    particles.Store_Velocity();
    particles.Store_Mass();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> OD_EXAMPLE<TV>::
~OD_EXAMPLE()
{
    structures.Delete_Pointers_And_Clean_Memory();
    delete &debug_particles;
    delete &particles;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void OD_EXAMPLE<TV>::
Write_Output_Files()
{
    Write_To_File(stream_type,viewer_dir.current_directory+"/particles",particles);
    Write_To_File(stream_type,viewer_dir.current_directory+"/restart_data",time);

    FILE_OSTREAM output;
    Safe_Open_Output(output,stream_type,viewer_dir.current_directory+"/structures");
    Write_Binary(output,structures.m);
    for(int k=0;k<structures.m;k++) structures(k)->Write_Structure(output);

    for(int i=0;i<particles.X.m;i++){
        Add_Debug_Particle(particles.X(i),VECTOR<T,3>(1,1,1));
        Debug_Particle_Set_Attribute<TV>("V",particles.V(i));}
    debug_particles.Write_Debug_Particles(stream_type,viewer_dir);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV> void OD_EXAMPLE<TV>::
Read_Output_Files()
{
    Read_From_File(viewer_dir.current_directory+"/particles",particles);
    Read_From_File(viewer_dir.current_directory+"/restart_data",time);

    FILE_ISTREAM input;
    Safe_Open_Input(input,viewer_dir.current_directory+"/structures");

    int m;Read_Binary(input,m);
    if(!structures.m){
        structures.Resize(m);
        for(int k=0;k<structures.m;k++) structures(k)=STRUCTURE<TV>::Create_Structure(input,particles);}
    else if(structures.m<=m){
        int old_number_of_structures=structures.m;structures.Resize(m);
        for(int k=0;k<old_number_of_structures;k++) structures(k)->Read_Structure(input);
        for(int k=old_number_of_structures;k<m;k++) structures(k)=STRUCTURE<TV>::Create_Structure(input,particles);}
    else{
        LOG::cout<<"Current number of structures ("<<structures.m<<") is greater than number in file ("<<m<<").";
        PHYSBAM_FATAL_ERROR();}
}

namespace PhysBAM{
template class OD_EXAMPLE<VECTOR<float,2> >;
template class OD_EXAMPLE<VECTOR<float,3> >;
template class OD_EXAMPLE<VECTOR<double,2> >;
template class OD_EXAMPLE<VECTOR<double,3> >;
}

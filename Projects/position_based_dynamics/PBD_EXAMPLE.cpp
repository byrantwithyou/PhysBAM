//#####################################################################
// Copyright 2015
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Read_Write/TYPED_STREAM.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include "PBD_EXAMPLE.h"
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PBD_EXAMPLE<TV>::
PBD_EXAMPLE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :stream_type(stream_type_input),debug_particles(*new DEBUG_PARTICLES<TV>),
    last_frame(100),
    write_substeps_level(-1),substeps_delay_frame(-1),viewer_dir("output"),
    restart(0),dt(0),time(0),frame_dt((T)1/24),min_dt(1e-8),max_dt(frame_dt),
    print_stats(false),solver_iterations(1),test_diff(false),threads(1)
{
    parse_args.Parse();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PBD_EXAMPLE<TV>::
~PBD_EXAMPLE()
{
    delete &debug_particles;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void PBD_EXAMPLE<TV>::
Write_Output_Files()
{
    Write_To_File(stream_type,viewer_dir.current_directory+"/inverse_mass",w);
    Write_To_File(stream_type,viewer_dir.current_directory+"/positions",X);
    Write_To_File(stream_type,viewer_dir.current_directory+"/velocities",V);
    Write_To_File(stream_type,viewer_dir.current_directory+"/restart_data",time);

    for(int i=0;i<X.m;i++){
        Add_Debug_Particle(X(i),VECTOR<T,3>(1,1,1));
        Debug_Particle_Set_Attribute<TV>("V",V(i));}
    debug_particles.Write_Debug_Particles(stream_type,viewer_dir);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV> void PBD_EXAMPLE<TV>::
Read_Output_Files()
{
    Read_From_File(viewer_dir.current_directory+"/inverse_mass",w);
    Read_From_File(viewer_dir.current_directory+"/positions",X);
    Read_From_File(viewer_dir.current_directory+"/velocities",V);
    Read_From_File(viewer_dir.current_directory+"/restart_data",time);
}

//#####################################################################
// Function Add_Constraints
//#####################################################################
template<class TV> int PBD_EXAMPLE<TV>::
Add_Constraints(PBD_CONSTRAINTS_BASE<TV>& c)
{
    return constraints.Append(&c);
}

namespace PhysBAM{
template class PBD_EXAMPLE<VECTOR<float,2> >;
template class PBD_EXAMPLE<VECTOR<float,3> >;
template class PBD_EXAMPLE<VECTOR<double,2> >;
template class PBD_EXAMPLE<VECTOR<double,3> >;
}

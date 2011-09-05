//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Dynamics/Particles/COMPRESSIBLE_FLUID_PARTICLES.h>
#include "Conversion_Tools/READ_AEROF_FLUID.h"

using namespace PhysBAM;
int main(int argc,char* argv[])
{
    typedef float T;
    typedef VECTOR<T,3> TV;
    typedef double RW;

    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-input_mesh","mesh");
    parse_args.Add_String_Argument("-result_dir","results");
    parse_args.Parse(argc,argv);
    std::string mesh_filename=parse_args.Get_String_Value("-input_mesh");
    std::string results_prefix=parse_args.Get_String_Value("-result_dir");

//##########################  INITIALIZATION  #########################
    TETRAHEDRON_MESH mesh;
    GEOMETRY_PARTICLES<TV> particles;
    TETRAHEDRALIZED_VOLUME<T> ale_fluid_mesh(mesh,particles);
    CONVERSION_TOOLS::Read_AeroF_Fluid_Mesh(mesh_filename,ale_fluid_mesh);

    COMPRESSIBLE_FLUID_PARTICLES<TV> compressible_fluid_data;
    compressible_fluid_data.array_collection->Resize(particles.array_collection->Size());
    CONVERSION_TOOLS::Read_AeroF_Fluid_Data(results_prefix,compressible_fluid_data);
//#####################################################################

    return 0;
}
//#####################################################################

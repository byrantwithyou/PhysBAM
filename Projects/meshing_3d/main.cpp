//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class main for 3D meshing w/ opengl visualization
//##################################################################### 
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Dynamics/Particles/DYNAMICS_PARTICLES_FORWARD.h>
#include "Generic/GENERIC_EXAMPLE.h"
#include "MESHING_DRIVER.h"
//#include "Sphere/SPHERE_EXAMPLE.h"
//#include "Two_Spheres/TWO_SPHERE_EXAMPLE.h"
//#######################################################################
// Function main
//#######################################################################
int main(int argc,char **argv)
{
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));

    std::string opt_out,opt_data;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-o",&opt_out,"dir","output directory");
    parse_args.Add("-d",&opt_data,"dir","data directory");
    parse_args.Set_Extra_Arguments(1,"<parameter file>");
    parse_args.Parse();

    std::string parameter_file=parse_args.Extra_Arg(0);
    GENERIC_EXAMPLE<T> example(stream_type,parameter_file);
    example.tetrahedral_meshing.output_directory=opt_out;
    example.data_directory=opt_data;

    MESHING_DRIVER<T> driver(example);
    driver.Execute_Main_Program();
    return 0;
}
//#######################################################################

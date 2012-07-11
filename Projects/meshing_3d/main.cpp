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

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_String_Argument("-o","","output directory","output directory");
    parse_args.Add_String_Argument("-d","","data directory","data directory");
    parse_args.Set_Extra_Arguments(1,"<parameter file>");
    parse_args.Parse();

    std::string parameter_file=parse_args.Extra_Arg(0);
    GENERIC_EXAMPLE<T> example(stream_type,parameter_file);

    if(parse_args.Is_Value_Set("-o")) example.tetrahedral_meshing.output_directory=parse_args.Get_String_Value("-o");
    if(parse_args.Is_Value_Set("-d")) example.data_directory=parse_args.Get_String_Value("-d");

    MESHING_DRIVER<T> driver(example);
    driver.Execute_Main_Program();
    return 0;
}
//#######################################################################

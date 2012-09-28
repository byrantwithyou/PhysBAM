//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include "Generic/GENERIC_RENDER_EXAMPLE.h"
#include "RAY_TRACING_DRIVER_WITH_PREVIEW.h"

using namespace PhysBAM;
int main(int argc, char *argv[]) 
{  
    std::string scene_filename;
    int frame_number=0;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Extra(&scene_filename,"scene file","scene file");
    parse_args.Extra(&frame_number,"frame number","frame number");
    parse_args.Parse();
    if(parse_args.unclaimed_arguments){parse_args.Print_Usage();exit(0);}

    GENERIC_RENDER_EXAMPLE<double,float> example(scene_filename,frame_number);
    RAY_TRACING_DRIVER_WITH_PREVIEW<double>(example).Execute_Main_Program();

    return 0;
}
//#####################################################################

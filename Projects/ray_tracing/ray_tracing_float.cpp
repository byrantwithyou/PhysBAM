//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Duc Nguyen, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include "Generic/GENERIC_RENDER_EXAMPLE.h"
#include "RAY_TRACING_DRIVER.h"
#include "RAY_TRACING_DRIVER_WITH_PREVIEW.h"
using namespace PhysBAM;
int main(int argc, char *argv[]) 
{  
    PROCESS_UTILITIES::Set_Backtrace(true);

    std::string scene_filename;
    int frame_number;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Extra(&scene_filename, "scene file", "scene file");
    parse_args.Extra(&frame_number, "frame number","frame number");
    parse_args.Parse();
    if(parse_args.unclaimed_arguments){parse_args.Print_Usage();exit(0);}

    GENERIC_RENDER_EXAMPLE<float,float> example(scene_filename,frame_number);
    RAY_TRACING_DRIVER<float>(example).Execute_Main_Program();

    return 0;
}
//#####################################################################

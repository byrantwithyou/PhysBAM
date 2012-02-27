//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include "LOCKING_TEST.h"

using namespace PhysBAM;

void Parse_Arguments(int argc,char *argv[]);

//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{
    LOCKING_TEST<double> test;
    test.Parse(argc,argv);

    // for(int i=1;i<argc;i++){
        // if(!strcmp(argv[i],"-affine_velicities")) affine_velicities=true;
        // else if (!strcmp(argv[i],"-bilinear_velocities")) affine_velicities=false;
        // else {LOG::cerr<<"Invalid argument"<<std::endl;exit(1);}}



    return 0;
}

//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include "LOCKING_TEST.h"

using namespace PhysBAM;

//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{
    LOCKING_TEST<double> test;
    test.Parse_Arguments(argc,argv);
    test.Compute_System_Matrix();
    test.Print_System_Matrix();

    return 0;
}

//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

using namespace PhysBAM;

typedef float RW;
typedef double T;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,TV::m> TV_INT;

int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Parse();

    return 0;
}


//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __STANDARD_TESTS_2D__
#define __STANDARD_TESTS_2D__

#include <Tools/Parsing/PARSE_ARGS.h>
#include "STOKES_MF_EXAMPLE.h"

namespace PhysBAM{

template<class TV> class STANDARD_TESTS;

template<class T>
class STANDARD_TESTS<VECTOR<T,2> >:public STOKES_MF_EXAMPLE<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef STOKES_MF_EXAMPLE<TV> BASE;
public:
    using BASE::test_number;using BASE::output_directory;
    using BASE::Add_Vertex;using BASE::Add_Edge;using BASE::Build_Grid;

    STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    void Initialize() override;

//#####################################################################
};
}

#endif

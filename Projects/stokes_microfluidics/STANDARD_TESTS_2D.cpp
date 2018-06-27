//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "STANDARD_TESTS_2D.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> STANDARD_TESTS<VECTOR<T,2> >::
STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :STOKES_MF_EXAMPLE<TV>(stream_type_input,parse_args)
{
    parse_args.Parse();
    output_directory=LOG::sprintf("Test_%i",test_number);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,2> >::
Initialize()
{
    switch(test_number)
    {
        case 1:{
            int v0=Add_Vertex(TV_INT(0,0));
            int v1=Add_Edge(v0,0,10);
            Add_Edge(v1,1,20);
            Add_Edge(v0,1,-10);
            Build_Grid(0.01);
        } break;
        default: PHYSBAM_FATAL_ERROR("test number not implemented");
    }
}
template class STANDARD_TESTS<VECTOR<float,2> >;
template class STANDARD_TESTS<VECTOR<double,2> >;
}

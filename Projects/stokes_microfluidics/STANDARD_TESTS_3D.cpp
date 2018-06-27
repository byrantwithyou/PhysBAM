//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "STANDARD_TESTS_3D.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> STANDARD_TESTS<VECTOR<T,3> >::
STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :STOKES_MF_EXAMPLE<TV>(stream_type_input,parse_args)
{
    parse_args.Parse();
    output_directory=LOG::sprintf("Test_3d_%i",test_number);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
Initialize()
{
    switch(test_number)
    {
        case 1:{
        } break;
        default: PHYSBAM_FATAL_ERROR("test number not implemented");
    }
}
template class STANDARD_TESTS<VECTOR<float,3> >;
template class STANDARD_TESTS<VECTOR<double,3> >;
}

//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "FLUIDS_COLOR_3D.h"

namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> FLUIDS_COLOR<VECTOR<T,3> >::
FLUIDS_COLOR(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :FLUIDS_COLOR_BASE<TV>(stream_type_input,parse_args)
{
    parse_args.Parse();

    if(!Initialize_Common_Example())
        Initialize_Example();

    After_Initialize_Example();
    parse_args.Parse();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> FLUIDS_COLOR<VECTOR<T,3> >::
~FLUIDS_COLOR()
{
}
//#####################################################################
// Function Initialize_Example
//#####################################################################
template<class T> void FLUIDS_COLOR<VECTOR<T,3> >::
Initialize_Example()
{
    switch(test_number){
        default: PHYSBAM_FATAL_ERROR("Missing test number");}
}
//#####################################################################
template class FLUIDS_COLOR<VECTOR<float,3> >;
template class FLUIDS_COLOR<VECTOR<double,3> >;
}

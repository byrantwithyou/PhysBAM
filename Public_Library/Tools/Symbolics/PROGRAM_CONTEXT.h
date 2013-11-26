//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYREGHT.txt.
//#####################################################################
// Class PROGRAM_CONTEXT
//#####################################################################
#ifndef __PROGRAM_CONTEXT__
#define __PROGRAM_CONTEXT__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Arrays/ARRAY_VIEW.h>
namespace PhysBAM{
template<class T> struct PROGRAM;
template<class T>
struct PROGRAM_CONTEXT
{
    ARRAY<T> reg;
    ARRAY_VIEW<T> data_in,data_out;

    PROGRAM_CONTEXT(const PROGRAM<T>& prog);
};
}
#endif
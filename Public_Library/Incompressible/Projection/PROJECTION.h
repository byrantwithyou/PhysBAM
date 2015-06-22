//#####################################################################
// Copyright 2004-2005, Ronald Fedkiw, Frank Losasso, Andy Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROJECTION
//#####################################################################
#ifndef __PROJECTION__
#define __PROJECTION__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Data_Structures/TRIPLE.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{

template<class T>
class PROJECTION:public NONCOPYABLE
{
public:
    T density;

protected: 
    bool use_non_zero_divergence;

public:
    PROJECTION()
        :use_non_zero_divergence(false)
    {
        Set_Density(1000);
    }

    virtual ~PROJECTION()
    {}

    void Set_Density(const T density_input=1000)
    {assert(density_input);density=density_input;}
//#####################################################################
};
}
#endif


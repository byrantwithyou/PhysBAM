//#####################################################################
// Copyright 2004-2005, Ronald Fedkiw, Frank Losasso, Andy Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROJECTION
//#####################################################################
#ifndef __PROJECTION__
#define __PROJECTION__

#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Data_Structures/TRIPLE.h>
namespace PhysBAM{

template<class T>
class PROJECTION
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

    PROJECTION(const PROJECTION&) = delete;
    void operator=(const PROJECTION&) = delete;
    
    virtual ~PROJECTION()
    {}

    void Set_Density(const T density_input=1000)
    {assert(density_input);density=density_input;}
//#####################################################################
};
}
#endif


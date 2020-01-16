 //#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUID_BOUNDARY_VECTOR_PB__
#define __FLUID_BOUNDARY_VECTOR_PB__
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include "FLUID_BOUNDARY_VECTOR.h"
namespace PhysBAM{

template<class TV>
struct FLUID_BOUNDARY_VECTOR_PB:public FLUID_BOUNDARY_VECTOR<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    HASHTABLE<FACE_INDEX<TV::m>,T> V;
 
    FLUID_BOUNDARY_VECTOR_PB()=default;
    virtual ~FLUID_BOUNDARY_VECTOR_PB()=default;
};
}
#endif

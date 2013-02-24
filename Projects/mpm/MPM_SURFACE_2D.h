//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_SURFACE_2D
//#####################################################################
#ifndef __MPM_SURFACE_2D__
#define __MPM_SURFACE_2D__

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include "MPM_SIMULATION.h"

namespace PhysBAM{

template<class TV>
class MPM_SURFACE_2D
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    const MPM_SIMULATION<TV>& sim;

    SEGMENTED_CURVE_2D<T> curve;
    ARRAY<TV> Xm;

    MPM_SURFACE_2D(const MPM_SIMULATION<TV>& sim);
    ~MPM_SURFACE_2D();

    void Initialize_With_A_Box(const T h,const RANGE<TV>& box);
    // void Build_Nearby_Phyxels_In_Material_Space(cosnt T h);
//#####################################################################
};
}
#endif

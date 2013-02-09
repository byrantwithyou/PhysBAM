//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_SIMULATION
//#####################################################################
#ifndef __MPM_SIMULATION__
#define __MPM_SIMULATION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX.h>
#include "MPM_PARTICLES.h"
#include "MPM_CONSTITUTIVE_MODEL.h"
#include "MPM_CUBIC_B_SPLINE.h"

namespace PhysBAM{

template<class TV>
class MPM_SIMULATION
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    MPM_PARTICLES<TV> particles;
    MPM_CONSTITUTIVE_MODEL<TV> constitutive_model;
    GRID<TV> grid;
    MPM_CUBIC_B_SPLINE<TV,3> grid_basis_function; 

    MPM_SIMULATION();
    ~MPM_SIMULATION();

    void Initialize();
//#####################################################################
};
}
#endif

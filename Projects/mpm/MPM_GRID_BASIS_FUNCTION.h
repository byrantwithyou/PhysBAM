//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_GRID_BASIS_FUNCTION
//#####################################################################
#ifndef __MPM_GRID_BASIS_FUNCTION__
#define __MPM_GRID_BASIS_FUNCTION__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>

namespace PhysBAM{

template<class TV>
class MPM_GRID_BASIS_FUNCTION
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    MPM_GRID_BASIS_FUNCTION();
    ~MPM_GRID_BASIS_FUNCTION();
//#####################################################################
};
}
#endif

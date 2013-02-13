//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_PARTICLES
//#####################################################################
#ifndef __MPM_PARTICLES__
#define __MPM_PARTICLES__

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <PhysBAM_Tools/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLES_PARTICLES_FORWARD.h>
#include "MPM_PARTICLES_FORWARD.h"

namespace PhysBAM{

template<class TV>
class MPM_PARTICLES:public CLONEABLE<MPM_PARTICLES<TV>,DEFORMABLE_PARTICLES<TV> > // X, V, mass
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<MPM_PARTICLES<TV>,DEFORMABLE_PARTICLES<TV> > BASE;
public:
    using BASE::X;using BASE::V;using BASE::mass;using BASE::Store_Velocity;using BASE::Store_Mass;using BASE::Remove_Array;using BASE::Get_Attribute_Index;using BASE::Remove_Array_Using_Index;

    ARRAY_VIEW<T> density;
    ARRAY_VIEW<T> volume;
    ARRAY_VIEW<MATRIX<T,TV::m> > Fe;
    ARRAY_VIEW<MATRIX<T,TV::m> > Fp;

    MPM_PARTICLES();
    virtual ~MPM_PARTICLES();

    void Initialize_X_As_A_Grid(const VECTOR<int,TV::m>& count,const RANGE<TV>& box);
//#####################################################################
};
}
#endif

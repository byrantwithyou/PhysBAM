//#####################################################################
// Copyright 2015, Peter Chen
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SMOKE_PARTICLES__
#define __SMOKE_PARTICLES__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Particles/PARTICLES.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>

namespace PhysBAM{

const ATTRIBUTE_ID ATTRIBUTE_ID_SMOKE_C(47);

template<class TV>
class SMOKE_PARTICLES:public CLONEABLE<SMOKE_PARTICLES<TV>,DEFORMABLE_PARTICLES<TV> >
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<SMOKE_PARTICLES<TV>,DEFORMABLE_PARTICLES<TV> > BASE;
public:
    using BASE::Add_Array;using BASE::Remove_Array;

    ARRAY_VIEW<MATRIX<T,TV::m> > C;

    SMOKE_PARTICLES();
    virtual ~SMOKE_PARTICLES();
//#####################################################################
};

}
#endif

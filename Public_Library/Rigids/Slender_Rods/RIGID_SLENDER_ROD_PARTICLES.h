//#####################################################################
// Copyright 2012, Steve Cook, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_SLENDER_ROD_PARTICLES
//#####################################################################
#ifndef __RIGID_SLENDER_ROD_PARTICLES__
#define __RIGID_SLENDER_ROD_PARTICLES__

#include <Core/Vectors/TWIST.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Clone/CLONEABLE.h>
#include <Tools/Particles/PARTICLES.h>
namespace PhysBAM{

template<class TV>
class RIGID_SLENDER_ROD_PARTICLES:public CLONEABLE<RIGID_SLENDER_ROD_PARTICLES<TV>,PARTICLES<TV> >
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<RIGID_SLENDER_ROD_PARTICLES<TV>,PARTICLES<TV> > BASE;
public:
    using BASE::Remove_Array;using BASE::Add_Array;

    ARRAY_VIEW<T> mass;
    ARRAY_VIEW<TWIST<TV> > twist;
    ARRAY_VIEW<TV> X;
    ARRAY_VIEW<TV> orientation;
    ARRAY_VIEW<T> length;

    RIGID_SLENDER_ROD_PARTICLES();
    virtual ~RIGID_SLENDER_ROD_PARTICLES();

    typename TV::SPIN Angular_Momentum(const int index) const
    {return twist(index).angular*Inertia_Scalar(index);}

    T Inertia_Scalar(const int index) const
    {return mass(index)*sqr(length(index))/12;}

    T One_Over_Inertia_Scalar(const int index) const
    {return 12/(mass(index)*sqr(length(index)));}

    TV Center(const int index) const
    {return X(index);}

//#####################################################################
};
}
#endif

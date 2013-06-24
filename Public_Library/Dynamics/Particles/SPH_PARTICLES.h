//#####################################################################
// Copyright 2004-2009, Michael Lentine, Frank Losasso, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPH_PARTICLES
//#####################################################################
#ifndef __SPH_PARTICLES__
#define __SPH_PARTICLES__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Clone/CLONEABLE.h>
#include <Tools/Particles/PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
namespace PhysBAM{

template<class TV>
class SPH_PARTICLES:public CLONEABLE<SPH_PARTICLES<TV>,GEOMETRY_PARTICLES<TV> >
{
    typedef CLONEABLE<SPH_PARTICLES<TV>,GEOMETRY_PARTICLES<TV> > BASE;
    typedef typename TV::SCALAR T;
public:
    using BASE::X;using BASE::V;

    SPH_PARTICLES();
    virtual ~SPH_PARTICLES();

    template<class T_INDICES>
    void Euler_Step(const T_INDICES& indices,const ARRAY<TV>& X,const T dt)
    {X.Subset(indices)+=dt*V.Subset(indices);}

    template<class T_INDICES>
    void Euler_Step(const T_INDICES& indices,const ARRAY<TV>& F,const ARRAY<T>& mass,const T dt)
    {V.Subset(indices)+=dt/mass.Subset(indices)*F.Subset(indices);}
};
}
#endif

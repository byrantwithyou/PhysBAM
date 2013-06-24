//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET_REMOVED_PARTICLES
//#####################################################################
#ifndef __PARTICLE_LEVELSET_REMOVED_PARTICLES__
#define __PARTICLE_LEVELSET_REMOVED_PARTICLES__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Clone/CLONEABLE.h>
#include <Dynamics/Particles/PARTICLE_LEVELSET_PARTICLES.h>
namespace PhysBAM{

template<class TV>
class PARTICLE_LEVELSET_REMOVED_PARTICLES:public CLONEABLE<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>,PARTICLE_LEVELSET_PARTICLES<TV> >
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>,PARTICLE_LEVELSET_PARTICLES<TV> > BASE;
public:
    using BASE::X;using BASE::V;

    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* next; //This should always be 0

    PARTICLE_LEVELSET_REMOVED_PARTICLES();
    ~PARTICLE_LEVELSET_REMOVED_PARTICLES();

    template<class T_INDICES>
    void Euler_Step(const T_INDICES& indices,const ARRAY<TV>& F,const ARRAY<T>& mass,const T dt)
    {V.Subset(indices)+=dt/mass.Subset(indices)*F.Subset(indices);}

    template<class T_INDICES>
    void Euler_Step(const T_INDICES& indices,const T dt)
    {X.Subset(indices)+=dt*V.Subset(indices);}

    void Euler_Step_Position(const T dt)
    {X+=dt*V;}

//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2007-2009, Jon Gretarsson, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMPRESSIBLE_FLUID_PARTICLES
//#####################################################################
#ifndef __COMPRESSIBLE_FLUID_PARTICLES__
#define __COMPRESSIBLE_FLUID_PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Particles/PARTICLES.h>
namespace PhysBAM{

template<class TV>
class COMPRESSIBLE_FLUID_PARTICLES:public CLONEABLE<COMPRESSIBLE_FLUID_PARTICLES<TV>,PARTICLES<TV> >
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<COMPRESSIBLE_FLUID_PARTICLES<TV>,PARTICLES<TV> > BASE;
public:

    ARRAY_VIEW<TV> X;
    ARRAY_VIEW<T> rho;
    ARRAY_VIEW<T> E;
    ARRAY_VIEW<T> phi;
    ARRAY_VIEW<TV> grad_phi;
    ARRAY_VIEW<TV> V;

    COMPRESSIBLE_FLUID_PARTICLES();
    virtual ~COMPRESSIBLE_FLUID_PARTICLES();

    template<class T_INDICES>
    void Euler_Step_Grad_Phi(const T_INDICES& indices,const ARRAY<TV>& X,const T dt)
    {X.Subset(indices)+=dt*grad_phi.Subset(indices);}

    template<class T_INDICES>
    void Euler_Step_Grad_Phi(const T_INDICES& indices,const ARRAY<TV>& F,const ARRAY<T>& mass,const T dt)
    {grad_phi.array.Subset(indices)+=dt/mass.Subset(indices)*F.Subset(indices);}

    template<class T_INDICES>
    void Euler_Step_Velocity(const T_INDICES& indices,const ARRAY<TV>& X,const T dt)
    {X.Subset(indices)+=dt*V.Subset(indices);}

    template<class T_INDICES>
    void Euler_Step_Velocity(const T_INDICES& indices,const ARRAY<TV>& F,const ARRAY<T>& mass,const T dt)
    {V.Subset(indices)+=dt/mass.Subset(indices)*F.Subset(indices);}
};
}
#endif

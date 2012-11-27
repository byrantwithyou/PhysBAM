//#####################################################################
// Copyright 2009, Nipun Kwatra, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GEOMETRY_PARTICLES
//#####################################################################
#ifndef __GEOMETRY_PARTICLES__
#define __GEOMETRY_PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Particles/PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
namespace PhysBAM{

template<class TV>
class GEOMETRY_PARTICLES:public CLONEABLE<GEOMETRY_PARTICLES<TV>,PARTICLES<TV> >
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<GEOMETRY_PARTICLES<TV>,PARTICLES<TV> > BASE;
public:
    using BASE::Remove_Array;using BASE::Add_Array;

    ARRAY_VIEW<TV> X,V;
    bool store_velocity;

    GEOMETRY_PARTICLES();
    virtual ~GEOMETRY_PARTICLES();

    void Store_Velocity(bool store=true)
    {store_velocity=store;if(store) Add_Array(ATTRIBUTE_ID_V,&V);else Remove_Array(ATTRIBUTE_ID_V);}

    template<class T_INDICES>
    void Euler_Step(const T_INDICES& indices,const T dt)
    {X.Subset(indices)+=dt*V.Subset(indices);}

    template<class T_INDICES>
    void Euler_Step(const T_INDICES& indices,const ARRAY<TV>& F,const ARRAY<T>& mass,const T dt)
    {V.Subset(indices)+=dt/mass.Subset(indices)*F.Subset(indices);}

    void Clone_Helper(const GEOMETRY_PARTICLES& particles)
    {Store_Velocity(particles.store_velocity);BASE::Clone_Helper(particles);}

//#####################################################################
};
}
#endif

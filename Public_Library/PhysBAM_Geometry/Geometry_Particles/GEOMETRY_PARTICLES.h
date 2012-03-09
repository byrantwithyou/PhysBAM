//#####################################################################
// Copyright 2009, Nipun Kwatra, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GEOMETRY_PARTICLES
//#####################################################################
#ifndef __GEOMETRY_PARTICLES__
#define __GEOMETRY_PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
namespace PhysBAM{

template<class TV>
class GEOMETRY_PARTICLES:public CLONEABLE<GEOMETRY_PARTICLES<TV>,POINT_CLOUD<TV> >
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<GEOMETRY_PARTICLES<TV>,POINT_CLOUD<TV> > BASE;
public:
    using BASE::array_collection;

    ARRAY_VIEW<TV> X,V;
    bool store_velocity;

    GEOMETRY_PARTICLES(ARRAY_COLLECTION* array_collection_input);
    GEOMETRY_PARTICLES();
    virtual ~GEOMETRY_PARTICLES();

    void Store_Velocity(bool store=true)
    {store_velocity=store;if(store) array_collection->Add_Array(ATTRIBUTE_ID_V,&V);else array_collection->Remove_Array(ATTRIBUTE_ID_V);}

    template<class T_INDICES>
    void Euler_Step(const T_INDICES& indices,const T dt)
    {X.Subset(indices)+=dt*V.Subset(indices);}

    template<class T_INDICES>
    void Euler_Step(const T_INDICES& indices,const ARRAY<TV>& F,const ARRAY<T>& mass,const T dt)
    {V.Subset(indices)+=dt/mass.Subset(indices)*F.Subset(indices);}

//#####################################################################
    void Initialize_Array_Collection();
//#####################################################################
};
}
#endif

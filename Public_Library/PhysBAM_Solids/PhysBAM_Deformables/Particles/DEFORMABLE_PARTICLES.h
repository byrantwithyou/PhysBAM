//#####################################################################
// Copyright 2004-2009, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_PARTICLES
//#####################################################################
#ifndef __DEFORMABLE_PARTICLES__
#define __DEFORMABLE_PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Particles/PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLES_PARTICLES_FORWARD.h>

namespace PhysBAM{

template<class TV> class SOFT_BINDINGS; 

template<class TV>
class DEFORMABLE_PARTICLES:public CLONEABLE<DEFORMABLE_PARTICLES<TV>,GEOMETRY_PARTICLES<TV> > // X, V
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<DEFORMABLE_PARTICLES<TV>,GEOMETRY_PARTICLES<TV> > BASE;
public:
    using BASE::X;using BASE::V;using BASE::Remove_Array;using BASE::Get_Attribute_Index;using BASE::Remove_Array_Using_Index;

    ARRAY_VIEW<T> mass;
    ARRAY_VIEW<T> one_over_mass;
    ARRAY_VIEW<T> effective_mass;
    ARRAY_VIEW<T> one_over_effective_mass;
    bool store_mass;

    DEFORMABLE_PARTICLES();
    virtual ~DEFORMABLE_PARTICLES();

    void Store_Mass(bool store=true)
    {store_mass=store;if(store) this->Add_Array(ATTRIBUTE_ID_MASS,&mass);else Remove_Array(ATTRIBUTE_ID_MASS);}

    T Min_Mass() const 
    {return mass.Size()?mass.Min():FLT_MAX;}

    T Max_Mass() const 
    {return mass.Size()?mass.Max():0;}

    void Euler_Step_Position(const T dt)
    {X+=dt*V;}

    template<class T_INDICES>
    void Euler_Step_Position(const T_INDICES& indices,const T dt)
    {X.Subset(indices)+=dt*V.Subset(indices);}

//#####################################################################
    TV Center_Of_Mass() const;
    void Compute_Auxiliary_Attributes(const SOFT_BINDINGS<TV>& soft_bindings);
    template<class T_INDICES> void Compute_Auxiliary_Attributes(const SOFT_BINDINGS<TV>& soft_bindings,const T_INDICES& indices,const bool copy_existing_elements=true);
//#####################################################################
};
}
#endif

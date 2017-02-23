//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_GRAVITY
//#####################################################################
#ifndef __RIGID_GRAVITY__
#define __RIGID_GRAVITY__

#include <Rigids/Forces_And_Torques/RIGID_POINTWISE_FORCE.h>
namespace PhysBAM{

template<class TV>
class RIGID_GRAVITY:public RIGID_POINTWISE_FORCE<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef RIGID_POINTWISE_FORCE<TV> BASE;
    typedef typename FORCE_ELEMENTS::ITERATOR ELEMENT_ITERATOR;
    TV gravity;
    
    using BASE::force_rigid_body_particles;using BASE::rigid_body_collection;
public:

    RIGID_GRAVITY(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,ARRAY<int>* influenced_rigid_body_particles_input=0,
        const TV& gravity_input=-(T)9.8*TV::Axis_Vector(1-(TV::m==1)))
        :BASE(rigid_body_collection_input,influenced_rigid_body_particles_input),gravity(gravity_input)
    {}

    virtual ~RIGID_GRAVITY()
    {}

    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const override
    {}

    int Velocity_Dependent_Forces_Size() const override
    {return 0;}

    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<T> aggregate,const T time) const override
    {}

    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const override
    {}

//#####################################################################
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const override;
    T Potential_Energy(const T time) const override;
//#####################################################################
};
}
#endif

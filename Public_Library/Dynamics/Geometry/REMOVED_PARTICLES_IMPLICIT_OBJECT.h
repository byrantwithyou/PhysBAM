//#####################################################################
// Copyright 2006, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REMOVED_PARTICLES_IMPLICIT_OBJECT
//#####################################################################
#ifndef __REMOVED_PARTICLES_IMPLICIT_OBJECT__
#define __REMOVED_PARTICLES_IMPLICIT_OBJECT__

#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Dynamics/Level_Sets/REMOVED_PARTICLES_PROCESSING.h>
namespace PhysBAM{

template<class TV>
class REMOVED_PARTICLES_IMPLICIT_OBJECT:public IMPLICIT_OBJECT<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef IMPLICIT_OBJECT<TV> BASE;
    using BASE::box;

    REMOVED_PARTICLES_PROCESSING<T>* particle_processing;
    LEVELSET_IMPLICIT_OBJECT<TV>* levelset;

    T blending_radius;

    REMOVED_PARTICLES_IMPLICIT_OBJECT(REMOVED_PARTICLES_PROCESSING<T>* particle_processing_input,LEVELSET_IMPLICIT_OBJECT<TV>* levelset_input)
        :particle_processing(particle_processing_input),levelset(levelset_input),blending_radius((T)0.05)
    {
        Update_Box();
    }

    ~REMOVED_PARTICLES_IMPLICIT_OBJECT()
    {delete particle_processing;delete levelset;}
    
    void Update_Box() override
    {box=particle_processing->particle_domain;box.Enlarge_To_Include_Box(levelset->box);}
    
    T operator()(const TV& location) const override
    {return particle_processing->Phi(location);}
    
    TV Normal(const TV& location,const int aggregate=-1) const override
    {assert((aggregate >= 1 && aggregate <= 6) || aggregate == -1);
    if(aggregate != -1)return box.Normal(aggregate);else return particle_processing->Normal(location);}
    
    T Integration_Step(const T phi) const override
    {return max(phi,particle_processing->tolerance);}
    
    virtual void Read(TYPED_ISTREAM input) override {PHYSBAM_FATAL_ERROR();}
    virtual void Write(TYPED_OSTREAM output) const override {PHYSBAM_FATAL_ERROR();}

//###########################################################################
};
}
#endif

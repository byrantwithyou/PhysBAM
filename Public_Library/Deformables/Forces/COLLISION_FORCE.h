//#####################################################################
// Copyright 2002-2008, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Igor Neverov, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COLLISION_FORCE
//#####################################################################
#ifndef __COLLISION_FORCE__
#define __COLLISION_FORCE__

#include <Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV>
class COLLISION_FORCE:public DEFORMABLES_FORCES<TV>
{
    typedef typename TV::SCALAR T;
public:
    T coefficient_of_friction;

    COLLISION_FORCE(DEFORMABLE_PARTICLES<TV>& particles_input)
        :DEFORMABLES_FORCES<TV>(particles_input),coefficient_of_friction((T)0.3)
    {}

    virtual ~COLLISION_FORCE() {}

    virtual void Apply_Friction(ARRAY_VIEW<TV> V,const T time) const=0;
};
}
#endif

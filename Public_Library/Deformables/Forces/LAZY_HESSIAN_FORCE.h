//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAZY_HESSIAN_FORCE
//#####################################################################
#ifndef __LAZY_HESSIAN_FORCE__
#define __LAZY_HESSIAN_FORCE__

#include <Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV>
class LAZY_HESSIAN_FORCE:public DEFORMABLES_FORCES<TV>
{
private:
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:

    LAZY_HESSIAN_FORCE(DEFORMABLE_PARTICLES<TV>& particles_input)
        :DEFORMABLES_FORCES<TV>(particles_input)
    {}

    virtual ~LAZY_HESSIAN_FORCE()
    {}

    virtual void Need_To_Recompute_Hessian(bool)=0;
//#####################################################################
};
}
#endif

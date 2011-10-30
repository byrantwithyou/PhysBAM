//#####################################################################
// Copyright 2011, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUINTIC_HEAVISIDE_TRANSITION
//#####################################################################
#ifndef __QUINTIC_HEAVISIDE_TRANSITION__
#define __QUINTIC_HEAVISIDE_TRANSITION__

#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/HEAVISIDE_TRANSITION.h>

namespace PhysBAM{

template<class T>
class QUINTIC_HEAVISIDE_TRANSITION: public HEAVISIDE_TRANSITION<T> 
{

private:

    inline T H_base (const T r) const
    {
        return (6*sqr(r)-15*r+10)*r*r*r;
    }

    inline T Hr_base (const T r) const
    {
        return 30*sqr(r)*sqr(1-r);
    }

    inline T Hrr_base (const T r) const
    {
        return 60*r*(2*sqr(r)-3*r+1);
    }
};
}

#endif

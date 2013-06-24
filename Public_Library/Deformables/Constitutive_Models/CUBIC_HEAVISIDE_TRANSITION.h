//#####################################################################
// Copyright 2011, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUBIC_HEAVISIDE_TRANSITION
//#####################################################################
#ifndef __CUBIC_HEAVISIDE_TRANSITION__
#define __CUBIC_HEAVISIDE_TRANSITION__

#include <Deformables/Constitutive_Models/HEAVISIDE_TRANSITION.h>

namespace PhysBAM{

template<class T>
class CUBIC_HEAVISIDE_TRANSITION: public HEAVISIDE_TRANSITION<T> 
{

private:

    inline T H_base (const T r) const
    {
        return sqr(r)*(-2*r+3);
    }

    inline T Hr_base (const T r) const
    {
        return 6*r*(1-r);
    }

    inline T Hrr_base (const T r) const
    {
        return -12*r+6;
    }
};
}

#endif

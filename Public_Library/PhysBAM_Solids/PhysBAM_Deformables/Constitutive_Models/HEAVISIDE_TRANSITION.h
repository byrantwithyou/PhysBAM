//#####################################################################
// Copyright 2011, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HEAVISIDE_TRANSITION
//#####################################################################
#ifndef __HEAVISIDE_TRANSITION__
#define __HEAVISIDE_TRANSITION__

namespace PhysBAM{

template<class T>
class HEAVISIDE_TRANSITION
{

private:

    T J_min,J_max,J_diff;

public:

    HEAVISIDE_TRANSITION () {}
    virtual ~HEAVISIDE_TRANSITION () {}

    inline void Initialize (const T input_J_min, const T input_J_max)
    {
        J_min  = input_J_min;
        J_max  = input_J_max;
        J_diff = J_max-J_min; 
    }
    
    inline T H (const T J) const
    {
        assert(J>J_min && J<J_max);
        T r = (J-J_min)/J_diff;
        return r*r*(-2*r+3);
    }

    inline T HJ (const T J) const
    {
        assert(J>J_min && J<J_max);
        T r = (J-J_min)/J_diff;
        return 6*r*(1-r)/J_diff;
    }
};
}
#endif

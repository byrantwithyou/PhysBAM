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

    T x_min,x_max,x_dif;
    
    virtual T H_base   (const T r) const = 0;
    virtual T Hr_base  (const T r) const = 0;
    virtual T Hrr_base (const T r) const = 0;

public:

    inline HEAVISIDE_TRANSITION ():
        x_min(0),x_max(1),x_dif(1){}

    inline HEAVISIDE_TRANSITION (const T input_x_min, const T input_x_max):
        x_min(input_x_min),x_max(input_x_max),x_dif(input_x_max-input_x_min)
    {
        assert(input_x_max>input_x_min);
    }

    inline ~HEAVISIDE_TRANSITION () {}
    
    inline void Initialize (const T input_x_min, const T input_x_max)
    {
        assert(input_x_max>input_x_min);
        x_min  = input_x_min;
        x_max  = input_x_max;
        x_dif = x_max-x_min; 
    }

    inline T H (const T x) const
    {
        assert(x>=x_min && x<=x_max);
        T r = (x-x_min)/x_dif;
        return this->H_base(r);
    }

    inline T Hx (const T x) const
    {
        assert(x>=x_min && x<=x_max);
        T r = (x-x_min)/x_dif;
        return this->Hr_base(r)/x_dif;
    }

    inline T Hxx (const T x) const
    {
        assert(x>=x_min && x<=x_max);
        T r = (x-x_min)/x_dif;
        return this->Hrr_base(r)/sqr(x_dif);
    }
};
}

#endif

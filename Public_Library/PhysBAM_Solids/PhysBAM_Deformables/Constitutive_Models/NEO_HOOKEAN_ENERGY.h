//#####################################################################
// Copyright 2011, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NEO_HOOKEAN_ENERGY
//#####################################################################
#ifndef __NEO_HOOKEAN_ENERGY__
#define __NEO_HOOKEAN_ENERGY__

namespace PhysBAM{

class NEO_HOOKEAN_ENERGY<T>
{

private:

    T mu;
    T lambda;

public:

    NEO_HOOKEAN_ENERGY () {}
    virtual ~NEO_HOOKEAN_ENERGY () {}

    inline void Initialize (const T input_mu, const T input_lambda)
    {
        mu = input_mu;
        lambda = input_lambda;
    }
    
    inline T E (const T x, const T y) const
    {
        assert(x>0 && y>0);
        return mu*(0.5*x*x + 0.5*y*y - 1 - log(xy)) + 0.5*lambda*sqr(log(xy)); 
    }

    inline T Ex (const T x, const T y) const
    {
        assert(x>0 && y>0);
        return mu*(x - 1/x) + lambda*log(xy)/x;
    }

    inline T Ey (const T x, const T y) const
    {
        assert(x>0 && y>0);
        return mu*(y - 1/y) + lambda*log(xy)/y;
    }

    inline T Exx (const T x, const T y) const
    {
        assert(x>0 && y>0);
        return mu*(1 + 1/(x*x)) + lambda*(1 - log(xy))/(x*x);
    }

    inline T Eyy (const T x, const T y) const
    {
        assert(x>0 && y>0);
        return mu*(1 + 1/(y*y)) + lambda*(1 - log(xy))/(y*y);
    }

    inline T Exy (const T x, const T y) const
    {
        assert(x>0 && y>0);
        return lambda/(x*y);
    }

    inline T Exxy (const T x, const T y) const
    {
        assert(x>0 && y>0);
        return -lambda/(x*x*y);
    }

    inline T Exyy (const T x, const T y) const
    {
        assert(x>0 && y>0);
        return -lambda/(x*y*y);
    }

};
}
#endif

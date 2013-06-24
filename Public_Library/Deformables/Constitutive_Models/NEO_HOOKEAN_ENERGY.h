//#####################################################################
// Copyright 2011, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NEO_HOOKEAN_ENERGY
//#####################################################################
#ifndef __NEO_HOOKEAN_ENERGY__
#define __NEO_HOOKEAN_ENERGY__

namespace PhysBAM{

template<class T>
class NEO_HOOKEAN_ENERGY
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
    
    // 2D Energy

    inline T E (const T x, const T y) const
    {
        assert(x>0 && y>0);
        return mu*((T)0.5*x*x + (T)0.5*y*y - 1 - log(x*y)) + (T)0.5*lambda*sqr(log(x*y)); 
    }

    inline T Ex (const T x, const T y) const
    {
        assert(x>0 && y>0);
        return mu*(x - 1/x) + lambda*log(x*y)/x;
    }

    inline T Ey (const T x, const T y) const
    {
        assert(x>0 && y>0);
        return mu*(y - 1/y) + lambda*log(x*y)/y;
    }

    inline T Exx (const T x, const T y) const
    {
        assert(x>0 && y>0);
        return mu*(1 + 1/(x*x)) + lambda*(1 - log(x*y))/(x*x);
    }

    inline T Eyy (const T x, const T y) const
    {
        assert(x>0 && y>0);
        return mu*(1 + 1/(y*y)) + lambda*(1 - log(x*y))/(y*y);
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

    // 3D Energy

    inline T E (const T x, const T y, const T z) const
    {
        assert(x>0 && y>0 && z>0);
        T log_J = log(x*y*z);
        return mu*(0.5*(sqr(x) + sqr(y) + sqr(z) - 3) - log_J) + 0.5*lambda*sqr(log_J); 
    }
};
}
#endif

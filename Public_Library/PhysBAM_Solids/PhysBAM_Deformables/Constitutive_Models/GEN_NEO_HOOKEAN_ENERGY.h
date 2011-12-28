//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GEN_NEO_HOOKEAN_ENERGY
//#####################################################################
#ifndef __GEN_NEO_HOOKEAN_ENERGY__
#define __GEN_NEO_HOOKEAN_ENERGY__
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/GENERAL_ENERGY.h>

namespace PhysBAM{

template<class T>
class GEN_NEO_HOOKEAN_ENERGY: public GENERAL_ENERGY<T>
{
public:
    T mu;
    T lambda;

    GEN_NEO_HOOKEAN_ENERGY () {}
    ~GEN_NEO_HOOKEAN_ENERGY () {}

    void Initialize(const T input_mu, const T input_lambda)
    {
        mu = input_mu;
        lambda = input_lambda;
    }

    T E(T x, T y, int simplex) const
    {
        T L=log(x*y);
        return mu*((T)0.5*x*x + (T)0.5*y*y - 1 - L) + (T)0.5*lambda*sqr(L); 
    }

    T Ex(T x, T y, int simplex) const
    {
        return mu*(x - 1/x) + lambda*log(x*y)/x;
    }

    T Exx(T x, T y, int simplex) const
    {
        return mu*(1 + 1/(x*x)) + lambda*(1 - log(x*y))/(x*x);
    }

    T Exy(T x, T y, int simplex) const
    {
        return lambda/(x*y);
    }

    T Exxx(T x, T y, int simplex) const
    {
        return -(2*mu + lambda*(3 - 2*log(x*y)))/(x*x*x);
    }

    T Exxy(T x, T y, int simplex) const
    {
        return -lambda/(x*x*y);
    }

    T Exxxx(T x, T y, int simplex) const
    {
        return (6*mu + lambda*(11 - 6*log(x*y)))/(x*x*x*x);
    }

    T Exxxy(T x, T y, int simplex) const
    {
        return 2*lambda/(x*x*x*y);
    }

    T Exxyy(T x, T y, int simplex) const
    {
        return lambda/(x*x*y*y);
    }

    T Ex_Ey_x_y(T x, T y, int simplex) const
    {
        return mu+(mu-lambda*log(x*y))/(x*y);
    }

    T E(T x, T y, T z, int simplex) const
    {
        T L=log(x*y*z);
        return mu*((T)0.5*x*x + (T)0.5*y*y + (T)0.5*z*z - (T)1.5 - L) + (T)0.5*lambda*sqr(L);
    }

    T Ex(T x, T y, T z, int simplex) const
    {
        T L=log(x*y*z);
        return mu*(x - 1/x) + lambda*L/x;
    }

    T Exx(T x, T y, T z, int simplex) const
    {
        T L=log(x*y*z);
        return mu + (mu - lambda*L + lambda)/sqr(x);
    }

    T Exy(T x, T y, T z, int simplex) const
    {
        return lambda/(x*y);
    }

    T Exxx(T x, T y, T z, int simplex) const
    {
        T L=log(x*y*z);
        return (-2*mu + 2*lambda*L - 3*lambda)/cube(x);
    }

    T Exxy(T x, T y, T z, int simplex) const
    {
        return -lambda/(x*x*y);
    }

    T Exyz(T x, T y, T z, int simplex) const
    {
        return 0;
    }

    T Exxxx(T x, T y, T z, int simplex) const
    {
        return (6*mu - 6*lambda*log(x*y*z) + 11*lambda)/sqr(sqr(x));
    }

    T Exxxy(T x, T y, T z, int simplex) const
    {
        return 2*lambda/(x*x*x*y);
    }

    T Exxyy(T x, T y, T z, int simplex) const
    {
        return lambda/(x*x*y*y);
    }

    T Exxyz(T x, T y, T z, int simplex) const
    {
        return 0;
    }

    T Ex_Ey_x_y(T x, T y, T z, int simplex) const
    {
        T L=log(x*y*z);
        return mu+(mu-lambda*L)/(x*y);
    }

    T Exz_Eyz_x_y(T x, T y, T z, int simplex) const
    {
        return -lambda/(x*y*z);
    }
};
}
#endif

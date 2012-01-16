//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NEO_J_INTERP_ENERGY
//#####################################################################
#ifndef __NEO_J_INTERP_ENERGY__
#define __NEO_J_INTERP_ENERGY__
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/GENERAL_ENERGY.h>

namespace PhysBAM{

template<class T>
class NEO_J_INTERP_ENERGY:public GENERAL_ENERGY<T>
{
public:
    T mu;
    T lambda;
    T J_min,J_max;
    T la_min;

    NEO_J_INTERP_ENERGY(T Ja=(T).4,T Jb=(T).6,T la=(T)-.1)
        :J_min(Ja),J_max(Jb),la_min(la)
    {}

    ~NEO_J_INTERP_ENERGY() {}

    void Initialize(const T input_mu, const T input_lambda)
    {
        mu = input_mu;
        lambda = input_lambda;
    }

    T f(T x) const {return x*x*(3-2*x);}
    T df(T x) const {return 6*x*(1-x);}
    T ddf(T x) const {return 6*(1-2*x);}
    T dddf(T x) const {return -12;}
    // T f(T x) const {return x*x*x;}
    // T df(T x) const {return 3*x*x;}
    // T ddf(T x) const {return 6*x;}
    // T dddf(T x) const {return 6;}

    T La(T J) const
    {
        if(J>J_max) return lambda;
        if(J<J_min) return la_min;
        return f((J-J_min)/(J_max-J_min))*(lambda-la_min)+la_min;
    }

    T dLa(T J) const
    {
        if(J>J_max) return 0;
        if(J<J_min) return 0;
        return df((J-J_min)/(J_max-J_min))*(lambda-la_min)/(J_max-J_min);
    }

    T ddLa(T J) const
    {
        if(J>J_max) return 0;
        if(J<J_min) return 0;
        return ddf((J-J_min)/(J_max-J_min))*(lambda-la_min)/sqr(J_max-J_min);
    }

    T dddLa(T J) const
    {
        if(J>J_max) return 0;
        if(J<J_min) return 0;
        return dddf((J-J_min)/(J_max-J_min))*(lambda-la_min)/cube(J_max-J_min);
    }

    T E(T x, T y, int simplex) const
    {
        T J=x*y,L=log(J);
        return mu*((T)0.5*x*x + (T)0.5*y*y - 1 - L) + (T)0.5*La(J)*sqr(L); 
    }

    T Ex(T x, T y, int simplex) const
    {
        T J=x*y,L=log(J);
        return mu*(x-1/x)+(T)0.5*dLa(J)*y*sqr(L)+La(J)*L/x;
    }

    T Exx(T x, T y, int simplex) const
    {
        T J=x*y,L=log(J);
        return mu*(1+1/sqr(x))+(T)0.5*ddLa(J)*sqr(y)*sqr(L)+2*dLa(J)*y*L/x+La(J)/sqr(x)-La(J)*L/sqr(x);
    }

    T Exy(T x, T y, int simplex) const
    {
        T J=x*y,L=log(J);
        return (T)0.5*ddLa(J)*J*sqr(L)+(T)0.5*dLa(J)*sqr(L)+2*dLa(J)*L+La(J)/J;
    }

    T Exxx(T x, T y, int simplex) const
    {
        PHYSBAM_FATAL_ERROR();
    }

    T Exxy(T x, T y, int simplex) const
    {
        T J=x*y,L=log(J);
        return (T)0.5*dddLa(J)*J*y*sqr(L)+ddLa(J)*y*sqr(L)+3*ddLa(J)*y*L+dLa(J)*L/x+3*dLa(J)/x-La(J)/(J*x);
    }

    T Exxxx(T x, T y, int simplex) const
    {
        PHYSBAM_FATAL_ERROR();
    }

    T Exxxy(T x, T y, int simplex) const
    {
        PHYSBAM_FATAL_ERROR();
    }

    T Exxyy(T x, T y, int simplex) const
    {
        PHYSBAM_FATAL_ERROR();
    }

    T Ex_Ey_x_y(T x, T y, int simplex) const
    {
        T J=x*y,L=log(J);
        return mu-(T)0.5*dLa(J)*sqr(L)+mu/J-La(J)*L/J;
    }

    T Exx_Exy_x_y(T x, T y, int simplex) const
    {
        PHYSBAM_FATAL_ERROR();
    }

    T Exxx_Exxy_x_y(T x, T y, int simplex) const
    {
        PHYSBAM_FATAL_ERROR();
    }

    T Exxy_Exyy_x_y(T x, T y, int simplex) const
    {
        PHYSBAM_FATAL_ERROR();
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
        PHYSBAM_FATAL_ERROR();
    }

    T Exxy(T x, T y, T z, int simplex) const
    {
        return -lambda/(x*x*y);
    }

    T Exyz(T x, T y, T z, int simplex) const
    {
        return 0;
    }

    T Exxyz(T x, T y, T z, int simplex) const
    {
        return 0;
    }

    T Exxxx(T x, T y, T z, int simplex) const
    {
        PHYSBAM_FATAL_ERROR();
    }

    T Exxxy(T x, T y, T z, int simplex) const
    {
        PHYSBAM_FATAL_ERROR();
    }

    T Exxyy(T x, T y, T z, int simplex) const
    {
        PHYSBAM_FATAL_ERROR();
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

    T Exx_Exy_x_y(T x, T y, T z, int simplex) const
    {
        PHYSBAM_FATAL_ERROR();
    }

    T Exxx_Exxy_x_y(T x, T y, T z, int simplex) const
    {
        PHYSBAM_FATAL_ERROR();
    }

    T Exxy_Exyy_x_y(T x, T y, T z, int simplex) const
    {
        PHYSBAM_FATAL_ERROR();
    }

    T Exzz_Eyzz_x_y(T x, T y, T z, int simplex) const
    {
        PHYSBAM_FATAL_ERROR();
    }

    T Exxz_Exyz_x_y(T x, T y, T z, int simplex) const
    {
        PHYSBAM_FATAL_ERROR();
    }
};
}
#endif

//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GENERAL_ENERGY
//#####################################################################
#ifndef __GENERAL_ENERGY__
#define __GENERAL_ENERGY__

namespace PhysBAM{

template<class TV>
class GENERAL_ENERGY
{
public:
    typedef typename TV::SCALAR T;
    T mu;
    T lambda;

    GENERAL_ENERGY () {}
    ~GENERAL_ENERGY () {}

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

    T Exxy(T x, T y, int simplex) const
    {
        return -lambda/(x*x*y);
    }

    T Ex_Ey_x_y(T x, T y, int simplex) const
    {
        return mu+(mu-lambda*log(x*y))/(x*y);
    }

    T Ey(T x, T y, int simplex) const {return Ex(y,x,simplex);}
    T Eyy(T x, T y, int simplex) const {return Exx(y,x,simplex);}
    T Exyy(T x, T y, int simplex) const {return Exxy(y,x,simplex);}

// TODO: Fix these.
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

    T Ex_Ey_x_y(T x, T y, T z, int simplex) const
    {
        T L=log(x*y*z);
        return mu+(mu-lambda*L)/(x*y);
    }

    T Exz_Eyz_x_y(T x, T y, T z, int simplex) const
    {
        return -lambda/(x*y*z);
    }

    T Ey(T x, T y, T z, int simplex) const {return Ex(y,x,z,simplex);}
    T Ez(T x, T y, T z, int simplex) const {return Ex(z,y,x,simplex);}

    T Eyy(T x, T y, T z, int simplex) const {return Exx(y,x,z,simplex);}
    T Ezz(T x, T y, T z, int simplex) const {return Exx(z,y,x,simplex);}
    T Exz(T x, T y, T z, int simplex) const {return Exy(x,z,y,simplex);}
    T Eyz(T x, T y, T z, int simplex) const {return Exy(y,z,x,simplex);}

    T Ex_Ez_x_z(T x, T y, T z, int simplex) const {return Ex_Ey_x_y(x,z,y,simplex);}
    T Ey_Ez_y_z(T x, T y, T z, int simplex) const {return Ex_Ey_x_y(y,z,x,simplex);}

    T Exy_Eyz_x_z(T x, T y, T z, int simplex) const {return Exz_Eyz_x_y(x,z,y,simplex);}
    T Exy_Exz_y_z(T x, T y, T z, int simplex) const {return Exz_Eyz_x_y(y,z,x,simplex);}

    T Exyy(T x, T y, T z, int simplex) const {return Exxy(y,x,z,simplex);}
    T Exxz(T x, T y, T z, int simplex) const {return Exxy(x,z,y,simplex);}
    T Exzz(T x, T y, T z, int simplex) const {return Exxy(z,x,y,simplex);}
    T Eyyz(T x, T y, T z, int simplex) const {return Exxy(y,z,x,simplex);}
    T Eyzz(T x, T y, T z, int simplex) const {return Exxy(z,y,x,simplex);}

    T Exyyz(T x, T y, T z, int simplex) const {return Exxyz(y,x,z,simplex);}
    T Exyzz(T x, T y, T z, int simplex) const {return Exxyz(z,y,x,simplex);}
};
}
#endif

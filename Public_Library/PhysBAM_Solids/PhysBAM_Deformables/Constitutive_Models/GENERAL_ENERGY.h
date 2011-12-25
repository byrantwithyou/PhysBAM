//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GENERAL_ENERGY
//#####################################################################
#ifndef __GENERAL_ENERGY__
#define __GENERAL_ENERGY__

namespace PhysBAM{

template<class T>
class GENERAL_ENERGY
{
public:

    GENERAL_ENERGY () {}
    virtual ~GENERAL_ENERGY () {}

    virtual void Initialize(const T input_mu, const T input_lambda)=0;

    virtual T E(T x, T y, int simplex) const=0;
    virtual T Ex(T x, T y, int simplex) const=0;
    virtual T Exx(T x, T y, int simplex) const=0;
    virtual T Exy(T x, T y, int simplex) const=0;
    virtual T Exxy(T x, T y, int simplex) const=0;
    virtual T Ex_Ey_x_y(T x, T y, int simplex) const=0; // (Ex-Ey)/(x-y)

    virtual T E(T x, T y, T z, int simplex) const=0;
    virtual T Ex(T x, T y, T z, int simplex) const=0;
    virtual T Exx(T x, T y, T z, int simplex) const=0;
    virtual T Exy(T x, T y, T z, int simplex) const=0;
    virtual T Exxy(T x, T y, T z, int simplex) const=0;
    virtual T Exyz(T x, T y, T z, int simplex) const=0;
    virtual T Exxyz(T x, T y, T z, int simplex) const=0;
    virtual T Ex_Ey_x_y(T x, T y, T z, int simplex) const=0; // (Ex-Ey)/(x-y)
    virtual T Exz_Eyz_x_y(T x, T y, T z, int simplex) const=0; // (Exz-Eyz)/(x-y)

    T Ey(T x, T y, int simplex) const {return Ex(y,x,simplex);}
    T Eyy(T x, T y, int simplex) const {return Exx(y,x,simplex);}
    T Exyy(T x, T y, int simplex) const {return Exxy(y,x,simplex);}

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

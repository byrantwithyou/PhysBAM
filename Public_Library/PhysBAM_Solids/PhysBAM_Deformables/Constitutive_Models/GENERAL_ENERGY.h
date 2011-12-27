//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GENERAL_ENERGY
//#####################################################################
#ifndef __GENERAL_ENERGY__
#define __GENERAL_ENERGY__
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

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
    virtual T Exxx(T x, T y, int simplex) const=0;
    virtual T Exxy(T x, T y, int simplex) const=0;
    virtual T Ex_Ey_x_y(T x, T y, int simplex) const=0; // (Ex-Ey)/(x-y)

    virtual T E(T x, T y, T z, int simplex) const=0;
    virtual T Ex(T x, T y, T z, int simplex) const=0;
    virtual T Exx(T x, T y, T z, int simplex) const=0;
    virtual T Exy(T x, T y, T z, int simplex) const=0;
    virtual T Exxx(T x, T y, T z, int simplex) const=0;
    virtual T Exxy(T x, T y, T z, int simplex) const=0;
    virtual T Exyz(T x, T y, T z, int simplex) const=0;
    virtual T Exxyz(T x, T y, T z, int simplex) const=0;
    virtual T Ex_Ey_x_y(T x, T y, T z, int simplex) const=0; // (Ex-Ey)/(x-y)
    virtual T Exz_Eyz_x_y(T x, T y, T z, int simplex) const=0; // (Exz-Eyz)/(x-y)

    T Ey(T x, T y, int simplex) const {return Ex(y,x,simplex);}
    T Eyy(T x, T y, int simplex) const {return Exx(y,x,simplex);}
    T Exyy(T x, T y, int simplex) const {return Exxy(y,x,simplex);}
    T Eyyy(T x, T y, int simplex) const {return Exxx(y,x,simplex);}

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

    T Eyyy(T x, T y, T z, int simplex) const {return Exxx(y,x,z,simplex);}
    T Ezzz(T x, T y, T z, int simplex) const {return Exxx(z,x,y,simplex);}
    T Exyy(T x, T y, T z, int simplex) const {return Exxy(y,x,z,simplex);}
    T Exxz(T x, T y, T z, int simplex) const {return Exxy(x,z,y,simplex);}
    T Exzz(T x, T y, T z, int simplex) const {return Exxy(z,x,y,simplex);}
    T Eyyz(T x, T y, T z, int simplex) const {return Exxy(y,z,x,simplex);}
    T Eyzz(T x, T y, T z, int simplex) const {return Exxy(z,y,x,simplex);}

    T Exyyz(T x, T y, T z, int simplex) const {return Exxyz(y,x,z,simplex);}
    T Exyzz(T x, T y, T z, int simplex) const {return Exxyz(z,y,x,simplex);}

    void Test(T x, T y, int simplex) const;
    void Test(T x, T y, T z, int simplex) const;

    T E(const VECTOR<T,2>& f, int simplex) const {return E(f.x,f.y,simplex);}

    T E(const VECTOR<T,3>& f, int simplex) const {return E(f.x,f.y,f.z,simplex);}

    VECTOR<T,2> dE(const VECTOR<T,2>& f, int simplex) const {return VECTOR<T,2>(Ex(f.x,f.y,simplex),Ey(f.x,f.y,simplex));}

    VECTOR<T,3> dE(const VECTOR<T,3>& f, int simplex) const {return VECTOR<T,3>(Ex(f.x,f.y,f.z,simplex),Ey(f.x,f.y,f.z,simplex),Ez(f.x,f.y,f.z,simplex));}

    SYMMETRIC_MATRIX<T,2> ddE(const VECTOR<T,2>& f, int simplex) const
    {return SYMMETRIC_MATRIX<T,2>(Exx(f.x,f.y,simplex),Exy(f.x,f.y,simplex),Eyy(f.x,f.y,simplex));}

    SYMMETRIC_MATRIX<T,3> ddE(const VECTOR<T,3>& f, int simplex) const
    {return SYMMETRIC_MATRIX<T,3>(Exx(f.x,f.y,f.z,simplex),Exy(f.x,f.y,f.z,simplex),Exz(f.x,f.y,f.z,simplex),
        Eyy(f.x,f.y,f.z,simplex),Eyz(f.x,f.y,f.z,simplex),Ezz(f.x,f.y,f.z,simplex));}

    void dddE(const VECTOR<T,2>& f, int simplex,SYMMETRIC_MATRIX<T,2>* sm) const
    {
        sm[0].x11=Exxx(f.x,f.y,simplex);
        sm[1].x11=sm[0].x21=Exxy(f.x,f.y,simplex);
        sm[1].x21=sm[0].x22=Exyy(f.x,f.y,simplex);
        sm[1].x22=Eyyy(f.x,f.y,simplex);
    }

    void dddE(const VECTOR<T,3>& f, int simplex,SYMMETRIC_MATRIX<T,3>* sm) const
    {
        sm[0].x11=Exxx(f.x,f.y,f.z,simplex);
        sm[0].x21=sm[1].x11=Exxy(f.x,f.y,f.z,simplex);
        sm[0].x31=sm[2].x11=Exxz(f.x,f.y,f.z,simplex);
        sm[0].x22=sm[1].x21=Exyy(f.x,f.y,f.z,simplex);
        sm[0].x32=sm[1].x31=sm[2].x21=Exyz(f.x,f.y,f.z,simplex);
        sm[0].x33=sm[2].x31=Exzz(f.x,f.y,f.z,simplex);
        sm[1].x22=Eyyy(f.x,f.y,f.z,simplex);
        sm[1].x32=sm[2].x22=Eyyz(f.x,f.y,f.z,simplex);
        sm[1].x33=sm[2].x32=Eyzz(f.x,f.y,f.z,simplex);
        sm[2].x33=Ezzz(f.x,f.y,f.z,simplex);
    }
};
}
#endif

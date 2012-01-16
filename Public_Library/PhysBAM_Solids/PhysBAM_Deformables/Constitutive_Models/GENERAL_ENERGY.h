//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GENERAL_ENERGY
//#####################################################################
#ifndef __GENERAL_ENERGY__
#define __GENERAL_ENERGY__
#include <PhysBAM_Tools/Matrices/MATRIX.h>
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
    virtual T Exxxx(T x, T y, int simplex) const=0;
    virtual T Exxxy(T x, T y, int simplex) const=0;
    virtual T Exxyy(T x, T y, int simplex) const=0;
    virtual T Ex_Ey_x_y(T x, T y, int simplex) const=0; // (Ex-Ey)/(x-y)
    virtual T Exxy_Exyy_x_y(T x, T y, int simplex) const=0; // (Exxy-Exyy)/(x-y)
    virtual T Exx_Eyy_x_y(T x, T y, int simplex) const=0;
    virtual T Exxx_Eyyy_x_y(T x, T y, int simplex) const=0;

    virtual T E(T x, T y, T z, int simplex) const=0;
    virtual T Ex(T x, T y, T z, int simplex) const=0;
    virtual T Exx(T x, T y, T z, int simplex) const=0;
    virtual T Exy(T x, T y, T z, int simplex) const=0;
    virtual T Exxx(T x, T y, T z, int simplex) const=0;
    virtual T Exxy(T x, T y, T z, int simplex) const=0;
    virtual T Exyz(T x, T y, T z, int simplex) const=0;
    virtual T Exxxx(T x, T y, T z, int simplex) const=0;
    virtual T Exxxy(T x, T y, T z, int simplex) const=0;
    virtual T Exxyy(T x, T y, T z, int simplex) const=0;
    virtual T Exxyz(T x, T y, T z, int simplex) const=0;

    virtual T Ex_Ey_x_y(T x, T y, T z, int simplex) const=0; // (Ex-Ey)/(x-y)
    virtual T Exz_Eyz_x_y(T x, T y, T z, int simplex) const=0; // (Exz-Eyz)/(x-y)
    virtual T Exxy_Exyy_x_y(T x, T y, T z, int simplex) const=0; // (Exxy-Exyy)/(x-y)
    virtual T Exzz_Eyzz_x_y(T x, T y, T z, int simplex) const=0; // (Exzz-Eyzz)/(x-y)
    virtual T Exx_Eyy_x_y(T x, T y, T z, int simplex) const=0;
    virtual T Exxx_Eyyy_x_y(T x, T y, T z, int simplex) const=0;
    virtual T Exxz_Eyyz_x_y(T x, T y, T z, int simplex) const=0;

    T Ey(T x, T y, int simplex) const {return Ex(y,x,simplex);}
    T Eyy(T x, T y, int simplex) const {return Exx(y,x,simplex);}
    T Exyy(T x, T y, int simplex) const {return Exxy(y,x,simplex);}
    T Eyyy(T x, T y, int simplex) const {return Exxx(y,x,simplex);}
    T Exyyy(T x, T y, int simplex) const {return Exxxy(y,x,simplex);}
    T Eyyyy(T x, T y, int simplex) const {return Exxxx(y,x,simplex);}

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

    T Exyy_Eyyz_x_z(T x, T y, T z, int simplex) const {return Exzz_Eyzz_x_y(x,z,y,simplex);}
    T Exxz_Exzz_x_z(T x, T y, T z, int simplex) const {return Exxy_Exyy_x_y(x,z,y,simplex);}
    T Exxy_Exxz_y_z(T x, T y, T z, int simplex) const {return Exzz_Eyzz_x_y(z,y,x,simplex);}
    T Eyyz_Eyzz_y_z(T x, T y, T z, int simplex) const {return Exxy_Exyy_x_y(z,y,x,simplex);}

    T Eyyy(T x, T y, T z, int simplex) const {return Exxx(y,x,z,simplex);}
    T Ezzz(T x, T y, T z, int simplex) const {return Exxx(z,x,y,simplex);}
    T Exyy(T x, T y, T z, int simplex) const {return Exxy(y,x,z,simplex);}
    T Exxz(T x, T y, T z, int simplex) const {return Exxy(x,z,y,simplex);}
    T Exzz(T x, T y, T z, int simplex) const {return Exxy(z,x,y,simplex);}
    T Eyyz(T x, T y, T z, int simplex) const {return Exxy(y,z,x,simplex);}
    T Eyzz(T x, T y, T z, int simplex) const {return Exxy(z,y,x,simplex);}

    T Eyyyy(T x, T y, T z, int simplex) const {return Exxxx(y,x,z,simplex);}
    T Ezzzz(T x, T y, T z, int simplex) const {return Exxxx(z,y,x,simplex);}
    T Exyyy(T x, T y, T z, int simplex) const {return Exxxy(y,x,z,simplex);}
    T Exxxz(T x, T y, T z, int simplex) const {return Exxxy(x,z,y,simplex);}
    T Exzzz(T x, T y, T z, int simplex) const {return Exxxy(z,x,y,simplex);}
    T Eyyyz(T x, T y, T z, int simplex) const {return Exxxy(y,z,x,simplex);}
    T Eyzzz(T x, T y, T z, int simplex) const {return Exxxy(z,y,x,simplex);}
    T Exxzz(T x, T y, T z, int simplex) const {return Exxyy(x,z,y,simplex);}
    T Eyyzz(T x, T y, T z, int simplex) const {return Exxyy(y,z,x,simplex);}
    T Exyyz(T x, T y, T z, int simplex) const {return Exxyz(y,x,z,simplex);}
    T Exyzz(T x, T y, T z, int simplex) const {return Exxyz(z,y,x,simplex);}

    T Eyy_Ezz_y_z(T x, T y, T z, int simplex) const {return Exx_Eyy_x_y(z,y,x,simplex);}
    T Exx_Ezz_x_z(T x, T y, T z, int simplex) const {return Exx_Eyy_x_y(x,z,y,simplex);}
    T Exyy_Exzz_y_z(T x, T y, T z, int simplex) const {return Exxz_Eyyz_x_y(y,z,x,simplex);}
    T Exxy_Eyzz_x_z(T x, T y, T z, int simplex) const {return Exxz_Eyyz_x_y(x,z,y,simplex);}
    T Eyyy_Ezzz_y_z(T x, T y, T z, int simplex) const {return Exxx_Eyyy_x_y(z,y,x,simplex);}
    T Exxx_Ezzz_x_z(T x, T y, T z, int simplex) const {return Exxx_Eyyy_x_y(x,z,y,simplex);}

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

    void ddddE(const VECTOR<T,2>& f, int simplex,SYMMETRIC_MATRIX<T,2>* sm) const
    {
        sm[0+2*0].x11=Exxxx(f.x,f.y,simplex);
        sm[0+2*0].x21=sm[1+2*0].x11=sm[0+2*1].x11=Exxxy(f.x,f.y,simplex);
        sm[1+2*1].x11=sm[0+2*1].x21=sm[1+2*0].x21=sm[0+2*0].x22=Exxyy(f.x,f.y,simplex);
        sm[1+2*1].x21=sm[1+2*0].x22=sm[0+2*1].x22=Exyyy(f.x,f.y,simplex);
        sm[1+2*1].x22=Eyyyy(f.x,f.y,simplex);
    }

    void ddddE(const VECTOR<T,3>& f, int simplex,SYMMETRIC_MATRIX<T,3>* sm) const
    {
        sm[0+3*0].x11=Exxxx(f.x,f.y,f.z,simplex);
        sm[1+3*0].x11=sm[0+3*0].x21=sm[0+3*1].x11=Exxxy(f.x,f.y,f.z,simplex);
        sm[2+3*0].x11=sm[0+3*2].x11=sm[0+3*0].x31=Exxxz(f.x,f.y,f.z,simplex);
        sm[0+3*1].x21=sm[1+3*0].x21=sm[0+3*0].x22=sm[1+3*1].x11=Exxyy(f.x,f.y,f.z,simplex);
        sm[1+3*2].x11=sm[2+3*0].x21=sm[0+3*2].x21=sm[1+3*0].x31=sm[0+3*0].x32=sm[0+3*1].x31=sm[2+3*1].x11=Exxyz(f.x,f.y,f.z,simplex);
        sm[0+3*0].x33=sm[0+3*2].x31=sm[2+3*0].x31=sm[2+3*2].x11=Exxzz(f.x,f.y,f.z,simplex);
        sm[0+3*1].x22=sm[1+3*0].x22=sm[1+3*1].x21=Exyyy(f.x,f.y,f.z,simplex);
        sm[0+3*1].x32=sm[0+3*2].x22=sm[1+3*0].x32=sm[1+3*1].x31=sm[1+3*2].x21=sm[2+3*0].x22=sm[2+3*1].x21=Exyyz(f.x,f.y,f.z,simplex);
        sm[0+3*1].x33=sm[0+3*2].x32=sm[1+3*0].x33=sm[1+3*2].x31=sm[2+3*0].x32=sm[2+3*1].x31=sm[2+3*2].x21=Exyzz(f.x,f.y,f.z,simplex);
        sm[0+3*2].x33=sm[2+3*0].x33=sm[2+3*2].x31=Exzzz(f.x,f.y,f.z,simplex);
        sm[1+3*1].x22=Eyyyy(f.x,f.y,f.z,simplex);
        sm[1+3*1].x32=sm[1+3*2].x22=sm[2+3*1].x22=Eyyyz(f.x,f.y,f.z,simplex);
        sm[1+3*1].x33=sm[1+3*2].x32=sm[2+3*1].x32=sm[2+3*2].x22=Eyyzz(f.x,f.y,f.z,simplex);
        sm[1+3*2].x33=sm[2+3*1].x33=sm[2+3*2].x32=Eyzzz(f.x,f.y,f.z,simplex);
        sm[2+3*2].x33=Ezzzz(f.x,f.y,f.z,simplex);
    }

    void Compute_it(const VECTOR<T,2>& f,int simplex,VECTOR<T,1>& g_is,VECTOR<T,0>& H_xit,VECTOR<T,1>& H_iitt,VECTOR<T,0>& T_xxit,VECTOR<T,0>& T_xiitt,VECTOR<T,1>& T_iiittt,VECTOR<T,1>& T_itit) const
    {
        g_is(1)=Ex_Ey_x_y(f.x,f.y,simplex);
        H_iitt(1)=Exx_Eyy_x_y(f.x,f.y,simplex);
        T_iiittt(1)=Exxx_Eyyy_x_y(f.x,f.y,simplex);
        T_itit(1)=Exxy_Exyy_x_y(f.x,f.y,simplex);
    }

    void Compute_it(const VECTOR<T,3>& f,int simplex,VECTOR<T,3>& g_is,VECTOR<T,3>& H_xit,VECTOR<T,3>& H_iitt,VECTOR<T,3>& T_xxit,VECTOR<T,3>& T_xiitt,VECTOR<T,3>& T_iiittt,VECTOR<T,3>& T_itit) const
    {
        g_is(1)=Ey_Ez_y_z(f.x,f.y,f.z,simplex);
        g_is(2)=Ex_Ez_x_z(f.x,f.y,f.z,simplex);
        g_is(3)=Ex_Ey_x_y(f.x,f.y,f.z,simplex);
        H_iitt(1)=Eyy_Ezz_y_z(f.x,f.y,f.z,simplex);
        H_iitt(2)=Exx_Ezz_x_z(f.x,f.y,f.z,simplex);
        H_iitt(3)=Exx_Eyy_x_y(f.x,f.y,f.z,simplex);
        H_xit(1)=Exy_Exz_y_z(f.x,f.y,f.z,simplex);
        H_xit(2)=Exy_Eyz_x_z(f.x,f.y,f.z,simplex);
        H_xit(3)=Exz_Eyz_x_y(f.x,f.y,f.z,simplex);
        T_xxit(1)=Exxy_Exxz_y_z(f.x,f.y,f.z,simplex);
        T_xxit(2)=Exyy_Eyyz_x_z(f.x,f.y,f.z,simplex);
        T_xxit(3)=Exzz_Eyzz_x_y(f.x,f.y,f.z,simplex);
        T_xiitt(1)=Exyy_Exzz_y_z(f.x,f.y,f.z,simplex);
        T_xiitt(2)=Exxy_Eyzz_x_z(f.x,f.y,f.z,simplex);
        T_xiitt(3)=Exxz_Eyyz_x_y(f.x,f.y,f.z,simplex);
        T_iiittt(1)=Eyyy_Ezzz_y_z(f.x,f.y,f.z,simplex);
        T_iiittt(2)=Exxx_Ezzz_x_z(f.x,f.y,f.z,simplex);
        T_iiittt(3)=Exxx_Eyyy_x_y(f.x,f.y,f.z,simplex);
        T_itit(1)=Eyyz_Eyzz_y_z(f.x,f.y,f.z,simplex);
        T_itit(2)=Exxz_Exzz_x_z(f.x,f.y,f.z,simplex);
        T_itit(3)=Exxy_Exyy_x_y(f.x,f.y,f.z,simplex);
    }
};
}
#endif


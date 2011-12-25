//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/GENERAL_ENERGY.h>
#include <cstdio>
using namespace PhysBAM;
template<class T>
static void Test2d(const char* sf, const char* var, const GENERAL_ENERGY<T>& ge,T (GENERAL_ENERGY<T>::*f)(T,T,int) const,T (GENERAL_ENERGY<T>::*df)(T,T,int) const,T x,T y,int simplex,T dx,T dy)
{
    T g0=(ge.*f)(x,y,simplex);
    T g1=(ge.*f)(x+dx,y+dy,simplex);
    T dg0=(ge.*df)(x,y,simplex);
    T dg1=(ge.*df)(x+dx,y+dy,simplex);
    T av=(dg0+dg1)/2;
    T dif=(g1-g0)/(dx+dy);
    char buff[1000];
    sprintf(buff, "%6s%s %8.5f %8.5f (%8.5f)\n", sf, var, av, dif, (av-dif));
    LOG::cout<<buff;
}
#define TEST2D(f,a) {T dx=0,dy=0;d##a=e;Test2d(#f,#a,*this,&GENERAL_ENERGY<T>::f,&GENERAL_ENERGY<T>::f##a,x,y,simplex,dx,dy);}
//#####################################################################
// Function Test
//#####################################################################
template<class T> void GENERAL_ENERGY<T>::
Test(T x, T y, int simplex) const
{
    PHYSBAM_ASSERT(sizeof(T)==sizeof(double));
    T e=(T)1e-6;
    char buff[1000];
    sprintf(buff, "F: %8.5f %8.5f (J=%8.5f)", x, y, x*y);
    LOG::cout<<buff<<std::endl;
    TEST2D(E,x);
    TEST2D(E,y);
    TEST2D(Ex,x);
    TEST2D(Ex,y);
    TEST2D(Ey,y);
    TEST2D(Exx,y);
    TEST2D(Exy,y);
}
//#####################################################################
// Function Test
//#####################################################################
template<class T> void GENERAL_ENERGY<T>::
Test(T x, T y, T z, int simplex) const
{
    PHYSBAM_ASSERT(sizeof(T)==sizeof(double));
/*
    virtual T E(T x, T y, T z, int simplex) const=0;
    virtual T Ex(T x, T y, T z, int simplex) const=0;
    virtual T Exx(T x, T y, T z, int simplex) const=0;
    virtual T Exy(T x, T y, T z, int simplex) const=0;
    virtual T Exxy(T x, T y, T z, int simplex) const=0;
    virtual T Exyz(T x, T y, T z, int simplex) const=0;
    virtual T Exxyz(T x, T y, T z, int simplex) const=0;
    virtual T Ex_Ey_x_y(T x, T y, T z, int simplex) const=0; // (Ex-Ey)/(x-y)
    virtual T Exz_Eyz_x_y(T x, T y, T z, int simplex) const=0; // (Exz-Eyz)/(x-y)

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
*/
}
template class GENERAL_ENERGY<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GENERAL_ENERGY<double>;
#endif

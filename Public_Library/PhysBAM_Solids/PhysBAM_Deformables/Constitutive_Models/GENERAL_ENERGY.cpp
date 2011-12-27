//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/GENERAL_ENERGY.h>
#include <cstdio>
using namespace PhysBAM;
//#####################################################################
// Function Test2d
//#####################################################################
template<class T> static void
Test2d(const char* sf, const char* var, const GENERAL_ENERGY<T>& ge,T (GENERAL_ENERGY<T>::*f)(T,T,int) const,T (GENERAL_ENERGY<T>::*df)(T,T,int) const,
    T x,T y,int simplex,T dx,T dy)
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
    TEST2D(Exx,x);
    TEST2D(Exx,y);
    TEST2D(Exy,y);
    TEST2D(Eyy,y);
}
//#####################################################################
// Function Test23
//#####################################################################
template<class T> static void
Test3d(const char* sf, const char* var, const GENERAL_ENERGY<T>& ge,T (GENERAL_ENERGY<T>::*f)(T,T,T,int) const,T (GENERAL_ENERGY<T>::*df)(T,T,T,int) const,
    T x,T y,T z,int simplex,T dx,T dy,T dz)
{
    T g0=(ge.*f)(x,y,z,simplex);
    T g1=(ge.*f)(x+dx,y+dy,z+dz,simplex);
    T dg0=(ge.*df)(x,y,z,simplex);
    T dg1=(ge.*df)(x+dx,y+dy,z+dz,simplex);
    T av=(dg0+dg1)/2;
    T dif=(g1-g0)/(dx+dy+dz);
    char buff[1000];
    sprintf(buff, "%6s%s %8.5f %8.5f (%8.5f)\n", sf, var, av, dif, (av-dif));
    LOG::cout<<buff;
}
#define TEST3D(f,a) {T dx=0,dy=0,dz=0;d##a=e;Test3d(#f,#a,*this,&GENERAL_ENERGY<T>::f,&GENERAL_ENERGY<T>::f##a,x,y,z,simplex,dx,dy,dz);}
//#####################################################################
// Function Test
//#####################################################################
template<class T> void GENERAL_ENERGY<T>::
Test(T x, T y, T z, int simplex) const
{
    PHYSBAM_ASSERT(sizeof(T)==sizeof(double));
    T e=(T)1e-6;
    char buff[1000];
    sprintf(buff, "F: %8.5f %8.5f %8.5f (J=%8.5f)", x, y, z, x*y*z);
    LOG::cout<<buff<<std::endl;
    TEST3D(E,x);
    TEST3D(E,y);
    TEST3D(E,z);
    TEST3D(Ex,x);
    TEST3D(Ex,y);
    TEST3D(Ey,y);
    TEST3D(Ex,z);
    TEST3D(Ey,z);
    TEST3D(Ez,z);
    TEST3D(Exx,x);
    TEST3D(Exx,y);
    TEST3D(Exy,y);
    TEST3D(Exx,z);
    TEST3D(Exy,z);
    TEST3D(Eyy,y);
    TEST3D(Eyy,z);
    TEST3D(Exz,z);
    TEST3D(Eyz,z);
    TEST3D(Ezz,z);
    TEST3D(Exxy,z);
    TEST3D(Exyy,z);
    TEST3D(Exyz,z);
}
template class GENERAL_ENERGY<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GENERAL_ENERGY<double>;
#endif

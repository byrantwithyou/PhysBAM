//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Neil Molino, Igor Neverov, Duc Nguyen, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RANDOM_NUMBERS
//#####################################################################
// See Bratley and Fox. 1988. Algorithm 659: Implementing Sobol's quasirandom sequence generator. ACM Trans. Math. Softw. 14, 88-100.
//#####################################################################
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/MATRIX_MXN.h>
#include <Tools/Matrices/ROTATION.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Matrices/UPPER_TRIANGULAR_MATRIX_0X0.h>
#include <Tools/Matrices/UPPER_TRIANGULAR_MATRIX_1X1.h>
#include <Tools/Matrices/UPPER_TRIANGULAR_MATRIX_2X2.h>
#include <Tools/Matrices/UPPER_TRIANGULAR_MATRIX_3X3.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Vectors/TWIST.h>
#include <Tools/Vectors/VECTOR.h>
#include <limits>
using ::std::log;
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T,class GENERATOR> RANDOM_NUMBERS<T,GENERATOR>::
RANDOM_NUMBERS(const unsigned int seed)
    :gaussian_iset(0)
{
    Set_Seed(seed);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class GENERATOR> RANDOM_NUMBERS<T,GENERATOR>::
~RANDOM_NUMBERS()
{}
//#####################################################################
// Function Set_Seed
//#####################################################################
template<class T,class GENERATOR> void RANDOM_NUMBERS<T,GENERATOR>::
Set_Seed(const unsigned int seed_input)
{
    random_number_generator.Set_Seed(seed_input);
}
//#####################################################################
// Function Get_Uniform_Integer
//#####################################################################
template<class T,class GENERATOR> int RANDOM_NUMBERS<T,GENERATOR>::
Get_Uniform_Integer(const int a,const int b)
{
    return min(b,(int)(a+(b+1-a)*Get_Number())); // in [a,b]
}
//#####################################################################
// Function Get_Uniform_Number
//#####################################################################
template<class T,class GENERATOR> T RANDOM_NUMBERS<T,GENERATOR>::
Get_Uniform_Number(const T a,const T b)
{
    STATIC_ASSERT((!std::numeric_limits<T>::is_integer));
    return a+(b-a)*Get_Number(); // in [a,b)
}
//#####################################################################
// Function Get_Uniform_Vector
//#####################################################################
template<class T,class GENERATOR> template<int d> VECTOR<T,d> RANDOM_NUMBERS<T,GENERATOR>::
Get_Uniform_Vector(const VECTOR<T,d>& v0,const VECTOR<T,d>& v1)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r(i)=Get_Uniform_Number(v0(i),v1(i));
    return r;
}
//#####################################################################
// Function Get_Uniform_Vector
//#####################################################################
template<class T,class GENERATOR> template<int d> VECTOR<T,d> RANDOM_NUMBERS<T,GENERATOR>::
Get_Uniform_Vector(const T a,const T b)
{
    VECTOR<T,d> r;
    Fill_Uniform(r,a,b);
    return r;
}
//#####################################################################
// Function Fill_Uniform
//#####################################################################
template<class T,class GENERATOR> void RANDOM_NUMBERS<T,GENERATOR>::
Fill_Uniform(T& x,const T a,const T b)
{
    x=Get_Uniform_Number(a,b);
}
//#####################################################################
// Function Fill_Uniform
//#####################################################################
template<class T,class GENERATOR> template<class T_VECTOR> void RANDOM_NUMBERS<T,GENERATOR>::
Fill_Uniform(ARRAY_BASE<T,T_VECTOR>& v,const T a,const T b)
{
    for(int i=0;i<v.Size();i++) v(i)=Get_Uniform_Number(a,b);
}
//#####################################################################
// Function Fill_Uniform_Matrix
//#####################################################################
template<class T,class GENERATOR> template<class T_MATRIX> void RANDOM_NUMBERS<T,GENERATOR>::
Fill_Uniform(MATRIX_BASE<T,T_MATRIX>& m,const T a,const T b)
{
    for(int i=0;i<m.Rows();i++) for(int j=0;j<m.Columns();j++) m(i,j)=Get_Uniform_Number(a,b);
}
//#####################################################################
// Function Fill_Uniform_Matrix
//#####################################################################
template<class T,class GENERATOR> template<int d> void RANDOM_NUMBERS<T,GENERATOR>::
Fill_Uniform(DIAGONAL_MATRIX<T,d>& m,const T a,const T b)
{
    for(int i=0;i<m.Rows();i++) m(i,i)=Get_Uniform_Number(a,b);
}
//#####################################################################
// Function Fill_Uniform_Matrix
//#####################################################################
template<class T,class GENERATOR> template<int d> void RANDOM_NUMBERS<T,GENERATOR>::
Fill_Uniform(SYMMETRIC_MATRIX<T,d>& m,const T a,const T b)
{
    for(int i=0;i<m.Rows();i++) for(int j=0;j<=i;j++) m(i,j)=Get_Uniform_Number(a,b);
}
//#####################################################################
// Function Fill_Uniform_Matrix
//#####################################################################
template<class T,class GENERATOR> template<int d> void RANDOM_NUMBERS<T,GENERATOR>::
Fill_Uniform(UPPER_TRIANGULAR_MATRIX<T,d>& m,const T a,const T b)
{
    for(int j=0;j<d;j++) for(int i=0;i<=j;i++) m(i,j)=Get_Uniform_Number(a,b);
}
//#####################################################################
// Function Fill_Uniform
//#####################################################################
template<class T,class GENERATOR> template<class TV> void RANDOM_NUMBERS<T,GENERATOR>::
Fill_Uniform(TWIST<TV>& m,const T a,const T b)
{
    Fill_Uniform(m.linear,a,b);
    Fill_Uniform(m.angular,a,b);
}
//#####################################################################
// Function Get_Uniform_Vector
//#####################################################################
template<class T,class GENERATOR> template<class TV> TV RANDOM_NUMBERS<T,GENERATOR>::
Get_Uniform_Vector(const RANGE<TV>& box)
{
    return Get_Uniform_Vector(box.min_corner,box.max_corner);
}
//#####################################################################
// Function Get_Gaussian
//#####################################################################
template<class T,class GENERATOR> T RANDOM_NUMBERS<T,GENERATOR>::
Get_Gaussian()
{
    T fac,rsq,v1,v2;
    if(gaussian_iset==0){
        do{v1=2*Get_Uniform_Number((T)0,(T)1)-1;v2=2*Get_Uniform_Number((T)0,(T)1)-1;rsq=sqr(v1)+sqr(v2);}while(rsq>=1 || rsq==0);
        fac=sqrt(-2*log(rsq)/rsq);gset=v1*fac;gaussian_iset=1;return v2*fac;}
    else{gaussian_iset=0;return gset;}
}
//#####################################################################
// Function Get_Vector_In_Unit_Sphere
//#####################################################################
template<class T,class GENERATOR> template<class TV> TV RANDOM_NUMBERS<T,GENERATOR>::
Get_Vector_In_Unit_Sphere()
{
    for(;;){
        TV v=Get_Uniform_Vector(RANGE<TV>::Centered_Box());
        if(v.Magnitude_Squared()<=1) return v;}
}
//#####################################################################
// Function Get_Direction
//#####################################################################
template<class T,class GENERATOR> template<class TV> TV RANDOM_NUMBERS<T,GENERATOR>::
Get_Direction()
{
    if(!TV::m) return TV();
    for(;;){
        TV v=Get_Uniform_Vector(RANGE<TV>::Centered_Box());
        typename TV::SCALAR magnitude_squared=v.Magnitude_Squared();
        if(magnitude_squared>0 && magnitude_squared<=1) return v/sqrt(magnitude_squared);}
}
//#####################################################################
// Function Get_Rotation_Helper
//#####################################################################
template<class T>
ROTATION<VECTOR<T,1> > Get_Rotation_Helper(const VECTOR<T,0>&)
{
    return ROTATION<VECTOR<T,1> >();
}
//#####################################################################
// Function Get_Rotation_Helper
//#####################################################################
template<class T>
ROTATION<VECTOR<T,2> > Get_Rotation_Helper(const VECTOR<T,2>& v)
{
    return ROTATION<VECTOR<T,2> >::From_Complex(COMPLEX<T>(v));
}
//#####################################################################
// Function Get_Rotation_Helper
//#####################################################################
template<class T>
ROTATION<VECTOR<T,3> > Get_Rotation_Helper(const VECTOR<T,4>& v)
{
    return ROTATION<VECTOR<T,3> >::From_Quaternion(QUATERNION<T>(v));
}
//#####################################################################
// Function Get_Rotation
//#####################################################################
template<class T,class GENERATOR> template<class TV> ROTATION<TV> RANDOM_NUMBERS<T,GENERATOR>::
Get_Rotation()
{
    return Get_Rotation_Helper(Get_Direction<VECTOR<T,2*TV::m-2> >());
}
//#####################################################################
// Function Get_Frame
//#####################################################################
template<class T,class GENERATOR> template<class TV> FRAME<TV> RANDOM_NUMBERS<T,GENERATOR>::
Get_Frame(const TV& v0,const TV& v1)
{
    TV v=Get_Uniform_Vector(v0,v1);
    return FRAME<TV>(v,Get_Rotation<TV>());
}
//#####################################################################
// Function Get_Twist
//#####################################################################
template<class T,class GENERATOR> template<class TV> TWIST<TV> RANDOM_NUMBERS<T,GENERATOR>::
Get_Twist(const T& a)
{
    TWIST<TV> tw;
    tw.Set_Vector(Get_Uniform_Vector<T,TWIST<TV>::m>(-a,a));
    return tw;
}
//#####################################################################
#define INSTANTIATION_HELPER_V(T,d) \
    template VECTOR<T,d> RANDOM_NUMBERS<T>::Get_Direction<VECTOR<T,d> >(); \
    template VECTOR<T,d> RANDOM_NUMBERS<T>::Get_Uniform_Vector<VECTOR<T,d> >(RANGE<VECTOR<T,d> > const&); \
    template VECTOR<T,d> RANDOM_NUMBERS<T>::Get_Vector_In_Unit_Sphere<VECTOR<T,d> >(); \
    template ROTATION<VECTOR<T,d> > RANDOM_NUMBERS<T>::Get_Rotation<VECTOR<T,d> >(); \
    template void RANDOM_NUMBERS<T>::Fill_Uniform<VECTOR<T,d> >(TWIST<VECTOR<T,d> >&,T,T);
#define INSTANTIATION_HELPER_V23(T,d) \
    template void RANDOM_NUMBERS<T>::Fill_Uniform<d>(DIAGONAL_MATRIX<T,d>&,T,T); \
    template void RANDOM_NUMBERS<T>::Fill_Uniform<d>(SYMMETRIC_MATRIX<T,d>&,T,T); \
    template void RANDOM_NUMBERS<T>::Fill_Uniform<d>(UPPER_TRIANGULAR_MATRIX<T,d>&,T,T) 
#define IHmn(T,m,n) template void RANDOM_NUMBERS<T>::Fill_Uniform<MATRIX<T,m,n> >(MATRIX_BASE<T,MATRIX<T,m,n> >&,T,T)
#define IHm(T,m) IHmn(T,m,1);IHmn(T,m,2);IHmn(T,m,3);IHmn(T,m,4);IHmn(T,m,5);IHmn(T,m,6)
#define IH(T) IHm(T,1);IHm(T,2);IHm(T,3);IHm(T,4);IHm(T,5);IHm(T,6)
#define INST(T,d) \
    template VECTOR<T,d> RANDOM_NUMBERS<T>::Get_Uniform_Vector<d>(VECTOR<T,d> const&,VECTOR<T,d> const&); \
    template void RANDOM_NUMBERS<T,MT19937<T> >::Fill_Uniform<VECTOR<T,d> >(ARRAY_BASE<T,VECTOR<T,d> >&,T,T);

#define INSTANTIATION_HELPER(T) \
    template class RANDOM_NUMBERS<T>; \
    INSTANTIATION_HELPER_V(T,1); \
    INSTANTIATION_HELPER_V(T,2); \
    INSTANTIATION_HELPER_V(T,3); \
    INSTANTIATION_HELPER_V23(T,0); \
    INSTANTIATION_HELPER_V23(T,1); \
    INSTANTIATION_HELPER_V23(T,2); \
    INSTANTIATION_HELPER_V23(T,3); \
    IH(T);                                                              \
    template void RANDOM_NUMBERS<T>::Fill_Uniform<MATRIX_MXN<T> >(MATRIX_BASE<T,MATRIX_MXN<T> >&,T,T); \
    template void RANDOM_NUMBERS<T>::Fill_Uniform<ARRAY<T> >(ARRAY_BASE<T,ARRAY<T> >&,T,T); \
    INST(T,0);INST(T,1);INST(T,2);INST(T,3); \
    INST(T,4);INST(T,5);INST(T,6);INST(T,7); \
    INST(T,8);INST(T,9);INST(T,10);INST(T,11); \
    INST(T,12);

INSTANTIATION_HELPER(float);
template void RANDOM_NUMBERS<float,MT19937<float> >::Fill_Uniform<ARRAY_VIEW<float> >(ARRAY_BASE<float, ARRAY_VIEW<float> >&, float, float);
INSTANTIATION_HELPER(double);
template void RANDOM_NUMBERS<double,MT19937<double> >::Fill_Uniform<ARRAY_VIEW<double> >(ARRAY_BASE<double, ARRAY_VIEW<double> >&, double, double);
#undef INST
}

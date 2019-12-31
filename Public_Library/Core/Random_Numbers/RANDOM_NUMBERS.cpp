//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Neil Molino, Igor Neverov, Duc Nguyen, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RANDOM_NUMBERS
//#####################################################################
// See Bratley and Fox. 1988. Algorithm 659: Implementing Sobol's quasirandom sequence generator. ACM Trans. Math. Softw. 14, 88-100.
//#####################################################################
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/VECTOR.h>
#include <limits>
using ::std::log;
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> RANDOM_NUMBERS<T>::
RANDOM_NUMBERS(const unsigned int seed)
    :gaussian_iset(0)
{
    Set_Seed(seed);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> RANDOM_NUMBERS<T>::
~RANDOM_NUMBERS()
{}
//#####################################################################
// Function Set_Seed
//#####################################################################
template<class T> void RANDOM_NUMBERS<T>::
Set_Seed(const unsigned int seed_input)
{
    random_number_generator.seed(seed_input);
}
//#####################################################################
// Function Get_Gaussian
//#####################################################################
template<class T> T RANDOM_NUMBERS<T>::
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
template<class T> template<class TV> TV RANDOM_NUMBERS<T>::
Get_Vector_In_Unit_Sphere()
{
    for(;;){
        TV v;
        Fill_Uniform(v,-1,1);
        if(v.Magnitude_Squared()<=1) return v;}
}
//#####################################################################
// Function Get_Direction
//#####################################################################
template<class T> template<class TV> TV RANDOM_NUMBERS<T>::
Get_Direction()
{
    if(TV::m==0) return TV();
    for(;;){
        TV v;
        Fill_Uniform(v,-1,1);
        T x=v.Magnitude_Squared();
        if(x>0 && x<=1) return v/sqrt(x);}
}
//#####################################################################
template class RANDOM_NUMBERS<float>;
template class RANDOM_NUMBERS<double>;
template VECTOR<double,0> RANDOM_NUMBERS<double>::Get_Direction<VECTOR<double,0> >();
template VECTOR<double,1> RANDOM_NUMBERS<double>::Get_Direction<VECTOR<double,1> >();
template VECTOR<double,2> RANDOM_NUMBERS<double>::Get_Direction<VECTOR<double,2> >();
template VECTOR<double,3> RANDOM_NUMBERS<double>::Get_Direction<VECTOR<double,3> >();
template VECTOR<double,4> RANDOM_NUMBERS<double>::Get_Direction<VECTOR<double,4> >();
template VECTOR<float,0> RANDOM_NUMBERS<float>::Get_Direction<VECTOR<float,0> >();
template VECTOR<float,1> RANDOM_NUMBERS<float>::Get_Direction<VECTOR<float,1> >();
template VECTOR<float,2> RANDOM_NUMBERS<float>::Get_Direction<VECTOR<float,2> >();
template VECTOR<float,3> RANDOM_NUMBERS<float>::Get_Direction<VECTOR<float,3> >();
template VECTOR<float,4> RANDOM_NUMBERS<float>::Get_Direction<VECTOR<float,4> >();
template VECTOR<double,2> RANDOM_NUMBERS<double>::Get_Vector_In_Unit_Sphere<VECTOR<double,2> >();
template VECTOR<float,2> RANDOM_NUMBERS<float>::Get_Vector_In_Unit_Sphere<VECTOR<float,2> >();
template VECTOR<double,1> RANDOM_NUMBERS<double>::Get_Vector_In_Unit_Sphere<VECTOR<double,1> >();
template VECTOR<float,1> RANDOM_NUMBERS<float>::Get_Vector_In_Unit_Sphere<VECTOR<float,1> >();
template VECTOR<double,3> RANDOM_NUMBERS<double>::Get_Vector_In_Unit_Sphere<VECTOR<double,3> >();
template VECTOR<float,3> RANDOM_NUMBERS<float>::Get_Vector_In_Unit_Sphere<VECTOR<float,3> >();
template VECTOR<double,0> RANDOM_NUMBERS<double>::Get_Vector_In_Unit_Sphere<VECTOR<double,0> >();
template VECTOR<float,0> RANDOM_NUMBERS<float>::Get_Vector_In_Unit_Sphere<VECTOR<float,0> >();
}

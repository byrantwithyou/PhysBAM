//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Neil Molino, Igor Neverov, Duc Nguyen, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RANDOM_NUMBERS
//#####################################################################
#ifndef __RANDOM_NUMBERS__
#define __RANDOM_NUMBERS__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Read_Write/READ_WRITE_FUNCTIONS.h>
#include <ctime>
#include <random>
namespace PhysBAM{

template<class TV> class RANGE;
template<class TV> class ROTATION;
template<class TV> class FRAME;
template<class T,int d> class VECTOR;
template<class T,int d> class DIAGONAL_MATRIX;
template<class T,int d> class SYMMETRIC_MATRIX;
template<class T,int d> class UPPER_TRIANGULAR_MATRIX;
template<class T,class T_MATRIX> class MATRIX_BASE;
class TYPED_ISTREAM;
class TYPED_OSTREAM;

template<class T>
class RANDOM_NUMBERS
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    int gaussian_iset; // Used to force Get_Gaussian to reset
    T gset; // used internally by Get_Gaussian
    std::mt19937 random_number_generator;

    T Get_Number()
    {return random_number_generator()*(1/(T)4294967296.0);} // in [0,1)

    int Get_Uniform_Integer(const int a,const int b) // in [a,b]
    {return random_number_generator()%(b+1-a)+a;}

    T Get_Uniform_Number(const T a,const T b)
    {return a+(b-a)*Get_Number();} // in [a,b)

    explicit RANDOM_NUMBERS(const unsigned int seed=time(0));
    RANDOM_NUMBERS(const RANDOM_NUMBERS&) = delete;
    void operator=(const RANDOM_NUMBERS&) = delete;
    virtual ~RANDOM_NUMBERS();

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,gaussian_iset,gset);input>>random_number_generator;}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,gaussian_iset,gset);output<<random_number_generator;}

    template<class T2,class T_ARRAY,class ID> void Random_Shuffle(ARRAY_BASE<T2,T_ARRAY,ID>& array)
    {T_ARRAY& derived=array.Derived();for(ID i(derived.Size()-1);i>=ID(0);i--) exchange(array(i),array(Get_Uniform_Integer(0,i)));}

    template<class T_OBJECT,class ...Args> auto
    Fill_Uniform(T_OBJECT& obj,Args&&...args) -> enable_if_t<sizeof((Random_Fill_Uniform(*this,obj,args...),1))>
    {Random_Fill_Uniform(*this,obj,args...);}

    template<class T_OBJECT,class ...Args> auto
    Fill_Uniform(T_OBJECT& obj,T a,T b) -> enable_if_t<sizeof((Random_Fill_Uniform(*this,obj,a,b),1))>
    {Random_Fill_Uniform(*this,obj,a,b);}

//#####################################################################
    void Set_Seed(const unsigned int seed_input=time(0));
    T Get_Gaussian();
    template<class TV> TV Get_Vector_In_Unit_Sphere();
    template<class TV> TV Get_Direction();
//#####################################################################
};
//#####################################################################
// Function Random_Fill_Uniform
//#####################################################################
template<class T> void
Random_Fill_Uniform(RANDOM_NUMBERS<T>& rand,T& x,const T a,const T b)
{
    x=rand.Get_Uniform_Number(a,b);
}
//#####################################################################
// Function Random_Fill_Uniform
//#####################################################################
template<class T,class T2,class T_ARRAY,class ID> void
Random_Fill_Uniform(RANDOM_NUMBERS<T>& rand,ARRAY_BASE<T2,T_ARRAY,ID>& array,const T a,const T b)
{
    for(auto&x:array) rand.Fill_Uniform(x,a,b);
}
}
#endif

//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SMOOTH_GEAR
//##################################################################### 
#ifndef __SMOOTH_GEAR__
#define __SMOOTH_GEAR__

#include <Tools/Math_Tools/constants.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV> class SMOOTH_GEAR;

template<class T>
class SMOOTH_GEAR<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef TV VECTOR_T;

    T r,s,co,ci,den;
    int n;
    TV Co,Ci;

    SMOOTH_GEAR(T R=1,T S=.4,int N=8);
    SMOOTH_GEAR(const TV& dimensions,int N=8);

    struct HELPER
    {
        T a,d,k,b,m,sd;
        bool flip,ui;
        TV Y,dY,P;
    };

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,r,s,n);Compute_Centers();}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,r,s,n);}

//#####################################################################
    void Compute_Centers();
    void Compute_Helper(const TV& X,HELPER& h) const;
    RANGE<TV> Bounding_Box() const;
    T Signed_Distance(const TV& X) const;
    TV Surface(const TV& X) const;
    TV Surface(const TV& X,const HELPER& h) const;
    TV Normal(const TV& X) const;
    TV Normal(const TV& X,const HELPER& h) const;
    TV Normal(const TV& X,const int aggregate) const;
    VECTOR<T,1> Principal_Curvatures(const TV& X) const;
    bool Lazy_Inside(const TV& X) const;
    bool Lazy_Outside(const TV& X) const;
    bool Inside(const TV& X,const T thickness_over_two=0) const;
    bool Outside(const TV& X,const T thickness_over_two=0) const;
    bool Boundary(const TV& X,const T thickness_over_two) const;
    static std::string Name();
//#####################################################################
};   

template<class T>
class SMOOTH_GEAR<VECTOR<T,3> >
{
    typedef VECTOR<T,3> TV;
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef TV VECTOR_T;typedef SMOOTH_GEAR<VECTOR<T,2> > GEAR;

    GEAR g;
    T w;

    SMOOTH_GEAR(T R=1,T S=.4,int N=8,T W=1);
    SMOOTH_GEAR(const TV& dimensions,int N=8);

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,w,g);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,w,g);}

//#####################################################################
    RANGE<TV> Bounding_Box() const;
    T Signed_Distance(const TV& X) const;
    TV Surface(const TV& X) const;
    TV Normal(const TV& X) const;
    TV Normal(const TV& X,const int aggregate) const;
    VECTOR<T,2> Principal_Curvatures(const TV& X) const;
    bool Lazy_Inside(const TV& X) const;
    bool Lazy_Outside(const TV& X) const;
    bool Inside(const TV& X,const T thickness_over_two=0) const;
    bool Outside(const TV& X,const T thickness_over_two=0) const;
    bool Boundary(const TV& X,const T thickness_over_two) const;
    static std::string Name();
//#####################################################################
};
}
#endif

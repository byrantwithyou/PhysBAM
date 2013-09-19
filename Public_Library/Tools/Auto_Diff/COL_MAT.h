//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COL_MAT
//##################################################################### 
#ifndef __COL_MAT__
#define __COL_MAT__

#include <Tools/Auto_Diff/COL_VEC.h>
#include <Tools/Auto_Diff/PRIMITIVE_MATRICES.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{

// 00 = zero, 10 = gen, 01 = id, 11 = scale

template<class TV,class H> struct COL_MAT;

template<int nn,int m,int s> struct CM_FLAGS
{
    enum WA{n=nn,mask=m,sp_mask=s};
    STATIC_ASSERT(mask<(1<<n) && sp_mask<(1<<n));

    CM_FLAGS<n,m,s> operator+() const;
    CM_FLAGS<n,m|s,s> operator-() const;

    template<int m2,int s2> typename H_ADD<CM_FLAGS,CM_FLAGS<n,m2,s2> >::CM operator+(CM_FLAGS<n,m2,s2> x) const;
    template<int m2,int s2> typename H_ADD<CM_FLAGS,CM_FLAGS<n,s2|m2,s2> >::CM operator-(CM_FLAGS<n,m2,s2> x) const;
    CM_FLAGS<n,m|s,s> operator*(double a) const;
    CM_FLAGS<n,m|s,s> operator/(double a) const;
    CV_FLAGS<n,0> operator*(V_FLAGS<0>) const;
    CV_FLAGS<n,m|s> operator*(V_FLAGS<1>) const;
    CM_FLAGS<n,0,0> operator*(M_FLAGS<0,0>) const;
    CM_FLAGS<n,m,s> operator*(M_FLAGS<0,1>) const;
    CM_FLAGS<n,m|s,0> operator*(M_FLAGS<1,0>) const;
    CM_FLAGS<n,m|s,s> operator*(M_FLAGS<1,1>) const;
    CM_FLAGS<n-1,(m>>1),(s>>1)> static Base();
    M_FLAGS<m&1,s&1> static First();
    template<int m2,int s2> CM_FLAGS<n+1,(m<<1)|m2,(s<<1)|s2> static Append(M_FLAGS<m2,s2>);
};
typedef CM_FLAGS<0,0,0> CM_EMPTY;

#define MK_CM(TV,FLAGS) COL_MAT<TV,decltype(FLAGS)>

template<class TV>
struct COL_MAT<TV,CM_EMPTY>
{
    typedef typename TV::SCALAR T;

    MATRIX<T,TV::m> operator()(int i) const {return MATRIX<T,TV::m>();};

    void Neg(const COL_MAT& u) {}

    void Scale(const COL_MAT& u,T a) {}

    void Add(const COL_MAT& u,const COL_MAT& v) {}

    void Sub(const COL_MAT& u,const COL_MAT& v) {}

    template<class T_MAT>
    void Times(const COL_MAT<TV,CM_EMPTY>& u,const T_MAT& v) {}

    void Fill_From(const COL_MAT& u) {}
};

// First index is the long index
template<class TV,class H>
struct COL_MAT:public MK_CM(TV,H::Base())
{
    typedef typename TV::SCALAR T;

    typedef MK_MAT(TV,H::First()) T_DATA;

    T_DATA x;
    typedef MK_CM(TV,H::Base()) BASE;
    const BASE& Base() const {return static_cast<const BASE&>(*this);}
    BASE& Base() {return static_cast<BASE&>(*this);}
    typedef COL_MAT<TV,CM_FLAGS<H::n,0,0> > T_EMPTY;

    MATRIX<T,TV::m> operator()(int i) const {if(i==0) return Cast_Helper(x);return Base()(i-1);}

    template<class H2>
    void Neg(const COL_MAT<TV,H2>& u)
    {x=-u.x;Base().Neg(u.Base());}

    template<class H2>
    void Scale(const COL_MAT<TV,H2>& u,T a)
    {x=u.x*a;Base().Scale(u.Base(),a);}

    template<class H2,class H3>
    void Add(const COL_MAT<TV,H2>& u,const COL_MAT<TV,H3>& v)
    {x=u.x+v.x;Base().Add(u.Base(),v.Base());}

    template<class H2,class H3>
    void Sub(const COL_MAT<TV,H2>& u,const COL_MAT<TV,H3>& v)
    {x=u.x-v.x;Base().Sub(u.Base(),v.Base());}

    template<class H2,class T_MAT>
    void Times(const COL_MAT<TV,H2>& u,const T_MAT& m)
    {x=u.x*m;Base().Times(u.Base(),m);}

    template<class H2>
    void Fill_From(const COL_MAT<TV,H2>& u)
    {Fill_From_Helper(x,u.x);Base().Fill_From(u.Base());}
};

template<class TV,class TV2> void
Outer_Product_Helper(COL_MAT<TV,CM_EMPTY>& r,const COL_VEC<TV,CV_EMPTY>& u,const TV2& v) {}

template<class TV,class H,class G,class TV2> void
Outer_Product_Helper(COL_MAT<TV,H>& r,const COL_VEC<TV,G>& u,const TV2& v)
{
    Outer_Product_Helper(r.x,u.x,v);
    Outer_Product_Helper(r.Base(),u.Base(),v);
}

template<class TV,class G,class H,class TV2>
void Times_Helper(COL_VEC<TV,G>& cv,const COL_MAT<TV,H>& cm,const TV2& v)
{cv.x=cm.x*v;Times_Helper(cv.Base(),cm.Base(),v);}

template<class TV,class MAT>
void Times_Helper(COL_VEC<TV,CV_EMPTY>& cv,const COL_MAT<TV,CM_EMPTY>& cm,const MAT& v) {}

template<class TV,class H2,class H,class MAT>
void Times_Helper(COL_MAT<TV,H2>& cv,const COL_MAT<TV,H>& cm,const MAT& v)
{cv.x=cm.x*v;Times_Helper(cv.Base(),cm.Base(),v);}

template<class TV,class MAT>
void Times_Helper(COL_MAT<TV,CM_EMPTY>& cv,const COL_MAT<TV,CM_EMPTY>& cm,const MAT& v) {}


}
}
#endif

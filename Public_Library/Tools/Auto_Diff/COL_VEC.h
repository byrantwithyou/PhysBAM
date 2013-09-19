//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COL_VEC
//##################################################################### 
#ifndef __COL_VEC__
#define __COL_VEC__

#include <Tools/Auto_Diff/PRIMITIVE_VECTORS.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{

template<class TV,class G> struct COL_VEC;

template<int nn,int m> struct CV_FLAGS
{
    enum WA{n=nn,mask=m};
    STATIC_ASSERT(mask<(1<<n));

    CV_FLAGS<n,m> operator+() const;
    CV_FLAGS<n,m> operator-() const;
    template<int p> CV_FLAGS<n,m|p> operator+(CV_FLAGS<n,p> x) const;
    template<int p> CV_FLAGS<n,m|p> operator-(CV_FLAGS<n,p> x) const;
    CV_FLAGS<n,m> operator*(double a) const;
    CV_FLAGS<n,m> operator/(double a) const;
    CV_FLAGS<n-1,(m>>1)> static Base();
    V_FLAGS<m&1> static First();
    template<int m2> CV_FLAGS<n+1,(m<<1)|m2> static Append(V_FLAGS<m2>);
};
template<int n,int m,int m2> CV_FLAGS<n,m|m2> Choose(CV_FLAGS<n,m>,CV_FLAGS<n,m2>) {return CV_FLAGS<n,m|m2>();}

typedef CV_FLAGS<0,0> CV_EMPTY;
#define F_CV(n,m) CV_FLAGS<n,m>()
#define MK_CV(TV,FLAGS) COL_VEC<TV,decltype(FLAGS)>

template<class TV>
struct COL_VEC<TV,CV_EMPTY>
{
    typedef CV_EMPTY FLAGS;
    typedef typename TV::SCALAR T;

    TV operator()(int i) const {return TV();};

    void Neg(const COL_VEC& u) {}
    void Scale(const COL_VEC& u,T a) {}
    void Add(const COL_VEC& u,const COL_VEC& v) {}
    void Sub(const COL_VEC& u,const COL_VEC& v) {}
    void Fill_From(const COL_VEC& u) {}
};

template<class TV,class G>
struct COL_VEC:public MK_CV(TV,G::Base())
{
    typedef G FLAGS;
    typedef typename TV::SCALAR T;
    typedef MK_VEC(TV,G::First()) T_DATA;
    typedef MK_CV(TV,G::Base()) BASE;
    const BASE& Base() const {return static_cast<const BASE&>(*this);}
    BASE& Base() {return static_cast<BASE&>(*this);}
    T_DATA x;
    typedef COL_VEC<TV,CV_FLAGS<G::n,0> > T_EMPTY;

    TV operator()(int i) const {if(i==0) return Cast_Helper(x);return Base()(i-1);}

    template<int i>
    void Set_Entry(const TV& v,VECTOR<int,i>*) {Base().Set_Entry(v,(VECTOR<int,i-1>*)0);}
    void Set_Entry(const TV& v,VECTOR<int,0>*) {x=v;}

    void Neg(const COL_VEC& u)
    {x=-u.x;Base().Neg(u.Base());}

    void Scale(const COL_VEC& u,T a)
    {x=u.x*a;Base().Scale(u.Base(),a);}

    template<class G2,class G3>
    void Add(const COL_VEC<TV,G2>& u,const COL_VEC<TV,G3>& v)
    {x=u.x+v.x;Base().Add(u.Base(),v.Base());}

    template<class G2,class G3>
    void Sub(const COL_VEC<TV,G2>& u,const COL_VEC<TV,G3>& v)
    {x=u.x-v.x;Base().Sub(u.Base(),v.Base());}

    template<class G2>
    void Fill_From(const COL_VEC<TV,G2>& u)
    {Fill_From_Helper(x,u.x);Base().Fill_From(u.Base());}
};
}
}
#endif

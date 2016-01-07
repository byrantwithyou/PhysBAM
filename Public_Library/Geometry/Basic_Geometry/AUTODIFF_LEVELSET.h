//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class AUTODIFF_LEVELSET
//##################################################################### 
#ifndef __AUTODIFF_LEVELSET__
#define __AUTODIFF_LEVELSET__

#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Tools/Math_Tools/constants.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T,class AH> void inline
Principal_Curvatures_Helper(VECTOR<T,0>& c,const AH& h)
{
    c=VECTOR<T,0>();
}
template<class T,class AH> void inline
Principal_Curvatures_Helper(VECTOR<T,1>& c,const AH& h)
{
    VECTOR<T,2> n;
    SYMMETRIC_MATRIX<T,2> sm;
    Get<0>(n,h.dx);
    Get<0,0>(sm,h.ddx);
    VECTOR<T,2> t=n.Perpendicular();
    c=VECTOR<T,1>(t.Dot(sm*t));
}
template<class T,class AH> void inline
Principal_Curvatures_Helper(VECTOR<T,2>& c,const AH& h)
{
    VECTOR<T,3> n;
    SYMMETRIC_MATRIX<T,3> sm;
    Get<0>(n,h.dx);
    Get<0,0>(sm,h.ddx);
    SYMMETRIC_MATRIX<T,3> P=(T)1-SYMMETRIC_MATRIX<T,3>::Outer_Product(n),M=SYMMETRIC_MATRIX<T,3>::Conjugate(P,sm);
    T trace=M.Trace();
    QUADRATIC<T> quadratic(-1,trace,sqr(M(0,2))+sqr(M(0,1))+sqr(M(1,2))-M(1,1)*M(0,0)-M(2,2)*M(0,0)-M(2,2)*M(1,1));
    quadratic.Compute_Roots();
    if(quadratic.roots == 0) c=(T).5*VECTOR<T,2>(trace,trace);
    else if(quadratic.roots == 1) c=VECTOR<T,2>(quadratic.root1,quadratic.root1);
    else c=VECTOR<T,2>(quadratic.root1,quadratic.root2);
}

template<class TV,class HELPER>
class AUTODIFF_LEVELSET
{
    typedef typename TV::SCALAR T;
public:
    typedef TV VECTOR_T;
    typedef DIFF_LAYOUT<T,TV::m> LAYOUT;

    const HELPER& Derived() const
    {return static_cast<const HELPER&>(*this);}

    T Signed_Distance(const TV& X) const
    {return Derived().Raw_Phi(From_Var<LAYOUT,0>(X)).x;}

    TV Surface(const TV& X) const
    {auto ad=Derived().Raw_Phi(Diff_From_Var<LAYOUT,0>(X));TV n;Get<0>(n,ad.dx);return X-ad.x*n;}

    TV Normal(const TV& X) const
    {TV n;Get<0>(n,Derived().Raw_Phi(Diff_From_Var<LAYOUT,0>(X)).dx);return n;}

    TV Normal(const TV& X,const int aggregate) const
    {return Normal(X);}

    SYMMETRIC_MATRIX<T,TV::m> Hessian(const TV& X) const
    {SYMMETRIC_MATRIX<T,TV::m> H;Get<0,0>(H,Derived().Raw_Phi(Hess_From_Var<LAYOUT,0>(X)).ddx);return H;}

    VECTOR<T,TV::m-1> Principal_Curvatures(const TV& X) const
    {VECTOR<T,TV::m-1> c;Principal_Curvatures_Helper(c,Derived().Raw_Phi(Hess_From_Var<LAYOUT,0>(X)));return c;}

    bool Lazy_Inside(const TV& X) const
    {return Signed_Distance(X)<0;}

    bool Lazy_Outside(const TV& X) const
    {return Signed_Distance(X)>0;}

    bool Inside(const TV& X,const T thickness_over_two=0) const
    {return Signed_Distance(X)<-thickness_over_two;}

    bool Outside(const TV& X,const T thickness_over_two=0) const
    {return Signed_Distance(X)>thickness_over_two;}

    bool Boundary(const TV& X,const T thickness_over_two) const
    {return abs(Signed_Distance(X))<=thickness_over_two;}
};
}
#endif

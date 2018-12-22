//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class AUTODIFF_LEVELSET
//##################################################################### 
#ifndef __AUTODIFF_LEVELSET__
#define __AUTODIFF_LEVELSET__

#include <Core/Math_Tools/constants.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
namespace PhysBAM{

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
    {
        auto h=Derived().Raw_Phi(Hess_From_Var<LAYOUT,0>(X));
        TV n;
        SYMMETRIC_MATRIX<T,TV::m> sm;
        Get<0>(n,h.dx);
        Get<0,0>(sm,h.ddx);
        return Compute_Principal_Curvatures(n,sm);
    }

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

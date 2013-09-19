//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HETERO_VECTOR
//##################################################################### 
#ifndef __HETERO_VECTOR__
#define __HETERO_VECTOR__

#include <Tools/Auto_Diff/COL_VEC.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{

#define MK_G(TV,FLAGS) GRADIENT<TV,decltype(FLAGS)>

template<class TV,class G>
struct GRADIENT
{
    typedef typename TV::SCALAR T;
    COL_VEC<TV,G> c;

    GRADIENT operator- () const {GRADIENT r;r.c.Neg(c);return r;}
    GRADIENT operator+ () const {return *this;}
    GRADIENT operator* (T a) const {GRADIENT r;r.c.Scale(c,a);return r;}
    GRADIENT operator/ (T a) const {return *this*(1/a);}

    template<class G2>
    MK_G(TV,G()+G2()) operator+ (const GRADIENT<TV,G2>& z) const
    {MK_G(TV,G()+G2()) r;r.c.Add(c,z.c);return r;}

    template<class G2>
    MK_G(TV,G()-G2()) operator- (const GRADIENT<TV,G2>& z) const
    {MK_G(TV,G()-G2()) r;r.c.Sub(c,z.c);return r;}

    TV operator()(int i) const {return c(i);}

    template<int i>
    void Set_Entry(const TV& v) {STATIC_ASSERT(G::mask&(1<<i));c.Set_Entry(v,(VECTOR<int,i>*)0);}

    template<class G2>
    void Fill_From(const GRADIENT<TV,G2>& z)
    {STATIC_ASSERT_SAME(G,decltype(Choose(G(),G2())));c.Fill_From(z.c);}
};

template<class TV,class G>
GRADIENT<TV,G> operator* (typename TV::SCALAR a,const GRADIENT<TV,G>& u)
{return u*a;}

typedef VECTOR<double,3> TV;
template struct GRADIENT<TV,CV_FLAGS<2,0> >;
template struct GRADIENT<TV,CV_FLAGS<2,1> >;
template struct GRADIENT<TV,CV_FLAGS<2,2> >;
template struct GRADIENT<TV,CV_FLAGS<2,3> >;
}
}
#endif

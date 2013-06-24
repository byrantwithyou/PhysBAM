//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_LEVELSET_MOVING_ELLIPSOID
//#####################################################################
#ifndef __ANALYTIC_LEVELSET_MOVING_ELLIPSOID__
#define __ANALYTIC_LEVELSET_MOVING_ELLIPSOID__
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Geometry/Analytic_Tests/ANALYTIC_IMPLICIT_SURFACE_LEVELSET.h>
#include <boost/function.hpp>

namespace PhysBAM{
template<class TV>
struct ANALYTIC_LEVELSET_MOVING_ELLIPSOID:public ANALYTIC_IMPLICIT_SURFACE_LEVELSET<TV>
{
    typedef typename TV::SCALAR T;
    TV cen;
    boost::function<TV(T t)> r;
    ANALYTIC_LEVELSET_MOVING_ELLIPSOID(TV cc,boost::function<TV(T t)> r,int c_i,int c_o): ANALYTIC_IMPLICIT_SURFACE_LEVELSET<TV>(c_i,c_o),cen(cc),r(r) {}
    virtual T f(const TV& X,T t) const {return ((X-cen)/r(t)).Magnitude_Squared()-1;}
    virtual TV df(const TV& X,T t) const {return (T)2*(X-cen)/sqr(r(t));}
    virtual MATRIX<T,TV::m> ddf(const TV& X,T t) const {return DIAGONAL_MATRIX<T,2>((T)2/sqr(r(t)));}
    virtual TV Closest_Point_Estimate(const TV& X,T t) const {return X;}
};
}
#endif

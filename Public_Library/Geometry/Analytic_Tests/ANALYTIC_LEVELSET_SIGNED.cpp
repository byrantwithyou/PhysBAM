//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SIGNED.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> ANALYTIC_LEVELSET_SIGNED<TV>::
ANALYTIC_LEVELSET_SIGNED(int color_i,int color_o)
    :color_i(color_i),color_o(color_o)
{
}
//#####################################################################
// Function phi
//#####################################################################
template<class TV> typename TV::SCALAR ANALYTIC_LEVELSET_SIGNED<TV>::
phi(const TV& X,T t,int& c) const
{
    T p=phi2(X,t);
    c=p>0?color_o:color_i;
    return abs(p);
}
//#####################################################################
// Function N
//#####################################################################
template<class TV> TV ANALYTIC_LEVELSET_SIGNED<TV>::
N(const TV& X,T t,int c) const
{
    TV n=N2(X,t);
    return c==color_i?n:-n;
}
//#####################################################################
// Function dist
//#####################################################################
template<class TV> typename TV::SCALAR ANALYTIC_LEVELSET_SIGNED<TV>::
dist(const TV& X,T t,int c) const
{
    if(c!=color_i && c!=color_o)
        return this->Large_Phi();
    T p=phi2(X,t);
    return c==color_i?p:-p;
}
namespace PhysBAM{
template struct ANALYTIC_LEVELSET_SIGNED<VECTOR<float,1> >;
template struct ANALYTIC_LEVELSET_SIGNED<VECTOR<float,2> >;
template struct ANALYTIC_LEVELSET_SIGNED<VECTOR<float,3> >;
template struct ANALYTIC_LEVELSET_SIGNED<VECTOR<double,1> >;
template struct ANALYTIC_LEVELSET_SIGNED<VECTOR<double,2> >;
template struct ANALYTIC_LEVELSET_SIGNED<VECTOR<double,3> >;
}

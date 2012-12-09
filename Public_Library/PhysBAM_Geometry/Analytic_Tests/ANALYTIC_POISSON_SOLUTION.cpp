//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_POISSON_SOLUTION.h>
using namespace PhysBAM;
//#####################################################################
// Function Test
//#####################################################################
template<class TV> void ANALYTIC_POISSON_SOLUTION<TV>::
Test(const TV& X) const
{
    RANDOM_NUMBERS<T> rand;
    TV dX;
    T e=1e-6;
    rand.Fill_Uniform(dX,-e,e);
    T u0=u(X),u1=u(X+dX);
    TV du0=du(X),du1=du(X+dX);
    T AL=0,L=Laplacian(X);
    for(int d=0;d<TV::m;d++){
        TV de=dX(d)*TV::Axis_Vector(d);
        AL+=(du(X+de)-du(X-de))(d)/(2*dX(d));}
    T erru=abs((du0+du1).Dot(dX)/2-(u1-u0))/e;
    LOG::cout<<"analytic poisson solution diff test "<<erru<<"   "<<(L-AL)<<std::endl;
}
namespace PhysBAM{
template class ANALYTIC_POISSON_SOLUTION<VECTOR<float,1> >;
template class ANALYTIC_POISSON_SOLUTION<VECTOR<float,2> >;
template class ANALYTIC_POISSON_SOLUTION<VECTOR<float,3> >;
template class ANALYTIC_POISSON_SOLUTION<VECTOR<double,1> >;
template class ANALYTIC_POISSON_SOLUTION<VECTOR<double,2> >;
template class ANALYTIC_POISSON_SOLUTION<VECTOR<double,3> >;
}

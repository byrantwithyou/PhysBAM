//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/LOG.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET.h>
using namespace PhysBAM;
//#####################################################################
// Destructor
//#####################################################################
template<class TV> ANALYTIC_LEVELSET<TV>::
~ANALYTIC_LEVELSET()
{
}
//#####################################################################
// Function Test
//#####################################################################
template<class TV> void ANALYTIC_LEVELSET<TV>::
Test(const RANGE<TV>& domain) const
{
    RANDOM_NUMBERS<T> rand;
    TV X=rand.Get_Uniform_Vector(domain),dX;
    T e=1e-6,t=rand.Get_Uniform_Number(0,1);
    rand.Fill_Uniform(dX,-e,e);
    int c0,c1;
    T l0=-phi(X,t,c0),l1=-phi(X+dX,t,c1);
    if(c1!=c0) l1=-l1;
    TV dl0=N(X,t,c0),dl1=N(X+dX,t,c0);
    T errl=abs((dl0+dl1).Dot(dX)/2-(l1-l0))/e;
    LOG::cout<<"analytic level set diff test "<<errl<<"   "<<X<<std::endl;
}
namespace PhysBAM{
template class ANALYTIC_LEVELSET<VECTOR<float,1> >;
template class ANALYTIC_LEVELSET<VECTOR<float,2> >;
template class ANALYTIC_LEVELSET<VECTOR<float,3> >;
template class ANALYTIC_LEVELSET<VECTOR<double,1> >;
template class ANALYTIC_LEVELSET<VECTOR<double,2> >;
template class ANALYTIC_LEVELSET<VECTOR<double,3> >;
}

//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include "ANALYTIC_POISSON_TEST.h"
using namespace PhysBAM;
template<class TV> void ANALYTIC_POISSON_TEST<TV>::
Test(const RANGE<TV>& domain) const
{
    RANDOM_NUMBERS<T> rand;
    analytic_levelset->Test(domain);
    int c=-4;
    TV X;
    for(int i=0;i<100 && c<0;i++){
        X=rand.Get_Uniform_Vector(domain);
        analytic_levelset->phi(X,0,c);}
    if(c<0){
        LOG::cout<<"Could not find nonnegative color for poisson solution test."<<std::endl;
        return;}
    for(int i=0;i<analytic_solution.m;i++)
        analytic_solution(i)->Test(X,0);
}
namespace PhysBAM{
template struct ANALYTIC_POISSON_TEST<VECTOR<float,1> >;
template struct ANALYTIC_POISSON_TEST<VECTOR<float,2> >;
template struct ANALYTIC_POISSON_TEST<VECTOR<float,3> >;
template struct ANALYTIC_POISSON_TEST<VECTOR<double,1> >;
template struct ANALYTIC_POISSON_TEST<VECTOR<double,2> >;
template struct ANALYTIC_POISSON_TEST<VECTOR<double,3> >;
}

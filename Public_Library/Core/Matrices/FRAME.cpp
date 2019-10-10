//#####################################################################
// Copyright 2019, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/FRAME.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
namespace PhysBAM{
//#####################################################################
// Function Random_Fill_Uniform
//#####################################################################
template<class T,class TV> void
Random_Fill_Uniform(RANDOM_NUMBERS<T>& rand,FRAME<TV>& f,const T a,const T b)
{
    rand.Fill_Uniform(f.t,a,b);
    rand.Fill_Uniform(f.r);
}
//#####################################################################
// Function Random_Fill_Uniform
//#####################################################################
template<class T,class TV> void
Random_Fill_Uniform(RANDOM_NUMBERS<T>& rand,FRAME<TV>& f,const TV& v0,const TV& v1)
{
    rand.Fill_Uniform(f.t,v0,v1);
    rand.Fill_Uniform(f.r);
}
template void Random_Fill_Uniform<double,VECTOR<double,1> >(
    RANDOM_NUMBERS<double>&,FRAME<VECTOR<double,1> >&,double,double);
template void Random_Fill_Uniform<double,VECTOR<double,2> >(
    RANDOM_NUMBERS<double>&,FRAME<VECTOR<double,2> >&,double,double);
template void Random_Fill_Uniform<double,VECTOR<double,3> >(
    RANDOM_NUMBERS<double>&,FRAME<VECTOR<double,3> >&,double,double);
template void Random_Fill_Uniform<float,VECTOR<float,1> >(
    RANDOM_NUMBERS<float>&,FRAME<VECTOR<float,1> >&,float,float);
template void Random_Fill_Uniform<float,VECTOR<float,2> >(
    RANDOM_NUMBERS<float>&,FRAME<VECTOR<float,2> >&,float,float);
template void Random_Fill_Uniform<float,VECTOR<float,3> >(
    RANDOM_NUMBERS<float>&,FRAME<VECTOR<float,3> >&,float,float);
}

//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/LOG.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_NEST.h>
using namespace PhysBAM;
//#####################################################################
// Destructor
//#####################################################################
template<class TV> ANALYTIC_LEVELSET_NEST<TV>::
~ANALYTIC_LEVELSET_NEST()
{
    delete al;
    sub_al.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function phi
//#####################################################################
template<class TV> typename TV::SCALAR ANALYTIC_LEVELSET_NEST<TV>::
phi(const TV& X,T t,int& c) const
{
    int id=0;
    T p=al->phi(X,t,id);
    return min(p,sub_al(id)->phi(X,t,c));
}
//#####################################################################
// Function N
//#####################################################################
template<class TV> TV ANALYTIC_LEVELSET_NEST<TV>::
N(const TV& X,T t,int c) const
{
    TV N;
    find(X,t,c,&N);
    return N;
}
//#####################################################################
// Function dist
//#####################################################################
template<class TV> typename TV::SCALAR ANALYTIC_LEVELSET_NEST<TV>::
dist(const TV& X,T t,int c) const
{
    return find(X,t,c,0);
}
//#####################################################################
// Function find
//#####################################################################
template<class TV> typename TV::SCALAR ANALYTIC_LEVELSET_NEST<TV>::
find(const TV& X,T t,int c,TV* N) const
{
    T out_dist=this->Large_Phi(),in_dist=-this->Large_Phi(),in_sign=0,out_sign=0;
    int in_color=-1,out_color=-1;
    ANALYTIC_LEVELSET<TV>* out_al=0,*in_al=0;
    bool in=false;
    for(int i=0;i<sub_al.m;i++){
        T d=sub_al(i)->dist(X,t,c),e=al->dist(X,t,i);
        if(e<=0 && d<=0) in=true;
        if(d<=0 && d>=in_dist){in_dist=d;in_al=sub_al(i);in_color=c;in_sign=1;}
        if(e>0 && d>0 && -e>=in_dist){in_dist=-e;in_al=al;in_color=i;in_sign=-1;}
        if(e<=0 && d<=0 && e>=in_dist){in_dist=e;in_al=al;in_color=i;in_sign=1;}
        if(d>0 && d<=out_dist){out_dist=d;out_al=sub_al(i);out_color=c;out_sign=1;}
        if(d<=0 && e>0 && e<=out_dist){out_dist=e;out_al=al;out_color=i;out_sign=1;}}
    if(in){
        if(N) *N=in_sign*in_al->N(X,t,in_color);
        return in_dist;}
    if(N) *N=out_sign*out_al->N(X,t,out_color);
    return out_dist;
}
namespace PhysBAM{
template class ANALYTIC_LEVELSET_NEST<VECTOR<float,1> >;
template class ANALYTIC_LEVELSET_NEST<VECTOR<float,2> >;
template class ANALYTIC_LEVELSET_NEST<VECTOR<float,3> >;
template class ANALYTIC_LEVELSET_NEST<VECTOR<double,1> >;
template class ANALYTIC_LEVELSET_NEST<VECTOR<double,2> >;
template class ANALYTIC_LEVELSET_NEST<VECTOR<double,3> >;
}

//#####################################################################
// Copyright 2012, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
namespace PhysBAM{

struct QUADRATURE_RULE_TRI
{
    bool use_center;
    double center_weight;
    int number3;
    double location3[3];
    double weight3[3];
    int number6;
    double location6a[1],location6b[1];
    double weight6[1];
};

extern const QUADRATURE_RULE_TRI quadrature_rule_tri[9];

struct QUADRATURE_RULE_SEG
{
    bool use_center;
    double center_weight;
    int number;
    double location[3];
    double weight[3];
};

extern const QUADRATURE_RULE_SEG quadrature_rule_seg[12];

// Returns a pair<...>(averaged func, area-weighted normal).
template<class TV,class FUNC> auto Quadrature_Integration(const VECTOR<TV,2>& X, int order, FUNC func)
{
    typedef typename TV::SCALAR T;
    typedef decltype(func(TV())) T2;
    T2 r=T2();
    TV a=X(1)-X(0);
    const QUADRATURE_RULE_SEG& rule=quadrature_rule_seg[order];
    if(rule.use_center) r+=(T)rule.center_weight*func(X.Sum()/2);
    for(int i=0;i<rule.number && i<3;i++){
        T2 s=T2();
        s+=func(X.x+(T)rule.location[i]*a);
        s+=func(X.y-(T)rule.location[i]*a); 
        r+=(T)rule.weight[i]*s;}
    return std::make_pair(r,a.Rotate_Clockwise_90());
}

// Returns a pair<...>(averaged func, area-weighted normal).
template<class TV,class FUNC> auto Quadrature_Integration(const VECTOR<TV,3>& X, int order, FUNC func)
{
    typedef typename TV::SCALAR T;
    typedef decltype(func(TV())) T2;
    T2 r=T2();
    const QUADRATURE_RULE_TRI& rule=quadrature_rule_tri[order];
    if(rule.use_center) r+=(T)rule.center_weight*func(X.Sum()/3);
    for(int i=0;i<rule.number3;i++){
        T2 s=T2();
        s+=func(X.x+(T)rule.location3[i]*(X.y-X.x)+(T)rule.location3[i]*(X.z-X.x));
        s+=func(X.y+(T)rule.location3[i]*(X.z-X.y)+(T)rule.location3[i]*(X.x-X.y));
        s+=func(X.z+(T)rule.location3[i]*(X.x-X.z)+(T)rule.location3[i]*(X.y-X.z));
        r+=(T)rule.weight3[i]*s;}
    for(int i=0;i<rule.number6;i++){
        T2 s=T2();
        s+=func(X.x+(T)rule.location6a[i]*(X.y-X.x)+(T)rule.location6b[i]*(X.z-X.x));
        s+=func(X.y+(T)rule.location6a[i]*(X.z-X.y)+(T)rule.location6b[i]*(X.x-X.y));
        s+=func(X.z+(T)rule.location6a[i]*(X.x-X.z)+(T)rule.location6b[i]*(X.y-X.z));
        s+=func(X.x+(T)rule.location6b[i]*(X.y-X.x)+(T)rule.location6a[i]*(X.z-X.x));
        s+=func(X.y+(T)rule.location6b[i]*(X.z-X.y)+(T)rule.location6a[i]*(X.x-X.y));
        s+=func(X.z+(T)rule.location6b[i]*(X.x-X.z)+(T)rule.location6a[i]*(X.y-X.z));
        r+=(T)rule.weight6[i]*s;}
    return std::make_pair(r,(X(1)-X(0)).Cross(X(2)-X(0))/2);
}
}

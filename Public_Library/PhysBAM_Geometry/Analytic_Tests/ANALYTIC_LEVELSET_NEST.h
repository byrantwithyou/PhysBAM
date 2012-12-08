//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_LEVELSET_NEST
//#####################################################################
#ifndef __ANALYTIC_LEVELSET_NEST__
#define __ANALYTIC_LEVELSET_NEST__
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET.h>

namespace PhysBAM{
template<class TV>
struct ANALYTIC_LEVELSET_NEST:public ANALYTIC_LEVELSET<TV>
{
    typedef typename TV::SCALAR T;
    ANALYTIC_LEVELSET<TV>* al;
    ARRAY<ANALYTIC_LEVELSET<TV>*> sub_al;

    ANALYTIC_LEVELSET_NEST(ANALYTIC_LEVELSET<TV>* ls): al(ls) {}
    ANALYTIC_LEVELSET_NEST* Add(ANALYTIC_LEVELSET<TV>* ls){sub_al.Append(ls);return this;}
    ~ANALYTIC_LEVELSET_NEST() {delete al;sub_al.Delete_Pointers_And_Clean_Memory();}
    virtual T phi(const TV& X,T t,int& c) const {int id=0;T p=al->phi(X,t,id);return min(p,sub_al(id)->phi(X,t,c));}
    virtual TV N(const TV& X,T t,int c) const {ANALYTIC_LEVELSET<TV>* ls=0;find(X,t,c,ls);return ls->N(X,t,c);}
    virtual T dist(const TV& X,T t,int c) const {ANALYTIC_LEVELSET<TV>* ls=0;return find(X,t,c,ls);}
    T find(const TV& X,T t,int c,ANALYTIC_LEVELSET<TV>*& ls) const
    {
        T out_dist=this->Large_Phi(),in_dist=-this->Large_Phi();
        ANALYTIC_LEVELSET<TV>* out_al=0,*in_al=0;
        bool in=false;
        for(int i=0;i<sub_al.m;i++){
            T d=sub_al(i)->dist(X,t,c),e=al->dist(X,t,i);
            if(e<=0 && d<=0) in=true;
            if(d<=0 && d>=in_dist){in_dist=d;in_al=sub_al(i);}
            if(e>0 && d>0 && -e>=in_dist){in_dist=-e;in_al=al;}
            if(e<=0 && d<=0 && e>=in_dist){in_dist=e;in_al=al;}
            if(d>0 && d<=out_dist){out_dist=d;out_al=sub_al(i);}
            if(d<=0 && e>0 && e<=out_dist){out_dist=e;out_al=al;}}
        if(in){ls=in_al;return in_dist;}
        ls=out_al;
        return out_dist;
    }
};
}
#endif

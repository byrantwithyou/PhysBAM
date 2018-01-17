//#####################################################################
// Copyright 2018.  Ounan Ding, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIFF_TEST
//#####################################################################
#ifndef __DIFF_TEST__
#define __DIFF_TEST__

#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Log/LOG.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <cmath>
#include <string>
namespace PhysBAM{
using ::std::abs;

template<class A,class B,class NEST> 
struct DIFF_TEST_HOLDER
{
    A var;
    B diff;
    NEST nest;
};

template<class T,class U>
auto Compute_Total_Diff_Mult(const T& a,const U& b){return a*b;}

template<class T,int m>
T Compute_Total_Diff_Mult(const VECTOR<T,m>& a,const VECTOR<T,m>& b){return a.Dot(b);}

template<class T,int m>
T Compute_Total_Diff_Magnitude(const VECTOR<T,m>& a){return a.Magnitude();}

float Compute_Total_Diff_Magnitude(const float& a){return std::abs(a);}
double Compute_Total_Diff_Magnitude(const double& a){return std::abs(a);}

int Make_Holder(){return 0;}

template<class A,class B,class ...Args>
auto Make_Holder(const A& a,const B& b,Args&&... args)
{
    auto nest=Make_Holder(args...);
    DIFF_TEST_HOLDER<A,B,decltype(nest)> h={a,b,nest};
    return h;
}

template<class T>
T Compute_Total_Diff(int,int)
{
    return T();
}

template<class T,class PACK>
T Compute_Total_Diff(const PACK& a0,const PACK& a1)
{
    return Compute_Total_Diff_Mult(a0.diff+a1.diff,a1.var-a0.var)/2+Compute_Total_Diff<T>(a0.nest,a1.nest);
}

int Compute_Total_Diff_Eps(int,int){return 0;}

template<class PACK>
auto Compute_Total_Diff_Eps(const PACK& a0,const PACK& a1)
{
    auto mine=Compute_Total_Diff_Magnitude(a1.var-a0.var);
    auto rec=Compute_Total_Diff_Eps(a0.nest,a1.nest);
    typedef decltype(mine+rec) TYPE;
    return std::max((TYPE)mine,(TYPE)rec);
}

template<class T,class ...Args> 
void Diff_Test(const char* name,int step,const T& Y,Args&&... args)
{
    if(step!=0 && step!=1) return;
    typedef decltype(Make_Holder(args...)) PACK;
    static HASHTABLE<std::string,PAIR<T,PACK> > hash;

    PAIR<T,PACK> cur={Y,Make_Holder(args...)};
    if(step==0) hash.Set(name,cur);
    else if(step==1)
    {
        auto& old=hash.Get(name);
        T dYa=cur.x-old.x;
        T dYb=Compute_Total_Diff<T>(old.y,cur.y);
        auto eps=Compute_Total_Diff_Eps(old.y,cur.y);
        LOG::printf("diff test %s: %.16P %.16P %P (delta %.16P)\n",name,dYa/eps,dYb/eps,Compute_Total_Diff_Magnitude(dYa-dYb)/eps,eps);
        hash.Delete(name);
    }
}
}
#endif

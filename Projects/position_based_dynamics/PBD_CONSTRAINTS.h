//#####################################################################
// Copyright 2015
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PBD_CONSTRAINTS__
#define __PBD_CONSTRAINTS__
#include <Core/Log/LOG.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV>
class PBD_CONSTRAINTS_BASE
{
    typedef typename TV::SCALAR T;
public:
    PBD_CONSTRAINTS_BASE();
    PBD_CONSTRAINTS_BASE(const PBD_CONSTRAINTS_BASE&) = delete;
    void operator=(const PBD_CONSTRAINTS_BASE&) = delete;
    virtual ~PBD_CONSTRAINTS_BASE();

    virtual void Project(ARRAY<TV>& P,const ARRAY<T>& w,int solver_iterations) const=0;
    virtual void Test_Diff(const ARRAY<TV>& P) const=0;
};

template<class TV,class DATA,int n,class FUNC>
class PBD_CONSTRAINTS:public PBD_CONSTRAINTS_BASE<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,n> IV;
public:
    ARRAY<IV> indices;
    ARRAY<T> stiffness;
    ARRAY<DATA> data;
    const FUNC compute_constraint;

    PBD_CONSTRAINTS(const FUNC& compute_constraint): compute_constraint(compute_constraint) {}

    ~PBD_CONSTRAINTS() {}

    void Project(ARRAY<TV>& P,const ARRAY<T>& w,int solver_iterations) const
    {
        for(int j=0;j<indices.m;j++){
            IV index=indices(j);
            T k=stiffness(j);
            if(k!=1) k=1-pow(1-k,solver_iterations);
            VECTOR<TV,n> grad_c;
            T c=compute_constraint(VECTOR<TV,n>(P.Subset(index)),data(j),grad_c);
            T s=0;
            for(int i=0;i<n;i++) s+=w(index(i))*grad_c(i).Magnitude();
            s=c/s;
            for(int i=0;i<n;i++) P(index(i))-=grad_c(i)*(w(index(i))*k*s);}
    }

    void Test_Diff(const ARRAY<TV>& P) const
    {
        T eps=1e-6;
        VECTOR<TV,n> grad_c0;
        VECTOR<TV,n> grad_c1;
        VECTOR<TV,n> dP;
        static RANDOM_NUMBERS<T> rnd;

        for(int j=0;j<indices.m;j++){
            IV index=indices(j);
            rnd.Fill_Uniform(dP,-eps,eps);
            T c0=compute_constraint(VECTOR<TV,n>(P.Subset(index)),data(j),grad_c0);
            T c1=compute_constraint(VECTOR<TV,n>(P.Subset(index))+dP,data(j),grad_c1);
            T x=c1-c0;
            T y=0;
            for(int i=0;i<n;i++)
                y+=(grad_c0(i)+grad_c1(i)).Dot(dP(i))/2;
            LOG::printf("DIFF %P: %.16g %.16g -> %.16g\n",index,x/eps,y/eps,(x-y)/max(max(abs(x),abs(y)),(T)1e-30));}
    }

    void Add_Constraint(IV index,DATA d=DATA(),T k=1)
    {
        indices.Append(index);
        data.Append(d);
        stiffness.Append(k);
    }
};
}
#endif

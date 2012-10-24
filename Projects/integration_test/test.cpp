//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

using namespace PhysBAM;

const int d=2;
typedef double T;
typedef VECTOR<T,d> TV;
typedef VECTOR<int,d> TV_INT;

T Compute(const ARRAY<TV>& x,const TV& n,TV& dE,SYMMETRIC_MATRIX<T,d>& ddE)
{
    TV average=x.Average();
    
    const SYMMETRIC_MATRIX<T,TV::m> nnt=SYMMETRIC_MATRIX<T,TV::m>::Outer_Product(n);
    const SYMMETRIC_MATRIX<T,TV::m> dn=(T)1-nnt;    

    T E=0;
        
    const T average_dot_n=average.Dot(n);
    for(int i=0;i<x.m;i++){
        const T xi_dot_n=x(i).Dot(n);
        dE+=(x(i)-average)*(2*xi_dot_n);
        E+=sqr(xi_dot_n);}
    dE=dn*dE;
    E-=sqr(average_dot_n)*x.m;
    
    for(int i=0;i<x.m;i++){
        const TV xma=x(i)-average;
        SYMMETRIC_MATRIX<T,TV::m> tmp=
            MATRIX<T,TV::m>::Outer_Product(n,xma).Symmetric_Part()*(-2)
            +(nnt*(T)3-1)*n.Dot(xma);
        ddE+=tmp*(x(i).Dot(n));
        ddE+=SYMMETRIC_MATRIX<T,TV::m>::Outer_Product(dn*x(i));}
    ddE-=SYMMETRIC_MATRIX<T,TV::m>::Outer_Product(dn*average)*x.m;
    ddE*=2;
    
    return E;
}

int main(int argc, char* argv[])
{
    for(int e=0;e<12;e++){
        const int axis=e>>2;
        const int offset1=e&1;
        const int offset2=(e>>1)&1;
        LOG::cout<<e<<"\t"<<axis<<"\t"<<offset1<<"\t"<<offset2<<"\t";
        LOG::cout<<VECTOR<int,2>(offset1,offset2).Insert(0,axis)<<std::endl;}


    // const int x_full=2;
    // RANDOM_NUMBERS<T> rand;

    // ARRAY<TV> x1(x_full);
    // TV n1;
    // TV dE1;
    // SYMMETRIC_MATRIX<T,d> ddE1;

    // x1(0)=TV(.2,.5);
    // x1(1)=TV(.7,.1);
    // n1=TV(1,1); n1.Normalize();

    // ARRAY<TV> x2(x1);
    // TV n2(n1);
    // TV dE2;
    // SYMMETRIC_MATRIX<T,d> ddE2;

    // const T e=1e-6;
    // TV dn;
    // rand.Fill_Uniform(dn,-e,e);
    // n2+=dn;//n2.Normalize();
    // dn=n2-n1;

    // T E1=Compute(x1,n1,dE1,ddE1);
    // T E2=Compute(x2,n2,dE2,ddE2);

    // LOG::cout<<"dE  "<<(E2-E1-(dE1.Dot(dn)+dE2.Dot(dn))*(T).5)/e<<std::endl;
    // LOG::cout<<"ddE "<<(dE2-dE1-(ddE1*dn+ddE2*dn)*(T).5).Max_Abs()/e<<std::endl;

    // LOG::cout<<E1<<" "<<dE1<<" "<<ddE1<<std::endl;
    // LOG::cout<<E2<<" "<<dE2<<" "<<ddE2<<std::endl;

    return 0;
}

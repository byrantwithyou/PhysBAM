//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Tools/Log/LOG.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Symbolics/STATIC_POLYNOMIAL.h>
#include <Tools/Vectors/VECTOR.h>

using namespace PhysBAM;

typedef double T;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,3> TV_INT;
typedef VECTOR<T,2> TV2;
typedef VECTOR<int,2> TV_INT2;

template<class T,int rank,int d>
void Integral_Test()
{
    typedef VECTOR<T,rank> TV;
    RANDOM_NUMBERS<T> random;
    VECTOR<TV,3> pts;
    random.Fill_Uniform(pts(0),-(T)1,(T)1);
    random.Fill_Uniform(pts(1),-(T)1,(T)1);
    random.Fill_Uniform(pts(2),-(T)1,(T)1);

    STATIC_POLYNOMIAL<T,rank,d> P;
    for(int i=0;i<10;i++){
        TV power;
        random.Fill_Uniform(power,0,(T)5);
        P.Set_Term(VECTOR<int,TV::m>(power),random.Get_Uniform_Number(-1,1));}
    std::cout<<P<<std::endl;

    int N=1000000;
    T i0=P.Integrate_Over_Primitive(pts),i1=0;
    for(int i=0;i<N;i++){
        VECTOR<T,2> uv;
        random.Fill_Uniform(uv,0,1);
        if(uv.Sum()>1) uv=(T)1-uv;
        TV pt=pts(0)+(pts(1)-pts(0))*uv.x+(pts(2)-pts(0))*uv.y;
        i1+=P.Value(pt);}
    i1*=TV::Cross_Product(pts(1)-pts(0),pts(2)-pts(0)).Magnitude()/(2*N);
    LOG::cout<<"integrals c "<<i0<<"  r "<<i1<<"  e "<<abs(i1-i0)<<std::endl;
}

template<class T,int rank,int d>
void Integral_Test2()
{
    typedef VECTOR<T,rank> TV;
    RANDOM_NUMBERS<T> random;
    VECTOR<TV,2> pts;
    random.Fill_Uniform(pts(0),-(T)1,(T)1);
    random.Fill_Uniform(pts(1),-(T)1,(T)1);

    STATIC_POLYNOMIAL<T,rank,d> P;
    for(int i=0;i<10;i++){
        TV power;
        random.Fill_Uniform(power,0,(T)5);
        P.Set_Term(VECTOR<int,TV::m>(power),random.Get_Uniform_Number(-1,1));}
    std::cout<<P<<"     ";

    int N=1000000;
    T i0=P.Integrate_Over_Primitive(pts),i1=0;
    for(int i=0;i<N;i++){
        T u;
        random.Fill_Uniform(u,0,1);
        TV pt=pts(0)+(pts(1)-pts(0))*u;
        i1+=P.Value(pt);}
    i1*=(pts(1)-pts(0)).Magnitude()/N;
    LOG::cout<<"integrals c "<<i0<<"  r "<<i1<<"  e "<<abs(i1-i0)<<std::endl;
}

template<class T,int rank,int d>
void Test_Basic()
{
    typedef VECTOR<int,rank> TV_INT;
    STATIC_POLYNOMIAL<T,rank,d> p[rank],c,q;
    ARRAY<STATIC_POLYNOMIAL<T,rank,d> > poly_list;

    c.Set_Term(TV_INT(),3);
    poly_list.Append(c);
    for(int i=0;i<rank;i++){
        p[i].Set_Term(TV_INT::Axis_Vector(i),2);
        poly_list.Append(p[i]);}
    for(int i=0;i<rank;i++){
        p[i].Set_Term(TV_INT::Axis_Vector(i),-5);
        poly_list.Append(p[i]);}
    for(int i=0;i<rank;i++){
        p[i].Set_Term(TV_INT::Axis_Vector(i),-3);
        poly_list.Append(p[i]);}

    for(int i=0;i<40;i++){
        int r=rand()%poly_list.m;
        int s=rand()%poly_list.m;
        poly_list.Append(poly_list(r)+poly_list(s));
        poly_list.Last().Compress_Size();
        LOG::cout<<poly_list(r)<<"  +  "<<poly_list(s)<<"  =  "<<poly_list.Last()<<"  "<<poly_list.Last().size<<"  "<<poly_list(r).size<<"  +  "<<poly_list(s).size<<std::endl;
        poly_list.Append(poly_list(r)-poly_list(s));
        poly_list.Last().Compress_Size();
        LOG::cout<<poly_list(r)<<"  -  "<<poly_list(s)<<"  =  "<<poly_list.Last()<<"  "<<poly_list.Last().size<<"  "<<poly_list(r).size<<"  -  "<<poly_list(s).size<<std::endl;
        c.Multiply(poly_list(r),poly_list(s));
        poly_list.Append(c);
        poly_list.Last().Compress_Size();
        LOG::cout<<poly_list(r)<<"  *  "<<poly_list(s)<<"  =  "<<poly_list.Last()<<"  "<<poly_list.Last().size<<"  "<<poly_list(r).size<<"  *  "<<poly_list(s).size<<std::endl;
        for(int i=0;i<rank;i++){
            poly_list.Append(c.Differentiate(i));
            LOG::cout<<c<<"  "<<c.size<<"  diff("<<i<<")  ->  "<<poly_list.Last()<<std::endl;}
        for(int i=0;i<rank;i++){
            poly_list.Append(q=c.Integrate(i));
            LOG::cout<<c<<"  int("<<i<<")  ->  "<<poly_list.Last()<<std::endl;}
        c.Exchange_Variables(0,2);
        LOG::cout<<c<<"  exchange"<<std::endl;
}
}

int main(int argc, char* argv[])
{
//    Test_Basic<double,3,10>();

//    Integral_Test<T,3,5>();
    Integral_Test2<T,2,5>();
//    Integral_Test<VECTOR<T,3> >();
    
    return 0;
}

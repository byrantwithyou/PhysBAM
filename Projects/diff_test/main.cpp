//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#include <Tools/Tensors/VEC_ID_TENSOR.h>
#include <Tools/Vectors/VECTOR.h>

using namespace PhysBAM;

typedef float RW;
typedef double T;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,TV::m> TV_INT;
RANDOM_NUMBERS<T> rnd;
static const int N=2;

T eps=1e-6;

struct HELPER
{
    SYMMETRIC_MATRIX<T,TV::m> A;
    MATRIX<T,TV::m> db[N],ddb[N];
    TV u,v,x[N],b,dc[N],ddc[N];
    T c;

    VEC_ID_SYM_TENSOR<T,2,TV::m> dA[N],ddA[N];
};

void Init(HELPER& z)
{
    for(int i=0;i<N;i++) rnd.Fill_Uniform(z.x[i],-1,1);

    rnd.Fill_Uniform(z.u,-1,1);
    rnd.Fill_Uniform(z.v,-1,1);
}

void Fill(HELPER& z)
{
    for(int i=0;i<N;i++) z.A+=SYMMETRIC_MATRIX<T,TV::m>::Outer_Product(z.x[i]);
    z.b=z.A*z.u;
    z.c=z.b.Dot(z.v);
}

void Fill_Diff(HELPER& z)
{
    IDENTITY_MATRIX<T,TV::m> id;

    for(int i=0;i<N;i++){
        z.dA[i].v=z.x[i];
        z.db[i]=Contract<1>(z.dA[i],z.u);
        z.dc[i]=z.db[i].Transpose_Times(z.v);
    }
}

void Fill_Hess(HELPER& z,TV dx[N])
{
    IDENTITY_MATRIX<T,TV::m> id;

    for(int i=0;i<N;i++){
        z.ddA[i].v=dx[i];
        z.ddb[i]=Contract<1>(z.ddA[i],z.u);
        z.ddc[i]=z.ddb[i].Transpose_Times(z.v);
    }
}

void Test(const char* name,TV dx[N],T a,T b,TV da[N],TV db[N])
{
    T x=b-a;
    T y=(db[0]+da[0]).Dot(dx[0])/2+(db[1]+da[1]).Dot(dx[1])/2;
    LOG::printf("DIFF %s: %.16g %.16g -> %.16g\n",name,x/eps,y/eps,(x-y)/max(max(abs(x),abs(y)),1e-30));
}

void Test(const char* name,TV dx[N],TV a,TV b,MATRIX<T,TV::m> da[N],MATRIX<T,TV::m> db[N])
{
    TV x=b-a;
    TV y=(db[0]+da[0])*dx[0]/2+(db[1]+da[1])*dx[1]/2;
    LOG::printf("DIFF %s: %.16g %.16g -> %.16g\n",name,x.Magnitude()/eps,y.Magnitude()/eps,(x-y).Magnitude()/max(max(x.Magnitude(),y.Magnitude()),1e-30));
}

template<class T_TENSOR>
void Test(const char* name,TV dx[N],MATRIX<T,TV::m> a,MATRIX<T,TV::m> b,T_TENSOR da[N],T_TENSOR db[N])
{
    MATRIX<T,TV::m> x=b-a;
    MATRIX<T,TV::m> y=Contract<N>(db[0]+da[0],dx[0])/2+Contract<N>(db[1]+da[1],dx[1])/2;
    LOG::printf("DIFF %s: %.16g %.16g -> %.16g\n",name,x.Frobenius_Norm()/eps,y.Frobenius_Norm()/eps,(x-y).Frobenius_Norm()/max(max(x.Frobenius_Norm(),y.Frobenius_Norm()),1e-30));
}

void Test(const char* name,TV dx[N],TV a[N],TV b[N],MATRIX<T,TV::m> da[N][N],MATRIX<T,TV::m> db[N][N])
{
    VECTOR<TV,2> x,y;
    for(int i=0;i<N;i++){
        x(i)=b[i]-a[i];
        y(i)=(db[i][0]+da[i][0])*dx[0]/2+(db[i][1]+da[i][1])*dx[1]/2;}
    LOG::printf("DIFF %s: %.16g %.16g -> %.16g\n",name,x.Flattened().Magnitude()/eps,y.Flattened().Magnitude()/eps,(x-y).Flattened().Magnitude()/max(max(x.Flattened().Magnitude(),y.Flattened().Magnitude()),1e-30));
}

void Test(const char* name,TV dx[N],TV a[N],TV b[N],TV da[N],TV db[N])
{
    VECTOR<TV,2> x,y;
    for(int i=0;i<N;i++){
        x(i)=b[i]-a[i];
        y(i)=(db[i]+da[i])/2;}
    LOG::printf("DIFF %s: %.16g %.16g -> %.16g\n",name,x.Flattened().Magnitude()/eps,y.Flattened().Magnitude()/eps,(x-y).Flattened().Magnitude()/max(max(x.Flattened().Magnitude(),y.Flattened().Magnitude()),1e-30));
}

template<class T_MATRIX0,class T_MATRIX1> typename enable_if<IS_MATRIX<T_MATRIX0>::value && IS_MATRIX<T_MATRIX1>::value>::type
Test(const char* name,TV dx[N],T_MATRIX0 a[N],T_MATRIX0 b[N],T_MATRIX1 da[N],T_MATRIX1 db[N])
{
    MATRIX<T,TV::m> x[N],y[N],z[N];
    for(int i=0;i<N;i++){
        x[i]=b[i]-a[i];
        y[i]=(db[i]+da[i])/2;
        z[i]=x[i]-y[i];}
    T Mx=sqrt(x[0].Frobenius_Norm_Squared()+x[1].Frobenius_Norm_Squared());
    T My=sqrt(y[0].Frobenius_Norm_Squared()+y[1].Frobenius_Norm_Squared());
    T Mz=sqrt(z[0].Frobenius_Norm_Squared()+z[1].Frobenius_Norm_Squared());
    LOG::printf("DIFF %s: %.16g %.16g -> %.16g\n",name,Mx/eps,My/eps,Mz/max(max(Mx,My),1e-30));
}

template<class T_TENSOR0,class T_TENSOR1> typename enable_if<IS_TENSOR<T_TENSOR0>::value && IS_TENSOR<T_TENSOR1>::value>::type
Test(const char* name,TV dx[N],T_TENSOR0 a[N],T_TENSOR0 b[N],T_TENSOR1 da[N],T_TENSOR1 db[N])
{
    TENSOR<T,TV::m> x[N],y[N],z[N];
    for(int i=0;i<N;i++){
        x[i]+=b[i]-a[i];
        y[i]+=(db[i]+da[i])/2;
        z[i]=x[i]-y[i];}
    T Mx=0,My=0,Mz=0;
    for(int i=0;i<N;i++)
        for(int j=0;j<TV::m;j++){
            Mx+=x[i].x[j].Frobenius_Norm_Squared();
            My+=y[i].x[j].Frobenius_Norm_Squared();
            Mz+=z[i].x[j].Frobenius_Norm_Squared();}
    Mx=sqrt(Mx);
    My=sqrt(My);
    Mz=sqrt(Mz);
    LOG::printf("DIFF %s: %.16g %.16g -> %.16g\n",name,Mx/eps,My/eps,Mz/max(max(Mx,My),1e-30));
}

#define TEST(x) Test(#x,dx,z0.x,z1.x,z0.d##x,z1.d##x);

int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Parse();

    HELPER z0;
    Init(z0);

    TV dx[N];
    for(int i=0;i<N;i++) rnd.Fill_Uniform(dx[i],-eps,eps);
    HELPER z1=z0;
    for(int i=0;i<N;i++) z1.x[i]+=dx[i];
        
    Fill(z0);
    Fill(z1);
        
    Fill_Diff(z0);
    Fill_Diff(z1);
        
    Fill_Hess(z0,dx);
    Fill_Hess(z1,dx);

    TEST(A);
    TEST(b);
    TEST(c);
    TEST(dA);
    TEST(db);
    TEST(dc);

    return 0;
}

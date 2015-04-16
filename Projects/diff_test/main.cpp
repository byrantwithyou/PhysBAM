//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#include <Tools/Tensors/VEC_ID_TENSOR_1.h>
#include <Tools/Vectors/VECTOR.h>

using namespace PhysBAM;

typedef float RW;
typedef double T;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,TV::m> TV_INT;
RANDOM_NUMBERS<T> rnd;

T eps=1e-6;

struct HELPER
{
    MATRIX<T,TV::m> SSn,Fn;
    TV dw[2],xn[2],xh[2];
    T dt,W,mu,la,mu_N,V0;

    MATRIX<T,TV::m> A,B,F_hat,G,S,H,K;
    T J,a,b,q,c,g,h,k,m,n,p,Phi;

    VEC_ID_TENSOR_1<T,TV::m,TV::m> dA[2],dB[2],dF_hat[2];
    TENSOR<T,TV::m,TV::m,TV::m> dS[2],dH[2],dK[2];
    TV B_bar[2],F_bar[2],K_bar[2],H_bar[2],dJ[2],da[2],db[2],dc[2],dg[2],dq[2],dh[2],dk[2],dm[2],dn[2],dp[2],dPhi[2];

    MATRIX<T,TV::m> ddJ[2][2],dda[2][2],ddb[2][2],ddc[2][2],ddq[2][2],ddk[2][2],ddm[2][2],ddn[2][2],ddp[2][2],ddPhi[2][2];
};

void Init(HELPER& z)
{
    rnd.Fill_Uniform(z.SSn,-1,1);
    z.SSn=z.SSn.Times_Transpose(z.SSn);

    rnd.Fill_Uniform(z.Fn,-1,1);
    z.Fn=z.Fn+3;

    for(int i=0;i<2;i++){
        rnd.Fill_Uniform(z.dw[i],-1,1);
        rnd.Fill_Uniform(z.xn[i],-1,1);
        rnd.Fill_Uniform(z.xh[i],-1,1);}

    rnd.Fill_Uniform(z.dt,.001,1);
    rnd.Fill_Uniform(z.W,.001,1);
    rnd.Fill_Uniform(z.mu,.001,1);
    rnd.Fill_Uniform(z.la,.001,1);
    rnd.Fill_Uniform(z.mu_N,.001,1);
    rnd.Fill_Uniform(z.V0,.001,1);
}

void Fill(HELPER& z)
{
    for(int i=0;i<2;i++) z.A+=MATRIX<T,TV::m>::Outer_Product(z.xh[i]-z.xn[i],z.dw[i]);
    z.B=z.A*z.SSn;
    z.F_hat=(z.A+1)*z.Fn;
    z.G=z.SSn+z.dt/z.W*((T)1-z.SSn);
    z.S=z.G+z.B.Twice_Symmetric_Part();
    z.H=z.F_hat.Inverse();
    z.J=z.F_hat.Determinant();
    z.a=z.la/2*sqr(z.J-1);
    z.b=z.mu*log(z.J);
    z.q=1/(2*sqr(z.dt))*(z.A.Frobenius_Norm_Squared()+z.A.Double_Contract(z.A.Transposed()));
    z.c=z.mu_N*z.Fn.Determinant()*z.q;
    z.g=z.S.Trace();
    z.K=z.S.Inverse();
    z.h=z.S.Determinant();
    z.k=pow(z.h,-1.0/TV::m);//
    z.m=pow(z.J,2.0/TV::m);
    z.n=z.k*z.g;
    z.p=z.mu/2*z.m*z.n;
    z.Phi=z.V0*(z.p-z.b+z.a+z.c);
}

void Fill_Diff(HELPER& z)
{
    IDENTITY_MATRIX<T,TV::m> id;

    for(int i=0;i<2;i++){
        z.dA[i].v=z.dw[i];
        z.B_bar[i]=z.SSn*z.dw[i];
        z.dB[i].v=z.B_bar[i];
        z.F_bar[i]=z.Fn.Transpose_Times(z.dw[i]);
        z.dF_hat[i].v=z.F_bar[i];
        z.dS[i]=Tensor_Product_1(id,z.B_bar[i])+Tensor_Product_0(id,z.B_bar[i]);
        z.H_bar[i]=z.H.Transpose_Times(z.F_bar[i]);
        z.dH[i]=Tensor_Product_1(z.H,-z.H_bar[i]);
        z.dJ[i]=z.J*z.H_bar[i];
        z.da[i]=z.la*(z.J-1)*z.dJ[i];
        z.db[i]=z.mu*z.H_bar[i];
        z.dq[i]=(z.A*z.dw[i]+z.A.Transpose_Times(z.dw[i]))/sqr(z.dt);
        z.dc[i]=z.mu_N*z.Fn.Determinant()*z.dq[i];
        z.dg[i]=z.B_bar[i]*2;
        z.K_bar[i]=z.K.Transpose_Times(z.B_bar[i]);
        z.dK[i]=-Tensor_Product_1(z.K,z.K_bar[i])-Tensor_Product_0(z.K.Transposed(),z.K_bar[i]);
        z.dh[i]=2*z.h*z.K_bar[i];
        z.dk[i]=-2*z.k/TV::m*z.K_bar[i];
        z.dm[i]=2*z.m/TV::m*z.H_bar[i];
        z.dn[i]=z.dk[i]*z.g+z.k*z.dg[i];
        z.dp[i]=z.mu/2*z.dm[i]*z.n+z.mu/2*z.m*z.dn[i];
        z.dPhi[i]=z.V0*(z.dp[i]-z.db[i]+z.da[i]+z.dc[i]);}
}

void Fill_Hess(HELPER& z)
{
    for(int j=0;j<2;j++){
        for(int k=0;k<2;k++){
            z.ddJ[j][k]=MATRIX<T,TV::m>::Outer_Product(z.H_bar[j],z.J*z.H_bar[k])-MATRIX<T,TV::m>::Outer_Product(z.J*z.H_bar[k],z.H_bar[j]);
            z.dda[j][k]=MATRIX<T,TV::m>::Outer_Product(z.la*z.dJ[j],z.dJ[k])+z.la*(z.J-1)*z.ddJ[j][k];
            z.ddb[j][k]=MATRIX<T,TV::m>::Outer_Product(-z.mu*z.H_bar[k],z.H_bar[j]);
            z.ddq[j][k]=(z.dw[j].Dot(z.dw[k])+MATRIX<T,TV::m>::Outer_Product(z.dw[k],z.dw[j]))/sqr(z.dt);
            z.ddc[j][k]=z.mu_N*z.Fn.Determinant()*z.ddq[j][k];
            z.ddk[j][k]=MATRIX<T,TV::m>::Outer_Product(4*z.k/(TV::m*TV::m)*z.K_bar[j],z.K_bar[k])+MATRIX<T,TV::m>::Outer_Product(2*z.k/TV::m*z.K_bar[k],z.K_bar[j])+2*z.k/TV::m*z.B_bar[j].Dot(z.K_bar[k])*z.K;
            z.ddm[j][k]=
                MATRIX<T,TV::m>::Outer_Product(4*z.m/(TV::m*TV::m)*z.H_bar[j],z.H_bar[k])-
                MATRIX<T,TV::m>::Outer_Product(2*z.m/TV::m*z.H_bar[k],z.H_bar[j]);
            z.ddn[j][k]=z.ddk[j][k]*z.g+MATRIX<T,TV::m>::Outer_Product(z.dk[j],z.dg[k])+MATRIX<T,TV::m>::Outer_Product(z.dg[j],z.dk[k]);
            z.ddp[j][k]=z.mu/2*(z.ddm[j][k]*z.n+MATRIX<T,TV::m>::Outer_Product(z.dm[j],z.dn[k])+MATRIX<T,TV::m>::Outer_Product(z.dn[j],z.dm[k])+z.ddn[j][k]*z.m);
            z.ddPhi[j][k]=z.V0*(z.ddp[j][k]-z.ddb[j][k]+z.dda[j][k]+z.ddc[j][k]);}}
}

void Test(const char* name,TV dx[2],T a,T b,TV da[2],TV db[2])
{
    T x=b-a;
    T y=(db[0]+da[0]).Dot(dx[0])/2+(db[1]+da[1]).Dot(dx[1])/2;
    LOG::printf("DIFF %s: %.16g %.16g -> %.16g\n",name,x/eps,y/eps,(x-y)/max(max(abs(x),abs(y)),1e-30));
}

void Test(const char* name,TV dx[2],TV a,TV b,MATRIX<T,TV::m> da[2],MATRIX<T,TV::m> db[2])
{
    TV x=b-a;
    TV y=(db[0]+da[0])*dx[0]/2+(db[1]+da[1])*dx[1]/2;
    LOG::printf("DIFF %s: %.16g %.16g -> %.16g\n",name,x.Magnitude()/eps,y.Magnitude()/eps,(x-y).Magnitude()/max(max(x.Magnitude(),y.Magnitude()),1e-30));
}

template<class T_TENSOR>
void Test(const char* name,TV dx[2],MATRIX<T,TV::m> a,MATRIX<T,TV::m> b,T_TENSOR da[2],T_TENSOR db[2])
{
    MATRIX<T,TV::m> x=b-a;
    MATRIX<T,TV::m> y=Contract_2(db[0]+da[0],dx[0])/2+Contract_2(db[1]+da[1],dx[1])/2;
    LOG::printf("DIFF %s: %.16g %.16g -> %.16g\n",name,x.Frobenius_Norm()/eps,y.Frobenius_Norm()/eps,(x-y).Frobenius_Norm()/max(max(x.Frobenius_Norm(),y.Frobenius_Norm()),1e-30));
}

void Test(const char* name,TV dx[2],TV a[2],TV b[2],MATRIX<T,TV::m> da[2][2],MATRIX<T,TV::m> db[2][2])
{
    VECTOR<TV,2> x,y;
    for(int i=0;i<2;i++){
        x(i)=b[i]-a[i];
        y(i)=(db[i][0]+da[i][0])*dx[0]/2+(db[i][1]+da[i][1])*dx[1]/2;}
    LOG::printf("DIFF %s: %.16g %.16g -> %.16g\n",name,x.Flattened().Magnitude()/eps,y.Flattened().Magnitude()/eps,(x-y).Flattened().Magnitude()/max(max(x.Flattened().Magnitude(),y.Flattened().Magnitude()),1e-30));
}

#define TEST(x) Test(#x,dx,z0.x,z1.x,z0.d##x,z1.d##x);

int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Parse();

    HELPER z0;
    Init(z0);

    TV dx[2];
    for(int i=0;i<2;i++) rnd.Fill_Uniform(dx[i],-eps,eps);
    HELPER z1=z0;
    for(int i=0;i<2;i++) z1.xh[i]+=dx[i];
        
    Fill(z0);
    Fill(z1);
        
    Fill_Diff(z0);
    Fill_Diff(z1);
        
    Fill_Hess(z0);
    Fill_Hess(z1);

    TEST(A);
    TEST(B);
    TEST(F_hat);
    TEST(S);
    TEST(H);
    TEST(J);
    TEST(a);
    TEST(b);
    TEST(q);
    TEST(c);
    TEST(g);
    TEST(K);
    TEST(h);
    TEST(k);
    TEST(m);
    TEST(n);
    TEST(p);
    TEST(Phi);

    TEST(dJ);
    TEST(da);
    TEST(db);
    TEST(dq);
    TEST(dc);
    TEST(dk);
    TEST(dm);
    TEST(dn);
    TEST(dp);
    TEST(dPhi);

    return 0;
}

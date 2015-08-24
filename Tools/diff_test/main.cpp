//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include "symbolic_tensor.h"

random_type rnd;

double eps=1e-6;

struct HELPER
{
    symbolic_tensor x,w,F,S;
    symbolic_tensor zeta,xi,V0,dt,J,visc;
    symbolic_tensor B,dB,ddB;
    symbolic_tensor A,dA,ddA;
    symbolic_tensor Fh,dFh,ddFh;
    symbolic_tensor Sh,dSh,ddSh;
    symbolic_tensor tm,dtm,ddtm;
    symbolic_tensor psi,dpsi,ddpsi;
    symbolic_tensor phi,dphi,ddphi;
    symbolic_tensor dpsi_dS,dpsi_dF,ddpsi_dSdS,ddpsi_dFdS,ddpsi_dFdF,ddpsi_dSdF;
};

void Init(HELPER& z)
{
    z.x.set("ia",{4,3});
    z.w.set("ia",{4,3});
    z.x.random(rnd,-1,1);
    z.w.random(rnd,-1,1);
    z.zeta.set("",{});
    z.zeta.random(rnd,1,2);
    z.xi.set("",{});
    z.xi.random(rnd,1,2);
    z.V0.set("",{});
    z.V0.random(rnd,1,2);
    z.dt.set("",{});
    z.dt.random(rnd,1,2);
    z.J.set("",{});
    z.J.random(rnd,1,2);
    z.visc.set("",{});
    z.visc.random(rnd,1,2);
    z.F.set("ab",{3,3});
    z.F.random(rnd,-1,1);
    z.S.set("ab",{3,3});
    z.S.random(rnd,-1,1);
    z.S.set("ab",z.S("ab")+z.S("ba"));
}

double c=.1;

void Fill(HELPER& z,const symbolic_tensor& dx)
{
    symbolic_tensor id;
    id.set_id("ij",3);
    z.B.set("ab",z.x("ib")*z.w("ia"));
    z.dB.set("abks",id("bs")*z.w("ka"));

    symbolic_tensor W("ab",z.w("mb")*dx("ma"));
    symbolic_tensor R("ab",W("ab")+W("ba"));
    symbolic_tensor S("ab",z.B("ab")+z.B("ba"));
    symbolic_tensor c=z.visc/(z.dt*z.dt);
    z.phi.set("",z.V0*z.J/2*(S("ab")*S("ab"))*c);
    z.dphi.set("ks",2*z.V0*z.J*S("is")*c*z.w("ki"));
    z.ddphi.set("ks",2*z.V0*z.J*R("is")*c*z.w("ki"));
}

void Test(const char* name,const symbolic_tensor& dx,const symbolic_tensor& a,
    const symbolic_tensor& b,const symbolic_tensor& da,const symbolic_tensor& db);
void Test_Hess(const char* name,const symbolic_tensor& dx,const symbolic_tensor& a,
    const symbolic_tensor& b,const symbolic_tensor& da,const symbolic_tensor& db);
#define TEST(x) Test(#x,dx,z0.x,z1.x,z0.d##x,z1.d##x);
#define TEST2(x) Test(#x,dx,z0.x,z1.x,z0.d##x,z1.d##x);Test_Hess("d"#x,dx,z0.d##x,z1.d##x,z0.dd##x,z1.dd##x);

int main(int argc, char* argv[])
{
    HELPER z0;
    Init(z0);

    symbolic_tensor dx("ks",{4,3});
    dx.random(rnd,-eps,eps);

    HELPER z1=z0;
    z1.x=z0.x+dx(z0.x.indices.c_str());

    Fill(z0,dx);
    Fill(z1,dx);

//    TEST2(A);
    // TEST2(Fh);
    // TEST2(Sh);
//    TEST2(psi);
    TEST2(phi);

    return 0;
}


void Test(const char* name,const symbolic_tensor& dx,const symbolic_tensor& a,
    const symbolic_tensor& b,const symbolic_tensor& da,const symbolic_tensor& db)
{
    std::string diff_ind;
    for(size_t i=0,k=0;k<da.indices.size();k++){
        if(a.indices[i]==da.indices[k]) i++;
        else diff_ind+=da.indices[k];}
    symbolic_tensor E(a.indices,b(a.indices)-a(a.indices));
    symbolic_tensor F(a.indices,(da(da.indices)+db(da.indices))*dx(diff_ind)/2);
    symbolic_tensor G(a.indices,E-F);
    double e=norm(E.x);
    double f=norm(F.x);
    double g=norm(G.x);
    printf("DIFF %s: %.16g %.16g -> %.16g\n",name,e/eps,f/eps,g/max(max(abs(e),abs(f)),1e-30));
}

void Test_Hess(const char* name,const symbolic_tensor& dx,const symbolic_tensor& a,
    const symbolic_tensor& b,const symbolic_tensor& da,const symbolic_tensor& db)
{
    symbolic_tensor E(a.indices,b(a.indices)-a(a.indices));
    symbolic_tensor F(a.indices,(da(a.indices)+db(a.indices))/2);
    symbolic_tensor G(a.indices,E-F);
    double e=norm(E.x);
    double f=norm(F.x);
    double g=norm(G.x);
    printf("DIFF %s: %.16g %.16g -> %.16g\n",name,e/eps,f/eps,g/max(max(abs(e),abs(f)),1e-30));
}

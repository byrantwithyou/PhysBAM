//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Vectors/VECTOR.h>

using namespace PhysBAM;

typedef float RW;
typedef double T;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,TV::m> TV_INT;

MATRIX<T,TV::m> Diff_Stencil(const ARRAY<T,FACE_INDEX<TV::m> >& u,const TV_INT& index,const TV& one_over_dX)
{
    MATRIX<T,TV::m> du;
    for(int a=0;a<TV::m;a++){
        FACE_INDEX<TV::m> f0(a,index),f1(f0);
        f1.index(a)++;
        du(a,a)=u(f1)-u(f0);
        for(int b=0;b<TV::m;b++)
            if(a!=b){
                FACE_INDEX<TV::m> f00(f0),f01(f0),f10(f1),f11(f1);
                f00.index(b)--;
                f10.index(b)--;
                f01.index(b)++;
                f11.index(b)++;
                du(a,b)=(u(f01)+u(f11)-u(f00)-u(f10))/4;}}
    return du*DIAGONAL_MATRIX<T,TV::m>(one_over_dX);
}

T Diff_Stencil_Transpose(const ARRAY<MATRIX<T,TV::m>,TV_INT>& S,const FACE_INDEX<TV::m>& face,const TV& one_over_dX)
{
    TV_INT i0=face.index,i1(i0);
    int a=face.axis;
    i1(a)--;
    T u=(S(i1)(a,a)-S(i0)(a,a))*one_over_dX(a);
    for(int b=0;b<TV::m;b++)
        if(a!=b){
            TV_INT i00(i0),i01(i0),i10(i1),i11(i1);
            i01(b)--;
            i11(b)--;
            i00(b)++;
            i10(b)++;
            u+=(S(i01)(a,b)+S(i11)(a,b)-S(i00)(a,b)-S(i10)(a,b))/4*one_over_dX(b);}
    return u;
}

T Potential_Energy(const GRID<TV>& grid,const ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT>& S,const ARRAY<T,FACE_INDEX<TV::m> >& u,T dt,T inv_Wi,T beta,const ARRAY<T,TV_INT>& V0)
{
    T alpha=dt/inv_Wi,gamma=1/(1+alpha),psi=0;
    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        MATRIX<T,TV::m> A=Diff_Stencil(u,it.index,grid.one_over_dX)*dt+1;
        T trace=gamma*MATRIX<T,TV::m>::Inner_Product(A*S(it.index),A)+TV::m*alpha*gamma;
        psi+=beta/2*V0(it.index)*trace;}
    return psi;
}

void Add_Diff(ARRAY<T,FACE_INDEX<TV::m> >& du,const GRID<TV>& grid,const ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT>& S,const ARRAY<T,FACE_INDEX<TV::m> >& u,T dt,T inv_Wi,T beta,const ARRAY<T,TV_INT>& V0)
{
    ARRAY<MATRIX<T,TV::m>,TV_INT> tmp(S.domain);
    T alpha=dt/inv_Wi,gamma=1/(1+alpha);
    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        MATRIX<T,TV::m> A=Diff_Stencil(u,it.index,grid.one_over_dX)*dt+1;
        tmp(it.index)=beta*V0(it.index)*gamma*A.Times_Transpose(S(it.index));}
    for(FACE_ITERATOR<TV> it(grid);it.Valid();it.Next())
        du(it.Full_Index())=Diff_Stencil_Transpose(tmp,it.Full_Index(),grid.one_over_dX);
}

void Add_Hess(ARRAY<T,FACE_INDEX<TV::m> >& du,const GRID<TV>& grid,const ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT>& S,const ARRAY<T,FACE_INDEX<TV::m> >& u,T dt,T inv_Wi,T beta,const ARRAY<T,TV_INT>& V0)
{
    ARRAY<MATRIX<T,TV::m>,TV_INT> tmp(S.domain);
    T alpha=dt/inv_Wi,gamma=1/(1+alpha);
    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        MATRIX<T,TV::m> Y=Diff_Stencil(u,it.index,grid.one_over_dX);
        tmp(it.index)=beta*V0(it.index)*gamma*Y.Times_Transpose(S(it.index));}
    for(FACE_ITERATOR<TV> it(grid);it.Valid();it.Next())
        du(it.Full_Index())=Diff_Stencil_Transpose(tmp,it.Full_Index(),grid.one_over_dX);
}

int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Parse();

    int ghost=1;
    T eps=1e-6,dt,inv_Wi,beta;
    RANDOM_NUMBERS<T> rand;
    rand.Fill_Uniform(dt,.1,1);
    rand.Fill_Uniform(inv_Wi,.1,1);
    rand.Fill_Uniform(beta,.1,1);

    GRID<TV> grid(TV_INT()+16,RANGE<TV>::Unit_Box(),true);
    ARRAY<T,FACE_INDEX<TV::m> > xn(grid,ghost),dx(grid,ghost),x0(grid,ghost),x1(grid,ghost),u0(grid,ghost),u1(grid,ghost);
    for(FACE_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        rand.Fill_Uniform(xn(it.Full_Index()),-1,1);
        rand.Fill_Uniform(x0(it.Full_Index()),-1,1);
        rand.Fill_Uniform(dx(it.Full_Index()),-eps,eps);}
    x1.Copy(1,x0,1,dx);
    u0.Copy(1/dt,x0,-1/dt,xn);
    u1.Copy(1/dt,x1,-1/dt,xn);

    ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT> S(grid.Cell_Indices(ghost));
    ARRAY<T,TV_INT> V0(grid.Cell_Indices(ghost));
    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        rand.Fill_Uniform(S(it.index),-1,1);
        rand.Fill_Uniform(V0(it.index),.1,1);}

    T E0=Potential_Energy(grid,S,u0,dt,inv_Wi,beta,V0);
    T E1=Potential_Energy(grid,S,u1,dt,inv_Wi,beta,V0);

    ARRAY<T,FACE_INDEX<TV::m> > dE0(grid,ghost),dE1(grid,ghost);
    Add_Diff(dE0,grid,S,u0,dt,inv_Wi,beta,V0);
    Add_Diff(dE1,grid,S,u1,dt,inv_Wi,beta,V0);

    ARRAY<T,FACE_INDEX<TV::m> > ddE(grid,ghost);
    Add_Hess(ddE,grid,S,dx,dt,inv_Wi,beta,V0);

    T dE0_dx=0,dE1_dx=0,r=0,s=0,t=0;
    for(FACE_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        T de0=dE0(it.Full_Index()),de1=dE1(it.Full_Index()),d=dx(it.Full_Index()),dde=ddE(it.Full_Index());
        dE0_dx+=de0*d;
        dE1_dx+=de1*d;
        r+=sqr(de1-de0);
        s+=sqr(dde);
        t+=sqr(de1-de0-dde);}
    r=sqrt(r)/eps;
    s=sqrt(s)/eps;
    t=sqrt(t)/eps;

    T a=(E1-E0)/eps,b=(dE0_dx+dE1_dx)/(2*eps),c=abs(a-b)/maxabs(a,b,1e-100);
    LOG::printf("%g %g %g (%g)\n",a,b,abs(a-b),c);
    LOG::printf("%g %g %g (%g)\n",r,s,t,t/maxabs(r,s,1e-100));

    return 0;
}


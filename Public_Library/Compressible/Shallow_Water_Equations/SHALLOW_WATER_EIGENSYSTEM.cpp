//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_EIGENSYSTEM
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Matrices/MATRIX.h>
#include <Compressible/Shallow_Water_Equations/SHALLOW_WATER_EIGENSYSTEM.h>
namespace PhysBAM{
//#####################################################################
// Function Flux
//#####################################################################
// F(U) for i in (-2,m+3) - 3 ghost cells
template<class TV> void SHALLOW_WATER_EIGENSYSTEM<TV>::
Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)
{
    for(int i=-3;i<m+3;i++){
        TV_DIMENSION UU=U(i),f=UU(1+a)/UU(0)*UU;
        f(1+a)+=(T).5*gravity*sqr(UU(0));
        F(i)=f;}
}
//#####################################################################
// Function Maximum_Magnitude_Eigenvalue
//#####################################################################
// maximum magnitude eigenvalue for F(U) at point cell
template<class TV> typename TV::SCALAR SHALLOW_WATER_EIGENSYSTEM<TV>::
Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell)
{
    T u=U_cell(1+a)/U_cell(0);
    T celerity=0;
    if(U_cell(0) >= 0) celerity=sqrt(gravity*U_cell(0));
    return maxabs(u-celerity,u+celerity);
}
//#####################################################################
// Function Eigenvalues
//#####################################################################
// eigenvalues for F(U) at flux i and and at points i and i+1
template<class TV> bool SHALLOW_WATER_EIGENSYSTEM<TV>::
Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,VECTOR<T,d>& lambda,VECTOR<T,d>& lambda_left,VECTOR<T,d>& lambda_right)
{
    bool is_hyperbolic=true;
    auto eig=[this,&is_hyperbolic](TV_DIMENSION UU)
        {
            T u=UU(1+a)/UU(0);
            TV_DIMENSION ev;
            ev.Fill(u);
            if(UU(0)>=0){
                T celerity=sqrt(gravity*UU(0));
                ev(0)-=celerity;
                ev(d-1)+=celerity;}
            else is_hyperbolic=false;
            return ev;
        };

    lambda_left=eig(U(i)); // eigenvalues on the left - at point i
    lambda_right=eig(U(i+1)); // eigenvalues on the right - at point i+1
    lambda=eig((T).5*(U(i)+U(i+1))); // eigenvalues in the center - at flux i

    return is_hyperbolic; // eigensystem is well defined (hyperbolic)
}
//#####################################################################
// Function Eigenvectors
//#####################################################################
// eigenvectors for F(U) at flux i between points i and i+1
template<class TV> void SHALLOW_WATER_EIGENSYSTEM<TV>::
Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R)
{
    // eigensystem in the center - at flux i
    TV_DIMENSION UU=(T).5*(U(i)+U(i+1));
    TV u=UU.Remove_Index(0)/UU(0);
    T celerity=sqrt(gravity*UU(0)),one_over_2_celerity=1/(2*celerity);

    MATRIX<T,d,d> l,r;
    l(0,0)=one_over_2_celerity*(u(a)+celerity);
    l(0,1+a)=-one_over_2_celerity;
    l(d-1,0)=-one_over_2_celerity*(u(a)-celerity);
    l(d-1,1+a)=one_over_2_celerity;
    TV_DIMENSION r0=u.Prepend(1);
    r.Set_Row(0,r0);
    r.Set_Row(d-1,r0);
    r(0,1+a)-=celerity;
    r(d-1,1+a)+=celerity;
    for(int i=0,j=1;i<TV::m;i++){
        if(i==a) continue;
        l(j,0)=-u(i);
        l(j,1+i)=1;
        r(j,1+i)=1;
        j++;}
    R=r;
    L=l;
}
//#####################################################################
template class SHALLOW_WATER_EIGENSYSTEM<VECTOR<float,1> >;
template class SHALLOW_WATER_EIGENSYSTEM<VECTOR<float,2> >;
template class SHALLOW_WATER_EIGENSYSTEM<VECTOR<double,1> >;
template class SHALLOW_WATER_EIGENSYSTEM<VECTOR<double,2> >;
}

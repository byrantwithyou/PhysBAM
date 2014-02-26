//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_NONCONSERVATIVE_FLUX_FORM_ENO  
//##################################################################### 
#ifndef __ADVECTION_NONCONSERVATIVE_FLUX_FORM_ENO__
#define __ADVECTION_NONCONSERVATIVE_FLUX_FORM_ENO__

#include <Tools/Grids_Uniform_Advection/ADVECTION_SEPARABLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input,class T2>
class ADVECTION_NONCONSERVATIVE_FLUX_FORM_ENO:public ADVECTION_SEPARABLE_UNIFORM<T_input,T2>
{
    typedef T_input T;

    int order;
public:

    ADVECTION_NONCONSERVATIVE_FLUX_FORM_ENO()
    {
        Set_Order();
    }

    void Set_Order(const int order_input=3)  PHYSBAM_OVERRIDE
    {order=order_input;assert(order >=1 && order <=3);}

//#####################################################################
    void Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx) PHYSBAM_OVERRIDE;
//#####################################################################
};  
//#####################################################################
// Function Advection_Solver
//#####################################################################
template<class T,class T2> void ADVECTION_NONCONSERVATIVE_FLUX_FORM_ENO<T,T2>::
Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx)
{
    T one_over_dx=1/dx,one_over_two_dx=.5*one_over_dx,one_over_three_dx=((T)1/3)*one_over_dx;
    ARRAY<T2,VECTOR<int,1> > D0(-2,m+3),D1(-2,m+2),D2(-2,m+1); // divided differences
    for(int i=-3;i<m+3;i++) D0(i)=Z(i); 
    if(order >= 2) for(int i=-3;i<m+2;i++) D1(i)=(D0(i+1)-D0(i))*one_over_two_dx;     
    if(order == 3) for(int i=-3;i<m+1;i++) D2(i)=(D1(i+1)-D1(i))*one_over_three_dx;

    T2 flux_left,flux_right; 
    if(order == 1) for(int i=0;i<m;i++){
        if(u(i) > 0){flux_left=D0(i-1);flux_right=D0(i);}
        else{flux_left=D0(i);flux_right=D0(i+1);}      
        u_Zx(i)=u(i)*(flux_right-flux_left)*one_over_dx;;}
    else if(order == 2) for(int i=0;i<m;i++){
        if(u(i) > 0){flux_left=ENO(dx,D0(i-1),D1(i-2),D1(i-1));flux_right=ENO(dx,D0(i),D1(i-1),D1(i));}
        else{flux_left=ENO(dx,D0(i),-D1(i),-D1(i-1));flux_right=ENO(dx,D0(i+1),-D1(i+1),-D1(i));}      
        u_Zx(i)=u(i)*(flux_right-flux_left)*one_over_dx;;}
    else if(order == 3) for(int i=0;i<m;i++){
        if(u(i) > 0){flux_left=ENO(dx,D0(i-1),D1(i-2),D1(i-1),D2(i-3),D2(i-2),D2(i-1));flux_right=ENO(dx,D0(i),D1(i-1),D1(i),D2(i-2),D2(i-2),D2(i));}
        else{flux_left=ENO(dx,D0(i),-D1(i),-D1(i-1),D2(i),D2(i-1),D2(i-2));flux_right=ENO(dx,D0(i+1),-D1(i+1),-D1(i),D2(i+1),D2(i),D2(i-1));}              
        u_Zx(i)=u(i)*(flux_right-flux_left)*one_over_dx;}
}
//#####################################################################
} 
#endif


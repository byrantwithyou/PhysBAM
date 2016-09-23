//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_HAMILTON_JACOBI_ENO  
//##################################################################### 
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_PDE/Advection/ADVECTION_HAMILTON_JACOBI_ENO.h>
using namespace PhysBAM;
//#####################################################################
// Function Advection_Solver
//#####################################################################
// Z is (-2,m+3), u and u_Zx are (1,m)
template<class TV,class T2> void ADVECTION_HAMILTON_JACOBI_ENO<TV,T2>::
Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx)
{
    return Advection_Solver(0,m,dx,Z,u,u_Zx);
}
template<class TV,class T2> void ADVECTION_HAMILTON_JACOBI_ENO<TV,T2>::
Advection_Solver(const int m_start,const int m_end,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx)
{
    int i;T one_over_dx=1/dx,one_over_two_dx=(T).5*one_over_dx,one_over_three_dx=((T)1/3)*one_over_dx;
    ARRAY<T2,VECTOR<int,1> > D0(m_start-3,m_end+2),D1(m_start-3,m_end+1),D2(m_start-3,m_end);
    for(i=m_start-3;i<m_end+2;i++) D0(i)=(Z(i+1)-Z(i))*one_over_dx; 
    if(order >= 2) for(i=m_start-3;i<m_end+1;i++) D1(i)=(D0(i+1)-D0(i))*one_over_two_dx;
    if(order == 3) for(i=m_start-3;i<m_end;i++) D2(i)=(D1(i+1)-D1(i))*one_over_three_dx;

    if(order == 1) for(i=m_start;i<m_end;i++){if(u(i) > 0) u_Zx(i)=u(i)*D0(i-1);else u_Zx(i)=u(i)*D0(i);}
    else if(order == 2) for(i=m_start;i<m_end;i++){if(u(i) > 0) u_Zx(i)=u(i)*ENO(dx,D0(i-1),D1(i-2),D1(i-1));else u_Zx(i)=u(i)*ENO(dx,D0(i),-D1(i),-D1(i-1));}
    else if(order == 3) for(i=m_start;i<m_end;i++){if(u(i) > 0)
        u_Zx(i)=u(i)*ENO(dx,D0(i-1),D1(i-2),D1(i-1),D2(i-3),D2(i-2),D2(i-1));else u_Zx(i)=u(i)*ENO(dx,D0(i),-D1(i),-D1(i-1),D2(i),D2(i-1),D2(i-2));}
}
//#####################################################################
namespace PhysBAM{
template class ADVECTION_HAMILTON_JACOBI_ENO<VECTOR<float,1>,float>;
template class ADVECTION_HAMILTON_JACOBI_ENO<VECTOR<float,2>,float>;
template class ADVECTION_HAMILTON_JACOBI_ENO<VECTOR<float,3>,float>;
template class ADVECTION_HAMILTON_JACOBI_ENO<VECTOR<double,1>,double>;
template class ADVECTION_HAMILTON_JACOBI_ENO<VECTOR<double,2>,double>;
template class ADVECTION_HAMILTON_JACOBI_ENO<VECTOR<double,3>,double>;
}

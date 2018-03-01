//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Grid_PDE/Advection/ADVECTION_CONSERVATIVE_ENO.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_input,class T2> ADVECTION_CONSERVATIVE_ENO<T_input,T2>::
ADVECTION_CONSERVATIVE_ENO()
{
    Set_Order();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_input,class T2> ADVECTION_CONSERVATIVE_ENO<T_input,T2>::
~ADVECTION_CONSERVATIVE_ENO()
{
}
//#####################################################################
// Function Set_Order
//#####################################################################
template<class T_input,class T2> void ADVECTION_CONSERVATIVE_ENO<T_input,T2>::
Set_Order(const int order_input)
{
    order=order_input;
    assert(order>=1 && order<=3);
}
//#####################################################################
// Function Advection_Solver
//#####################################################################
// finds (uZ)_x with Local Lax Friedrichs
template<class T,class T2> void ADVECTION_CONSERVATIVE_ENO<T,T2>::
Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx)
{
    T one_over_dx=1/dx,one_over_two_dx=(T)(.5*one_over_dx),one_over_three_dx=((T)1/3)*one_over_dx;
    RANGE<VECTOR<int,1> > R;
    R.min_corner.x=-2;
    R.max_corner.x=m+3;
    ARRAY<T2,VECTOR<int,1> > DZ0(R),DZ1(R),DZ2(R); // divided differences
    for(int i=-3;i<m+3;i++) DZ0(i)=Z(i); 
    if(order >= 2) for(int i=-3;i<m+2;i++) DZ1(i)=(DZ0(i+1)-DZ0(i))*one_over_two_dx;     
    if(order == 3) for(int i=-3;i<m+1;i++) DZ2(i)=(DZ1(i+1)-DZ1(i))*one_over_three_dx;
    ARRAY<T2,VECTOR<int,1> > DUZ0(R),DUZ1(R),DUZ2(R); // divided differences
    for(int i=-3;i<m+3;i++) DUZ0(i)=u(i)*Z(i); 
    if(order >= 2) for(int i=-3;i<m+2;i++) DUZ1(i)=(DUZ0(i+1)-DUZ0(i))*one_over_two_dx;     
    if(order == 3) for(int i=-3;i<m+1;i++) DUZ2(i)=(DUZ1(i+1)-DUZ1(i))*one_over_three_dx;

    ARRAY<T2,VECTOR<int,1> > flux(VECTOR<int,1>()+m); // flux is to the right of each point 
    if(order == 1) for(int i=0;i<m;i++){
        T2 alpha=maxabs(u(i),u(i+1));
        T2 flux_left=DUZ0(i)+alpha*DZ0(i);
        T2 flux_right=DUZ0(i+1)-alpha*DZ0(i+1);
        flux(i)=(T)(.5*(flux_left+flux_right));}
    else if(order == 2) for(int i=0;i<m;i++){
        T2 alpha=maxabs(u(i),u(i+1));
        T2 flux_left=ENO(dx,DUZ0(i)+alpha*DZ0(i),DUZ1(i-1)+alpha*DZ1(i-1),DUZ1(i)+alpha*DZ1(i));
        T2 flux_right=ENO(dx,DUZ0(i+1)-alpha*DZ0(i+1),-(DUZ1(i+1)-alpha*DZ1(i+1)),-(DUZ1(i)-alpha*DZ1(i)));
        flux(i)=(T)(.5*(flux_left+flux_right));}
    else if(order == 3) for(int i=0;i<m;i++){
        T2 alpha=maxabs(u(i),u(i+1));
        T2 flux_left=ENO(dx,DUZ0(i)+alpha*DZ0(i),DUZ1(i-1)+alpha*DZ1(i-1),DUZ1(i)+alpha*DZ1(i),DUZ2(i-2)+alpha*DZ2(i-2),DUZ2(i-1)+alpha*DZ2(i-1),DUZ2(i)+alpha*DZ2(i));
        T2 flux_right=ENO(dx,DUZ0(i+1)-alpha*DZ0(i+1),-(DUZ1(i+1)-alpha*DZ1(i+1)),-(DUZ1(i)-alpha*DZ1(i)),DUZ2(i+1)-alpha*DZ2(i+1),DUZ2(i)-alpha*DZ2(i),DUZ2(i-1)-alpha*DZ2(i-1));
        flux(i)=(T)(.5*(flux_left+flux_right));}

    for(int i=0;i<m;i++) u_Zx(i)=(flux(i)-flux(i-1))*one_over_dx;
}
//#####################################################################
namespace PhysBAM{
template class ADVECTION_CONSERVATIVE_ENO<VECTOR<float,1>,float>;
template class ADVECTION_CONSERVATIVE_ENO<VECTOR<float,2>,float>;
template class ADVECTION_CONSERVATIVE_ENO<VECTOR<float,3>,float>;
template class ADVECTION_CONSERVATIVE_ENO<VECTOR<double,1>,double>;
template class ADVECTION_CONSERVATIVE_ENO<VECTOR<double,2>,double>;
template class ADVECTION_CONSERVATIVE_ENO<VECTOR<double,3>,double>;
}

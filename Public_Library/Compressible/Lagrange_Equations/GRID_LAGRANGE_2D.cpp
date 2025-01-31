//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Math_Tools/sqr.h>
#include <Core/Vectors/VECTOR_3D.h>
#include <Geometry/Basic_Geometry/POLYGON.h>
#include <Compressible/Lagrange_Equations/GRID_LAGRANGE_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T> void GRID_LAGRANGE_2D<T>::
Euler_Step(const ARRAY<T,TV_INT>& u,const ARRAY<T,TV_INT>& v,const T dt)
{       
    for(int i=0;i<m;i++) for(int j=0;j<n;j++){x(i,j)+=dt*u(i,j);y(i,j)+=dt*v(i,j);}
}
//#####################################################################
// Function Get_Lengths
//#####################################################################
// L0 is (0,m-1) by (0,n) while L1 is (0,m) by (0,n-1)
template<class T> void GRID_LAGRANGE_2D<T>::
Get_Lengths(ARRAY<T,TV_INT>& L0,ARRAY<T,TV_INT>& L1)
{
    int i,j;
    for(i=0;i<m-1;i++) for(j=0;j<n;j++){VECTOR<T,3> v1(x(i,j),y(i,j),0),v2(x(i+1,j),y(i+1,j),0),v3=v2-v1;L0(i,j)=v3.Magnitude();}
    for(i=0;i<m;i++) for(j=0;j<n-1;j++){VECTOR<T,3> v1(x(i,j),y(i,j),0),v2(x(i,j+1),y(i,j+1),0),v3=v2-v1;L1(i,j)=v3.Magnitude();}
}
//#####################################################################
// Function Get_Areas
//#####################################################################
// A is (0,m-1) by (0,n-1)
// points listed counter clockwise
template<class T> void GRID_LAGRANGE_2D<T>::
Get_Areas(ARRAY<T,TV_INT>& A)
{
    POLYGON<TV> quadrilateral(4);

    for(int i=0;i<m-1;i++) for(int j=0;j<n-1;j++){
        quadrilateral.X(0)=TV(x(i,j),y(i,j));
        quadrilateral.X(1)=TV(x(i+1,j),y(i+1,j));
        quadrilateral.X(2)=TV(x(i+1,j+1),y(i+1,j+1));
        quadrilateral.X(3)=TV(x(i,j+1),y(i,j+1));
        A(i,j)=quadrilateral.Area();}  
}
//#####################################################################
// Function Get_Normals
//#####################################################################
// located at the midpoints, defined as (-y,x)/length or (y,-x)/length
// N0 points up as (0,m-1) by (0,n), N1 points right as (0,m) by (0,n-1)
template<class T> void GRID_LAGRANGE_2D<T>::
Get_Normals(ARRAY<T,TV_INT>& N1_x,ARRAY<T,TV_INT>& N1_y,ARRAY<T,TV_INT>& N2_x,ARRAY<T,TV_INT>& N2_y)
{
    int i,j;
    ARRAY<T,TV_INT> L0(TV_INT(m-1,n)),L1(TV_INT(m,n-1));
    Get_Lengths(L0,L1);
    for(i=0;i<m-1;i++) for(j=0;j<n;j++){N1_x(i,j)=-(y(i+1,j)-y(i,j))/L0(i,j);N1_y(i,j)=(x(i+1,j)-x(i,j))/L0(0,j);}
    for(i=0;i<m;i++) for(j=0;j<n-1;j++){N2_x(i,j)=(y(i,j+1)-y(i,j))/L1(i,j);N2_y(i,j)=-(x(i,j+1)-x(i,j))/L1(i,j);}
}
//#####################################################################
// Function Get_Centers
//#####################################################################
// C is (0,m-1) by (0,n-1)
template<class T> void GRID_LAGRANGE_2D<T>::
Get_Centers(ARRAY<T,TV_INT>& C_x,ARRAY<T,TV_INT>& C_y)
{
    for(int i=0;i<m-1;i++) for(int j=0;j<n-1;j++){
        C_x(i,j)=(x(i,j)+x(i+1,j)+x(i,j+1)+x(i+1,j+1))/4;C_y(i,j)=(y(i,j)+y(i+1,j)+y(i,j+1)+y(i+1,j+1))/4;}
}
//#####################################################################
// Function Get_Midpoints
//#####################################################################
// M0 is (0,m-1) by (0,n) and M1 is (0,m) by (0,n-1)
template<class T> void GRID_LAGRANGE_2D<T>::
Get_Midpoints(ARRAY<T,TV_INT>& M1_x,ARRAY<T,TV_INT>& M1_y,ARRAY<T,TV_INT>& M2_x,ARRAY<T,TV_INT>& M2_y)
{
    int i,j;
    for(i=0;i<m-1;i++) for(j=0;j<n;j++){M1_x(i,j)=(x(i,j)+x(i+1,j))/2;M1_y(i,j)=(y(i,j)+y(i+1,j))/2;}
    for(i=0;i<m;i++) for(j=0;j<n-1;j++){M2_x(i,j)=(x(i,j)+x(i,j+1))/2;M2_y(i,j)=(y(i,j)+y(i,j+1))/2;}
}
//#####################################################################
// Function Get_Sub_Zone_Lengths
//#####################################################################
// LL0 left, LL1 right, LL2 bottom, LL3 top: all (0,m-1) by (0,n-1)
template<class T> void GRID_LAGRANGE_2D<T>::
Get_Sub_Zone_Lengths(ARRAY<T,TV_INT>& LL0,ARRAY<T,TV_INT>& LL1,ARRAY<T,TV_INT>& LL2,ARRAY<T,TV_INT>& LL3)
{
    ARRAY<T,TV_INT> C_x(TV_INT(m-1,n-1)),C_y(TV_INT(m-1,n-1));Get_Centers(C_x,C_y);
    ARRAY<T,TV_INT> M1_x(TV_INT(m-1,n)),M1_y(TV_INT(m-1,n)),M2_x(TV_INT(m,n-1)),M2_y(TV_INT(m,n-1));Get_Midpoints(M1_x,M1_y,M2_x,M2_y);  
    
    for(int i=0;i<m-1;i++) for(int j=0;j<n-1;j++){
        VECTOR<T,3> v1(M2_x(i,j),M2_y(i,j),0),v2(C_x(i,j),C_y(i,j),0),v3=v2-v1;LL0(i,j)=v3.Magnitude();
        VECTOR<T,3> v4(C_x(i,j),C_y(i,j),0),v5(M2_x(i+1,j),M2_y(i+1,j),0),v6=v5-v4;LL1(i,j)=v6.Magnitude();
        VECTOR<T,3> v7(M1_x(i,j),M1_y(i,j),0),v8(C_x(i,j),C_y(i,j),0),v9=v8-v7;LL2(i,j)=v9.Magnitude();
        VECTOR<T,3> v10(C_x(i,j),C_y(i,j),0),v11(M1_x(i,j+1),M1_y(i,j+1),0),v12=v11-v10;LL3(i,j)=v12.Magnitude();}
}
//#####################################################################
// Function Get_Sub_Zone_Areas
//#####################################################################
// AA0 lower left, AA1 lower right, AA2 upper left, AA3 upper right: all (0,m-1) by (0,n-1)
// points listed counter clockwise
template<class T> void GRID_LAGRANGE_2D<T>::
Get_Sub_Zone_Areas(ARRAY<T,TV_INT>& AA0,ARRAY<T,TV_INT>& AA1,ARRAY<T,TV_INT>& AA2,ARRAY<T,TV_INT>& AA3)
{
    POLYGON<TV> quadrilateral(4);
    
    ARRAY<T,TV_INT> C_x(TV_INT(m-1,n-1)),C_y(TV_INT(m-1,n-1));Get_Centers(C_x,C_y);
    ARRAY<T,TV_INT> M1_x(TV_INT(m-1,n)),M1_y(TV_INT(m-1,n)),M2_x(TV_INT(m,n-1)),M2_y(TV_INT(m,n-1));Get_Midpoints(M1_x,M1_y,M2_x,M2_y);
    for(int i=0;i<m-1;i++) for(int j=0;j<n-1;j++){
        quadrilateral.X(0)=TV(x(i,j),y(i,j));
        quadrilateral.X(1)=TV(M1_x(i,j),M1_y(i,j));
        quadrilateral.X(2)=TV(C_x(i,j),C_y(i,j));
        quadrilateral.X(3)=TV(M2_x(i,j),M2_y(i,j));
        AA0(i,j)=quadrilateral.Area();
        quadrilateral.X(0)=TV(M1_x(i,j),M1_y(i,j));
        quadrilateral.X(1)=TV(x(i+1,j),y(i+1,j));
        quadrilateral.X(2)=TV(M2_x(i+1,j),M2_y(i+1,j));
        quadrilateral.X(3)=TV(C_x(i,j),C_y(i,j));
        AA1(i,j)=quadrilateral.Area();
        quadrilateral.X(0)=TV(M2_x(i,j),M2_y(i,j));
        quadrilateral.X(1)=TV(C_x(i,j),C_y(i,j));
        quadrilateral.X(2)=TV(M1_x(i,j+1),M1_y(i,j+1));
        quadrilateral.X(3)=TV(x(i,j+1),y(i,j+1));
        AA2(i,j)=quadrilateral.Area();
        quadrilateral.X(0)=TV(C_x(i,j),C_y(i,j));
        quadrilateral.X(1)=TV(M2_x(i+1,j),M2_y(i+1,j));
        quadrilateral.X(2)=TV(x(i+1,j+1),y(i+1,j+1));
        quadrilateral.X(3)=TV(M1_x(i,j+1),M1_y(i,j+1));
        AA3(i,j)=quadrilateral.Area();}  
}
//#####################################################################
// Function Get_Sub_Zone_Normals
//#####################################################################
// NN0 left, NN1 right, NN2 bottom, NN3 top: all (0,m-1) by (0,n-1)
// NN0 and NN1 point up, NN2 and NN3 point right
template<class T> void GRID_LAGRANGE_2D<T>::
Get_Sub_Zone_Normals(ARRAY<T,TV_INT>& NN1_x,ARRAY<T,TV_INT>& NN1_y,ARRAY<T,TV_INT>& NN2_x,ARRAY<T,TV_INT>& NN2_y,ARRAY<T,TV_INT>& NN3_x,ARRAY<T,TV_INT>& NN3_y,ARRAY<T,TV_INT>& NN4_x,ARRAY<T,TV_INT>& NN4_y)
{
    ARRAY<T,TV_INT> C_x(TV_INT(m-1,n-1)),C_y(TV_INT(m-1,n-1));
    Get_Centers(C_x,C_y);
    ARRAY<T,TV_INT> M1_x(TV_INT(m-1,n)),M1_y(TV_INT(m-1,n)),M2_x(TV_INT(m,n-1)),M2_y(TV_INT(m,n-1));
    Get_Midpoints(M1_x,M1_y,M2_x,M2_y);
    ARRAY<T,TV_INT> L0(TV_INT(m-1,n-1)),L1(TV_INT(m-1,n-1)),L2(TV_INT(m-1,n-1)),L3(TV_INT(m-1,n-1));
    Get_Sub_Zone_Lengths(L0,L1,L2,L3);
    for(int i=0;i<m-1;i++) for(int j=0;j<n-1;j++){
        NN1_x(i,j)=-(C_y(i,j)-M2_y(i,j))/L0(i,j);NN1_y(i,j)=(C_x(i,j)-M2_x(i,j))/L0(i,j);
        NN2_x(i,j)=-(M2_y(i+1,j)-C_y(i,j))/L1(i,j);NN2_y(i,j)=(M2_x(i+1,j)-C_x(i,j))/L1(i,j);
        NN3_x(i,j)=(C_y(i,j)-M1_y(i,j))/L2(i,j);NN3_y(i,j)=-(C_x(i,j)-M1_x(i,j))/L2(i,j);
        NN4_x(i,j)=(M1_y(i,j)-C_y(i,j))/L3(i,j);NN4_y(i,j)=-(M1_x(i,j)-C_x(i,j))/L3(i,j);}    
}
//#####################################################################
namespace PhysBAM{
template class GRID_LAGRANGE_2D<float>;
template class GRID_LAGRANGE_2D<double>;
}

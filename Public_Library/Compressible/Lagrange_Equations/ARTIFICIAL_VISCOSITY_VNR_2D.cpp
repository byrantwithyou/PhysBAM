//#####################################################################
// Copyright 2002, Ronald Fedkiw
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Math_Tools/max.h>
#include <Core/Math_Tools/maxabs.h>
#include <Core/Math_Tools/min.h>
#include <Core/Math_Tools/sqr.h>
#include <Compressible/Lagrange_Equations/ARTIFICIAL_VISCOSITY_VNR_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Artificial_Viscosity
//#####################################################################
template<class T> void ARTIFICIAL_VISCOSITY_VNR_2D<T>::
Get_Artificial_Viscosity(EOS<T>& eos,GRID_LAGRANGE_2D<T>& grid,const ARRAY<T,VECTOR<int,2> >& mass,const ARRAY<T,VECTOR<int,2> >& u,const ARRAY<T,VECTOR<int,2> >& v,
                         const ARRAY<T,VECTOR<int,2> >& energy,ARRAY<T,VECTOR<int,2> >& Q0,ARRAY<T,VECTOR<int,2> >& Q1,ARRAY<T,VECTOR<int,2> >& Q2,ARRAY<T,VECTOR<int,2> >& Q3)
{
    int i,j;
    int m=grid.m,n=grid.n; 
    
    // get grid information
    ARRAY<T,VECTOR<int,2> > AA0(TV_INT(m-1,n-1)),AA1(TV_INT(m-1,n-1)),AA2(TV_INT(m-1,n-1)),AA3(TV_INT(m-1,n-1));
    grid.Get_Sub_Zone_Areas(AA0,AA1,AA2,AA3);
    ARRAY<T,VECTOR<int,2> > NN1_x(TV_INT(m-1,n-1)),NN1_y(TV_INT(m-1,n-1)),NN2_x(TV_INT(m-1,n-1)),NN2_y(TV_INT(m-1,n-1)),NN3_x(TV_INT(m-1,n-1)),NN3_y(TV_INT(m-1,n-1)),NN4_x(TV_INT(m-1,n-1)),NN4_y(TV_INT(m-1,n-1));
    grid.Get_Sub_Zone_Normals(NN1_x,NN1_y,NN2_x,NN2_y,NN3_x,NN3_y,NN4_x,NN4_y);                  

    // find the density at each node
    ARRAY<T,VECTOR<int,2> > density(TV_INT(m,n)),M_node(TV_INT(m,n)),A_node(TV_INT(m,n));
    for(i=0;i<m-1;i++) for(j=0;j<n-1;j++){
        T mass_over_4=mass(i,j)/4;
        M_node(i,j)+=mass_over_4;M_node(i+1,j)+=mass_over_4;M_node(i,j+1)+=mass_over_4;M_node(i+1,j+1)+=mass_over_4;
        A_node(i,j)+=AA0(i,j);A_node(i+1,j)+=AA1(i,j);A_node(i,j+1)+=AA2(i,j);A_node(i+1,j+1)+=AA3(i,j);}
    for(i=0;i<m;i++) for(j=0;j<n;j++)
        density(i,j)=M_node(i,j)/A_node(i,j);
    
    // find jumps in the velocity
    ARRAY<T,VECTOR<int,2> > u_jump1(TV_INT(m-1,n)),v_jump1(TV_INT(m-1,n)),velocity_jump1(TV_INT(m-1,n)),V1_x(TV_INT(m-1,n)),V1_y(TV_INT(m-1,n));
    for(i=0;i<m-1;i++) for(j=0;j<n;j++){
        u_jump1(i,j)=u(i+1,j)-u(i,j);v_jump1(i,j)=v(i+1,j)-v(i,j);velocity_jump1(i,j)=sqrt(sqr(u_jump1(i,j))+sqr(v_jump1(i,j)));
        if(abs(velocity_jump1(i,j)) <= 1e-8*maxabs(u(i,j),v(i,j),u(i+1,j),v(i+1,j))){V1_x(i,j)=0;V1_y(i,j)=0;}
        else{V1_x(i,j)=u_jump1(i,j)/velocity_jump1(i,j);V1_y(i,j)=v_jump1(i,j)/velocity_jump1(i,j);}}
    ARRAY<T,VECTOR<int,2> > u_jump2(TV_INT(m,n-1)),v_jump2(TV_INT(m,n-1)),velocity_jump2(TV_INT(m,n-1)),V2_x(TV_INT(m,n-1)),V2_y(TV_INT(m,n-1));
    for(i=0;i<m;i++) for(j=0;j<n-1;j++){
        u_jump2(i,j)=u(i,j+1)-u(i,j);v_jump2(i,j)=v(i,j+1)-v(i,j);velocity_jump2(i,j)=sqrt(sqr(u_jump2(i,j))+sqr(v_jump2(i,j)));
        if(abs(velocity_jump2(i,j)) <= 1e-8*maxabs(u(i,j),v(i,j),u(i,j+1),v(i,j+1))){V2_x(i,j)=0;V2_y(i,j)=0;}
        else{V2_x(i,j)=u_jump2(i,j)/velocity_jump2(i,j);V2_y(i,j)=v_jump2(i,j)/velocity_jump2(i,j);}} 

    // compute artificial viscosities
    for(i=0;i<m-1;i++) for(j=0;j<n-1;j++){
        // bottom edge
        if(V1_x(i,j)*NN3_x(i,j)+V1_y(i,j)*NN3_y(i,j) >= 0) Q0(i,j)=0;
        else{
            T density_ave=2*density(i,j)*density(i+1,j)/(density(i,j)+density(i+1,j));
            Q0(i,j)=constant*density_ave*sqr(velocity_jump1(i,j));
            if(limiter){
                T delta_x=grid.x(i+1,j)-grid.x(i,j),delta_y=grid.y(i+1,j)-grid.y(i,j),delta=sqrt(sqr(delta_x)+sqr(delta_y));
                T ux_center=velocity_jump1(i,j)/delta;
                T r_left=1,r_right=1;
                if(ux_center != 0){
                    T delta_x_direction=delta_x/delta,delta_y_direction=delta_y/delta;
                    if(i != 1){
                        T delta_x_left=grid.x(i,j)-grid.x(i-1,j),delta_y_left=grid.y(i,j)-grid.y(i-1,j);
                        T delta_left=delta_x_left*delta_x_direction+delta_y_left*delta_y_direction;
                        if(delta_left != 0) r_left=(u_jump1(i-1,j)*V1_x(i,j)+v_jump1(i-1,j)*V1_y(i,j))/delta_left/ux_center;}
                    if(i != m-1){
                        T delta_x_right=grid.x(i+2,j)-grid.x(i+1,j),delta_y_right=grid.y(i+2,j)-grid.y(i+1,j);
                        T delta_right=delta_x_right*delta_x_direction+delta_y_right*delta_y_direction;
                        if(delta_right != 0) r_right=(u_jump1(i+1,j)*V1_x(i,j)+v_jump1(i+1,j)*V1_y(i,j))/delta_right/ux_center;}}
                T psi=max((T)0,min((r_left+r_right)/2,2*r_left,2*r_right,(T)1));
                Q0(i,j)=(1-psi)*Q0(i,j);}}
        // top edge
        if(V1_x(i,j+1)*NN4_x(i,j)+V1_y(i,j+1)*NN4_y(i,j) >= 0) Q1(i,j)=0;
        else{
            T density_ave=2*density(i,j+1)*density(i+1,j+1)/(density(i,j+1)+density(i+1,j+1));
            Q1(i,j)=constant*density_ave*sqr(velocity_jump1(i,j+1));
            if(limiter){
                T delta_x=grid.x(i+1,j+1)-grid.x(i,j+1),delta_y=grid.y(i+1,j+1)-grid.y(i,j+1),delta=sqrt(sqr(delta_x)+sqr(delta_y));
                T ux_center=velocity_jump1(i,j+1)/delta;
                T r_left=1,r_right=1;
                if(ux_center != 0){
                    T delta_x_direction=delta_x/delta,delta_y_direction=delta_y/delta;
                    if(i != 1){
                        T delta_x_left=grid.x(i,j+1)-grid.x(i-1,j+1),delta_y_left=grid.y(i,j+1)-grid.y(i-1,j+1);
                        T delta_left=delta_x_left*delta_x_direction+delta_y_left*delta_y_direction;
                        if(delta_left != 0) r_left=(u_jump1(i-1,j+1)*V1_x(i,j+1)+v_jump1(i-1,j+1)*V1_y(i,j+1))/delta_left/ux_center;}
                    if(i != m-1){
                        T delta_x_right=grid.x(i+2,j+1)-grid.x(i+1,j+1),delta_y_right=grid.y(i+2,j+1)-grid.y(i+1,j+1);
                        T delta_right=delta_x_right*delta_x_direction+delta_y_right*delta_y_direction;
                        if(delta_right != 0) r_right=(u_jump1(i+1,j+1)*V1_x(i,j+1)+v_jump1(i+1,j+1)*V1_y(i,j+1))/delta_right/ux_center;}}
                T psi=max((T)0,min((r_left+r_right)/2,2*r_left,2*r_right,(T)1));
                Q1(i,j)=(1-psi)*Q1(i,j);}}
        // left edge
        if(V2_x(i,j)*NN1_x(i,j)+V2_y(i,j)*NN1_y(i,j) >= 0) Q2(i,j)=0;
        else{
            T density_ave=2*density(i,j)*density(i,j+1)/(density(i,j)+density(i,j+1));
            Q2(i,j)=constant*density_ave*sqr(velocity_jump2(i,j));
            if(limiter){
                T delta_x=grid.x(i,j+1)-grid.x(i,j),delta_y=grid.y(i,j+1)-grid.y(i,j),delta=sqrt(sqr(delta_x)+sqr(delta_y));
                T ux_center=velocity_jump2(i,j)/delta;
                T r_bottom=1,r_top=1;
                if(ux_center != 0){
                    T delta_x_direction=delta_x/delta,delta_y_direction=delta_y/delta;
                    if(j != 1){
                        T delta_x_bottom=grid.x(i,j)-grid.x(i,j-1),delta_y_bottom=grid.y(i,j)-grid.y(i,j-1);
                        T delta_bottom=delta_x_bottom*delta_x_direction+delta_y_bottom*delta_y_direction;
                        if(delta_bottom != 0) r_bottom=(u_jump2(i,j-1)*V2_x(i,j)+v_jump2(i,j-1)*V2_y(i,j))/delta_bottom/ux_center;}
                    if(j != n-1){
                        T delta_x_top=grid.x(i,j+2)-grid.x(i,j+1),delta_y_top=grid.y(i,j+2)-grid.y(i,j+1);
                        T delta_top=delta_x_top*delta_x_direction+delta_y_top*delta_y_direction;
                        if(delta_top != 0) r_top=(u_jump2(i,j+1)*V2_x(i,j)+v_jump2(i,j+1)*V2_y(i,j))/delta_top/ux_center;}}
                T psi=max((T)0,min((r_bottom+r_top)/2,2*r_bottom,2*r_top,(T)1));
                Q2(i,j)=(1-psi)*Q2(i,j);}}
        // right edge
        if(V2_x(i+1,j)*NN2_x(i,j)+V2_y(i+1,j)*NN2_y(i,j) >= 0) Q3(i,j)=0;
        else{
            T density_ave=2*density(i+1,j)*density(i+1,j+1)/(density(i+1,j)+density(i+1,j+1));
            Q3(i,j)=constant*density_ave*sqr(velocity_jump2(i+1,j));
            if(limiter){
                T delta_x=grid.x(i+1,j+1)-grid.x(i+1,j),delta_y=grid.y(i+1,j+1)-grid.y(i+1,j),delta=sqrt(sqr(delta_x)+sqr(delta_y));
                T ux_center=velocity_jump2(i+1,j)/delta;
                T r_bottom=1,r_top=1;
                if(ux_center != 0){
                    T delta_x_direction=delta_x/delta,delta_y_direction=delta_y/delta;
                    if(j != 1){
                        T delta_x_bottom=grid.x(i+1,j)-grid.x(i+1,j-1),delta_y_bottom=grid.y(i+1,j)-grid.y(i+1,j-1);
                        T delta_bottom=delta_x_bottom*delta_x_direction+delta_y_bottom*delta_y_direction;
                        if(delta_bottom != 0) r_bottom=(u_jump2(i+1,j-1)*V2_x(i+1,j)+v_jump2(i+1,j-1)*V2_y(i+1,j))/delta_bottom/ux_center;}
                    if(j != n-1){
                        T delta_x_top=grid.x(i+1,j+2)-grid.x(i+1,j+1),delta_y_top=grid.y(i+1,j+2)-grid.y(i+1,j+1);
                        T delta_top=delta_x_top*delta_x_direction+delta_y_top*delta_y_direction;
                        if(delta_top != 0) r_top=(u_jump2(i+1,j+1)*V2_x(i+1,j)+v_jump2(i+1,j+1)*V2_y(i+1,j))/delta_top/ux_center;}}
                T psi=max((T)0,min((r_bottom+r_top)/2,2*r_bottom,2*r_top,(T)1));
                Q3(i,j)=(1-psi)*Q3(i,j);}}}
}
//#####################################################################
namespace PhysBAM{
template class ARTIFICIAL_VISCOSITY_VNR_2D<float>;
template class ARTIFICIAL_VISCOSITY_VNR_2D<double>;
}

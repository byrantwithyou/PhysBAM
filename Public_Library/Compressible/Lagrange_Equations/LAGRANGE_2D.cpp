//#####################################################################
// Copyright 2002, Ronald Fedkiw
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Math_Tools/max.h>
#include <Core/Math_Tools/maxabs.h>
#include <Core/Math_Tools/min.h>
#include <Compressible/Lagrange_Equations/LAGRANGE_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
// providing time is optional, mass is never updated
template<class T> void LAGRANGE_2D<T>::
Euler_Step(const T dt,const T time)
{  
    int i,j;
    int m=grid.m,n=grid.n;
    
    // find nodal masses
    ARRAY<T,VECTOR<int,2> > M_node(TV_INT(m,n));
    for(i=0;i<m-1;i++) for(j=0;j<n-1;j++){
        T mass_over_4=mass(i,j)/4;
        M_node(i,j)+=mass_over_4;M_node(i+1,j)+=mass_over_4;M_node(i,j+1)+=mass_over_4;M_node(i+1,j+1)+=mass_over_4;}
    
    // get grid information
    ARRAY<T,VECTOR<int,2> > L0(TV_INT(m-1,n)),L1(TV_INT(m,n-1));grid.Get_Lengths(L0,L1);
    ARRAY<T,VECTOR<int,2> > N1_x(TV_INT(m-1,n)),N1_y(TV_INT(m-1,n)),N2_x(TV_INT(m,n-1)),N2_y(TV_INT(m,n-1));grid.Get_Normals(N1_x,N1_y,N2_x,N2_y);
    ARRAY<T,VECTOR<int,2> > LL0(TV_INT(m-1,n-1)),LL1(TV_INT(m-1,n-1)),LL2(TV_INT(m-1,n-1)),LL3(TV_INT(m-1,n-1));
    grid.Get_Sub_Zone_Lengths(LL0,LL1,LL2,LL3); 
    ARRAY<T,VECTOR<int,2> > AA0(TV_INT(m-1,n-1)),AA1(TV_INT(m-1,n-1)),AA2(TV_INT(m-1,n-1)),AA3(TV_INT(m-1,n-1));
    grid.Get_Sub_Zone_Areas(AA0,AA1,AA2,AA3);
    ARRAY<T,VECTOR<int,2> > NN1_x(TV_INT(m-1,n-1)),NN1_y(TV_INT(m-1,n-1)),NN2_x(TV_INT(m-1,n-1)),NN2_y(TV_INT(m-1,n-1)),
                      NN3_x(TV_INT(m-1,n-1)),NN3_y(TV_INT(m-1,n-1)),NN4_x(TV_INT(m-1,n-1)),NN4_y(TV_INT(m-1,n-1));
    grid.Get_Sub_Zone_Normals(NN1_x,NN1_y,NN2_x,NN2_y,NN3_x,NN3_y,NN4_x,NN4_y); 

    // forces F, and heating rates H = -F dot V 
    ARRAY<T,VECTOR<int,2> > F_x(TV_INT(m,n)),F_y(TV_INT(m,n)),H(TV_INT(m-1,n-1));

    // pressue
    for(i=0;i<m-1;i++) for(j=0;j<n-1;j++){
        T mass_over_4=mass(i,j)/4;
        T pressure,force,force_x,force_y;
        // subzone 1
        pressure=eos.p(mass_over_4/AA0(i,j),energy(i,j));
        force=pressure*L0(i,j)/2;force_x=-force*N1_x(i,j);force_y=-force*N1_y(i,j);
        F_x(i,j)+=force_x;F_y(i,j)+=force_y;H(i,j)+=-force_x*u(i,j)-force_y*v(i,j);
        force=pressure*L1(i,j)/2;force_x=-force*N2_x(i,j);force_y=-force*N2_y(i,j);
        F_x(i,j)+=force_x;F_y(i,j)+=force_y;H(i,j)+=-force_x*u(i,j)-force_y*v(i,j);
        // subzone 2
        pressure=eos.p(mass_over_4/AA1(i,j),energy(i,j));
        force=pressure*L0(i,j)/2;force_x=-force*N1_x(i,j);force_y=-force*N1_y(i,j);
        F_x(i+1,j)+=force_x;F_y(i+1,j)+=force_y;H(i,j)+=-force_x*u(i+1,j)-force_y*v(i+1,j);
        force=pressure*L1(i+1,j)/2;force_x=force*N2_x(i+1,j);force_y=force*N2_y(i+1,j);
        F_x(i+1,j)+=force_x;F_y(i+1,j)+=force_y;H(i,j)+=-force_x*u(i+1,j)-force_y*v(i+1,j);
        // subzone 3
        pressure=eos.p(mass_over_4/AA2(i,j),energy(i,j));
        force=pressure*L0(i,j+1)/2;force_x=force*N1_x(i,j+1);force_y=force*N1_y(i,j+1);
        F_x(i,j+1)+=force_x;F_y(i,j+1)+=force_y;H(i,j)+=-force_x*u(i,j+1)-force_y*v(i,j+1);
        force=pressure*L1(i,j)/2;force_x=-force*N2_x(i,j);force_y=-force*N2_y(i,j);
        F_x(i,j+1)+=force_x;F_y(i,j+1)+=force_y;H(i,j)+=-force_x*u(i,j+1)-force_y*v(i,j+1);
        // subzone 4
        pressure=eos.p(mass_over_4/AA3(i,j),energy(i,j));
        force=pressure*L0(i,j+1)/2;force_x=force*N1_x(i,j+1);force_y=force*N1_y(i,j+1);
        F_x(i+1,j+1)+=force_x;F_y(i+1,j+1)+=force_y;H(i,j)+=-force_x*u(i+1,j+1)-force_y*v(i+1,j+1);
        force=pressure*L1(i+1,j)/2;force_x=force*N2_x(i+1,j);force_y=force*N2_y(i+1,j);
        F_x(i+1,j+1)+=force_x;F_y(i+1,j+1)+=force_y;H(i,j)+=-force_x*u(i+1,j+1)-force_y*v(i+1,j+1);}

    // boundary conditions - external forces 
    for(i=0;i<m;i++) for(j=0;j<n;j++){F_x(i,j)+=external_force_x(i,j);F_y(i,j)+=external_force_y(i,j);}

    // artificial viscosity  
    ARRAY<T,VECTOR<int,2> > Q0(TV_INT(m-1,n-1)),Q1(TV_INT(m-1,n-1)),Q2(TV_INT(m-1,n-1)),Q3(TV_INT(m-1,n-1));
    artificial_viscosity->Get_Artificial_Viscosity(eos,grid,mass,u,v,energy,Q0,Q1,Q2,Q3);
    // compute directions of the velocity jumps
    ARRAY<T,VECTOR<int,2> > V1_x(TV_INT(m-1,n)),V1_y(TV_INT(m-1,n));
    for(i=0;i<m-1;i++) for(j=0;j<n;j++){
        T u_jump=u(i+1,j)-u(i,j),v_jump=v(i+1,j)-v(i,j),magnitude=sqrt(sqr(u_jump)+sqr(v_jump));
        if(abs(magnitude) <= 1e-8*maxabs(u(i,j),v(i,j),u(i+1,j),v(i+1,j))){V1_x(i,j)=0;V1_y(i,j)=0;}
        else{V1_x(i,j)=u_jump/magnitude;V1_y(i,j)=v_jump/magnitude;}}
    ARRAY<T,VECTOR<int,2> > V2_x(TV_INT(m,n-1)),V2_y(TV_INT(m,n-1));
    for(i=0;i<m;i++) for(j=0;j<n-1;j++){
        T u_jump=u(i,j+1)-u(i,j),v_jump=v(i,j+1)-v(i,j),magnitude=sqrt(sqr(u_jump)+sqr(v_jump));
        if(abs(magnitude) <= 1e-8*maxabs(u(i,j),v(i,j),u(i,j+1),v(i,j+1))){V2_x(i,j)=0;V2_y(i,j)=0;}
        else{V2_x(i,j)=u_jump/magnitude;V2_y(i,j)=v_jump/magnitude;}}
    // compute artificial viscosities
    T fraction,force,force_x,force_y;
    for(i=0;i<m-1;i++) for(j=0;j<n-1;j++){
        // bottom
        fraction=V1_x(i,j)*NN3_x(i,j)+V1_y(i,j)*NN3_y(i,j);
        force=Q0(i,j)*LL2(i,j)*fraction;force_x=force*V1_x(i,j);force_y=force*V1_y(i,j);
        F_x(i,j)+=-force_x;F_y(i,j)+=-force_y;H(i,j)+=force_x*u(i,j)+force_y*v(i,j);
        F_x(i+1,j)+=force_x;F_y(i+1,j)+=force_y;H(i,j)+=-force_x*u(i+1,j)-force_y*v(i+1,j);
        // top
        fraction=V1_x(i,j+1)*NN4_x(i,j)+V1_y(i,j+1)*NN4_y(i,j);
        force=Q1(i,j)*LL3(i,j)*fraction;force_x=force*V1_x(i,j+1);force_y=force*V1_y(i,j+1);
        F_x(i,j+1)+=-force_x;F_y(i,j+1)+=-force_y;H(i,j)+=force_x*u(i,j+1)+force_y*v(i,j+1);
        F_x(i+1,j+1)+=force_x;F_y(i+1,j+1)+=force_y;H(i,j)+=-force_x*u(i+1,j+1)-force_y*v(i+1,j+1);
        // left
        fraction=V2_x(i,j)*NN1_x(i,j)+V2_y(i,j)*NN1_y(i,j);
        force=Q2(i,j)*LL0(i,j)*fraction;force_x=force*V2_x(i,j);force_y=force*V2_y(i,j);
        F_x(i,j)+=-force_x;F_y(i,j)+=-force_y;H(i,j)+=force_x*u(i,j)+force_y*v(i,j);
        F_x(i,j+1)+=force_x;F_y(i,j+1)+=force_y;H(i,j)+=-force_x*u(i,j+1)-force_y*v(i,j+1);
        // right
        fraction=V2_x(i+1,j)*NN2_x(i,j)+V2_y(i+1,j)*NN2_y(i,j);
        force=Q3(i,j)*LL1(i,j)*fraction;force_x=force*V2_x(i+1,j);force_y=force*V2_y(i+1,j);
        F_x(i+1,j)+=-force_x;F_y(i+1,j)+=-force_y;H(i,j)+=force_x*u(i+1,j)+force_y*v(i+1,j);
        F_x(i+1,j+1)+=force_x;F_y(i+1,j+1)+=force_y;H(i,j)+=-force_x*u(i+1,j+1)-force_y*v(i+1,j+1);}

    // material strength
    if(material_strength){
        T delta_x,delta_y,delta,delta_x_direction,delta_y_direction,S,S_x,S_y;
        for(i=0;i<m-1;i++) for(j=0;j<n-1;j++){
            // bottom
            delta_x=grid.x(i+1,j)-grid.x(i,j);delta_y=grid.y(i+1,j)-grid.y(i,j);delta=sqrt(sqr(delta_x)+sqr(delta_y));
            delta_x_direction=delta_x/delta;delta_y_direction=delta_y/delta;
            S=-stiffness*(L0(i,j)/L1_not(i,j)-1);S_x=S*delta_x_direction;S_y=S*delta_y_direction;
            F_x(i,j)+=-S_x;F_y(i,j)+=-S_y;H(i,j)+=S_x*u(i,j)+S_y*v(i,j);
            F_x(i+1,j)+=S_x;F_y(i+1,j)+=S_y;H(i,j)+=-S_x*u(i+1,j)-S_y*v(i+1,j);
            // top
            delta_x=grid.x(i+1,j+1)-grid.x(i,j+1);delta_y=grid.y(i+1,j+1)-grid.y(i,j+1);delta=sqrt(sqr(delta_x)+sqr(delta_y));
            delta_x_direction=delta_x/delta;delta_y_direction=delta_y/delta;
            S=-stiffness*(L0(i,j+1)/L1_not(i,j+1)-1);S_x=S*delta_x_direction;S_y=S*delta_y_direction;
            F_x(i,j+1)+=-S_x;F_y(i,j+1)+=-S_y;H(i,j)+=S_x*u(i,j+1)+S_y*v(i,j+1);
            F_x(i+1,j+1)+=S_x;F_y(i+1,j+1)+=S_y;H(i,j)+=-S_x*u(i+1,j+1)-S_y*v(i+1,j+1);
            // left
            delta_x=grid.x(i,j+1)-grid.x(i,j);delta_y=grid.y(i,j+1)-grid.y(i,j);delta=sqrt(sqr(delta_x)+sqr(delta_y));
            delta_x_direction=delta_x/delta;delta_y_direction=delta_y/delta;
            S=-stiffness*(L1(i,j)/L2_not(i,j)-1);S_x=S*delta_x_direction;S_y=S*delta_y_direction;
            F_x(i,j)+=-S_x;F_y(i,j)+=-S_y;H(i,j)+=S_x*u(i,j)+S_y*v(i,j);
            F_x(i,j+1)+=S_x;F_y(i,j+1)+=S_y;H(i,j)+=-S_x*u(i,j+1)-S_y*v(i,j+1);
            // right
            delta_x=grid.x(i+1,j+1)-grid.x(i+1,j);delta_y=grid.y(i+1,j+1)-grid.y(i+1,j);delta=sqrt(sqr(delta_x)+sqr(delta_y));
            delta_x_direction=delta_x/delta;delta_y_direction=delta_y/delta;
            S=-stiffness*(L1(i+1,j)/L2_not(i+1,j)-1);S_x=S*delta_x_direction;S_y=S*delta_y_direction;
            F_x(i+1,j)+=-S_x;F_y(i+1,j)+=-S_y;H(i,j)+=S_x*u(i+1,j)+S_y*v(i+1,j);
            F_x(i+1,j+1)+=S_x;F_y(i+1,j+1)+=S_y;H(i,j)+=-S_x*u(i+1,j+1)-S_y*v(i+1,j+1);}}
        
    // update grid
    grid.Euler_Step(u,v,dt);

    // update velocity
    for(i=0;i<m;i++) for(j=0;j<n;j++){u(i,j)+=dt*F_x(i,j)/M_node(i,j);v(i,j)+=dt*F_y(i,j)/M_node(i,j);}

    // boundary conditions - external velocities 
    for(i=0;i<m;i++) for(j=0;j<n;j++) if(fixed_velocity(i,j)){u(i,j)=external_u(i,j);v(i,j)=external_v(i,j);}

    // update interal energy 
    for(i=0;i<m-1;i++) for(j=0;j<n-1;j++) energy(i,j)+=dt*H(i,j)/mass(i,j);
}
//#####################################################################
// Function CFL
//#####################################################################
// includes the effects of artificial viscosity
template<class T> T LAGRANGE_2D<T>::
CFL()
{
    int i,j;
    int m=grid.m,n=grid.n;

    // grid
    ARRAY<T,VECTOR<int,2> > length(TV_INT(m-1,n-1));
    ARRAY<T,VECTOR<int,2> > L0(TV_INT(m-1,n)),L1(TV_INT(m,n-1));grid.Get_Lengths(L0,L1);
    ARRAY<T,VECTOR<int,2> > LL0(TV_INT(m-1,n-1)),LL1(TV_INT(m-1,n-1)),LL2(TV_INT(m-1,n-1)),LL3(TV_INT(m-1,n-1));
    grid.Get_Sub_Zone_Lengths(LL0,LL1,LL2,LL3);
    for(i=0;i<m-1;i++) for(j=0;j<n-1;j++) length(i,j)=min(L0(i,j),L0(i,j+1),L1(i,j),L1(i+1,j),LL0(i,j)+LL1(i,j),LL2(i,j)+LL3(i,j));

    // artificial viscosity
    ARRAY<T,VECTOR<int,2> > Q(TV_INT(m-1,n-1)),Q0(TV_INT(m-1,n-1)),Q1(TV_INT(m-1,n-1)),Q2(TV_INT(m-1,n-1)),Q3(TV_INT(m-1,n-1));
    artificial_viscosity->Get_Artificial_Viscosity(eos,grid,mass,u,v,energy,Q0,Q1,Q2,Q3);
    for(i=0;i<m-1;i++) for(j=0;j<n-1;j++) Q(i,j)=max(Q0(i,j),Q1(i,j),Q2(i,j),Q3(i,j));

    // material strength
    ARRAY<T,VECTOR<int,2> > S(TV_INT(m-1,n-1));
    if(material_strength){
        T S0,S1,S2,S3;
        for(i=0;i<m-1;i++) for(j=0;j<n-1;j++){
            S0=-stiffness*(L0(i,j)/L1_not(i,j)-1);S1=-stiffness*(L0(i,j+1)/L1_not(i,j+1)-1);
            S2=-stiffness*(L1(i,j)/L2_not(i,j)-1);S3=-stiffness*(L1(i+1,j)/L2_not(i+1,j)-1);
            S(i,j)=maxabs(S0/LL2(i,j),S1/LL3(i,j),S2/LL0(i,j),S3/LL1(i,j));}}
    
    // determine the time step
    ARRAY<T,VECTOR<int,2> > timestep(TV_INT(m-1,n-1));
    ARRAY<T,VECTOR<int,2> > AA0(TV_INT(m-1,n-1)),AA1(TV_INT(m-1,n-1)),AA2(TV_INT(m-1,n-1)),AA3(TV_INT(m-1,n-1));
    grid.Get_Sub_Zone_Areas(AA0,AA1,AA2,AA3);
    for(i=0;i<m-1;i++) for(j=0;j<n-1;j++){
        T mass_over_4=mass(i,j)/4;
        T density0=mass_over_4/AA0(i,j),density1=mass_over_4/AA1(i,j),
                   density2=mass_over_4/AA2(i,j),density3=mass_over_4/AA3(i,j);
        T density=min(density0,density1,density2,density3);
        T P=max(eos.p(density0,energy(i,j)),eos.p(density1,energy(i,j)),
                               eos.p(density2,energy(i,j)),eos.p(density3,energy(i,j)))+Q(i,j)+abs(S(i,j));
        T E=eos.e_From_p_And_rho(P,density);
        timestep(i,j)=eos.c(density,E)/length(i,j);}

    T dt_convect=timestep.Max();
    return 1/dt_convect;
}
//#####################################################################
namespace PhysBAM{
template class LAGRANGE_2D<float>;
template class LAGRANGE_2D<double>;
}

//#####################################################################
// Copyright 2002-2007 Doug Enright, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_EULER_EQUATIONS_OBLIQUE_ST_2D
//#####################################################################
//
//#####################################################################
// Enright - September 10, 2003
//#####################################################################
#ifndef __BOUNDARY_EULER_EQUATIONS_OBLIQUE_ST_2D__
#define __BOUNDARY_EULER_EQUATIONS_OBLIQUE_ST_2D__

#include "math.h"

namespace PhysBAM{

template<class T_GRID>
class BOUNDARY_EULER_EQUATIONS_OBLIQUE_ST_2D:public BOUNDARY_UNIFORM<T_GRID,VECTOR<typename T_GRID::SCALAR,T_GRID::dimension+2> >
{
    typedef typename T_GRID::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE::template REBIND<VECTOR<T,4> >::TYPE T_ARRAYS_T2;
    typedef VECTOR<T,2> TV;
public:
    T angle;
    bool left_constant_extrapolation;
    bool right_constant_extrapolation;

    BOUNDARY_EULER_EQUATIONS_OBLIQUE_ST_2D(T angle_input,const bool left_constant_extrapolation_input=false,const bool right_constant_extrapolation_input=false)
        :angle(angle_input)
    {
        if(left_constant_extrapolation_input) Set_Left_Constant_Extrapolation();else left_constant_extrapolation=false;
        if(right_constant_extrapolation_input) Set_Right_Constant_Extrapolation();else right_constant_extrapolation=false;  
    }

    void Set_Left_Constant_Extrapolation()
    {left_constant_extrapolation=true;}
    
    void Set_Right_Constant_Extrapolation()
    {right_constant_extrapolation=true;}

//#####################################################################
    void Fill_Ghost_Cells(const GRID<TV>& grid,const T_ARRAYS_T2& u,T_ARRAYS_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells) PHYSBAM_OVERRIDE;
    void Apply_Boundary_Condition(const GRID<TV>& grid,T_ARRAYS_T2& u,const T time) PHYSBAM_OVERRIDE;
//#####################################################################
};
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class T_GRID> void BOUNDARY_EULER_EQUATIONS_OBLIQUE_ST_2D<T_GRID>::
Fill_Ghost_Cells(const GRID<TV>& grid,const T_ARRAYS_T2& u,T_ARRAYS_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells)
{
    //if(u.length == 1){Default();return;}

    int m=grid.counts.x,n=grid.counts.y;
    int i,j;

    //amount of shifting to do . . 
    int shift = (int)floor(tan(angle)+1.e-16); //somewhat of a hack . . .

    T_ARRAYS_T2::Put(u,u_ghost); // interior

    //just keep this the same for now . . .
    if (left_constant_extrapolation) for(j=1;j<=n;j++) u_ghost(-2,j)=u_ghost(-1,j)=u_ghost(0,j)=u_ghost(1,j);
    else
        for(j=1;j<=n;j++) for(i=-2;i<=0;i++){ // left
            T rho=u_ghost(2-i,j)(1);
            T u_velocity=-u_ghost(2-i,j)(2)/u_ghost(2-i,j)(1);
            T v_velocity=u_ghost(2-i,j)(3)/u_ghost(2-i,j)(1);
            T e=u_ghost(2-i,j)(4)/u_ghost(2-i,j)(1)-(sqr(-u_velocity)+sqr(v_velocity))/2;
            u_ghost(i,j)(1)=rho;
            u_ghost(i,j)(2)=rho*u_velocity;
            u_ghost(i,j)(3)=rho*v_velocity;
            u_ghost(i,j)(4)=rho*(e+(sqr(u_velocity)+sqr(v_velocity))/2);}
    if (right_constant_extrapolation) for(j=1;j<=n;j++) u_ghost(m+3,j)=u_ghost(m+2,j)=u_ghost(m+1,j)=u_ghost(m,j);
    else
        for(j=1;j<=n;j++) for(i=m+1;i<=m+3;i++){ // right
            T rho=u_ghost(2*m-i,j)(1);
            T u_velocity=-u_ghost(2*m-i,j)(2)/u_ghost(2*m-i,j)(1);
            T v_velocity=u_ghost(2*m-i,j)(3)/u_ghost(2*m-i,j)(1);
            T e=u_ghost(2*m-i,j)(4)/u_ghost(2*m-i,j)(1)-(sqr(-u_velocity)+sqr(v_velocity))/2;
            u_ghost(i,j)(1)=rho;
            u_ghost(i,j)(2)=rho*u_velocity;
            u_ghost(i,j)(3)=rho*v_velocity;
            u_ghost(i,j)(4)=rho*(e+(sqr(u_velocity)+sqr(v_velocity))/2);}

    int shift_index = (shift != 0) ? (n-1)/shift : 0; //always n-1 grids pts separating two points . . . 
    for(i=1;i<=m;i++) {
        for(j=-2;j<=0;j++) u_ghost(i,j)=u_ghost(clamp(i+shift_index,1,m),j+n-1); //bottom
        for(j=n+1;j<=n+3;j++) u_ghost(i,j)=u_ghost(clamp(i-shift_index,1,m),j-n+1); //top
    }
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class T_GRID> void BOUNDARY_EULER_EQUATIONS_OBLIQUE_ST_2D<T_GRID>::
Apply_Boundary_Condition(const GRID<TV>& grid,T_ARRAYS_T2& u,const T time) 
{
    //if(u.length == 1){Default();return;}
    
    int m=grid.counts.x,n=grid.counts.y;
    int i,j;
    int shift = (int)floor(tan(angle)+1.e-16);

    if (!left_constant_extrapolation) 
        for(j=1;j<=n;j++){
            // left wall
            T rho=u(1,j)(1);
            T u_velocity=u(1,j)(2)/u(1,j)(1);
            T v_velocity=u(1,j)(3)/u(1,j)(1);
            T e=u(1,j)(4)/u(1,j)(1)-(sqr(u_velocity)+sqr(v_velocity))/2;
            u_velocity=0; 
            u(1,j)(1)=rho;
            u(1,j)(2)=rho*u_velocity;
            u(1,j)(3)=rho*v_velocity;
            u(1,j)(4)=rho*(e+(sqr(u_velocity)+sqr(v_velocity))/2);}

    if (!right_constant_extrapolation)
        for(j=1;j<=n;j++){
            // right wall
            T rho=u(m,j)(1);
            T u_velocity=u(m,j)(2)/u(m,j)(1);
            T v_velocity=u(m,j)(3)/u(m,j)(1);
            T e=u(m,j)(4)/u(m,j)(1)-(sqr(u_velocity)+sqr(v_velocity))/2;
            u_velocity=0;
            u(m,j)(1)=rho;
            u(m,j)(2)=rho*u_velocity;
            u(m,j)(3)=rho*v_velocity;
            u(m,j)(4)=rho*(e+(sqr(u_velocity)+sqr(v_velocity))/2);}

    int shift_index = (shift != 0) ? (n-1)/shift : 0; //assuming square grid . . . 
    for(i=1;i<=m;i++) u(i,n) = u(clamp(i-shift_index,1,m),1);
}
//#####################################################################
}
#endif

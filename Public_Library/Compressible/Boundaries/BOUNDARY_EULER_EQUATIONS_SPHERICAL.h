//#####################################################################
// Copyright 2002-2003, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_EULER_EQUATIONS_SPHERICAL
//#####################################################################
//
// Assumes a MAC grid straddling the origin on the left, applying reflection.
// Extrapolation is applied to the right.  
//
//#####################################################################
#ifndef __BOUNDARY_EULER_EQUATIONS_SPHERICAL__
#define __BOUNDARY_EULER_EQUATIONS_SPHERICAL__

#include <Core/Log/DEBUG_UTILITIES.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
namespace PhysBAM{

template<class T_input>
class BOUNDARY_EULER_EQUATIONS_SPHERICAL:public BOUNDARY<VECTOR<T_input,1>,VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T,1> TV;typedef VECTOR<T,3> TV_DIMENSION;typedef VECTOR<int,1> TV_INT;typedef ARRAYS_ND_BASE<TV_DIMENSION,VECTOR<int,1> > T_ARRAYS_DIMENSION_BASE;
public:
    BOUNDARY_EULER_EQUATIONS_SPHERICAL() = default;
    ~BOUNDARY_EULER_EQUATIONS_SPHERICAL() = default;

//#####################################################################
    void Fill_Ghost_Cells(const GRID<TV>& grid,const T_ARRAYS_DIMENSION_BASE& u,T_ARRAYS_DIMENSION_BASE& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) const override;
    void Apply_Boundary_Condition(const GRID<TV>& grid,T_ARRAYS_DIMENSION_BASE& u,const T time) const override {} // do nothing
//#####################################################################
};
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class T> void BOUNDARY_EULER_EQUATIONS_SPHERICAL<T>::
Fill_Ghost_Cells(const GRID<TV>& grid,const T_ARRAYS_DIMENSION_BASE& u,T_ARRAYS_DIMENSION_BASE& u_ghost,const T dt,const T time,const int number_of_ghost_cells) const
{
    int m=grid.counts.x;
    int i;

    T_ARRAYS_DIMENSION_BASE::Put(u,u_ghost); // interior

    // left
    for(i=-3;i<0;i++){ 
        T rho=u_ghost(-i)(0);
        T u_velocity=-u_ghost(-i)(1)/u_ghost(-i)(0);
        T e=u_ghost(-i)(2)/u_ghost(-i)(0)-sqr(u_velocity)/2;
        u_ghost(i)(0)=rho;
        u_ghost(i)(1)=rho*u_velocity;
        u_ghost(i)(2)=rho*(e+sqr(u_velocity)/2);}
        
     // right 
    u_ghost(m+2)=u_ghost(m+1)=u_ghost(m)=u_ghost(m-1);
}
//#####################################################################
}
#endif

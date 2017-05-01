//#####################################################################
// Copyright 2002-2007, Doug Enright, Jon Gretarsson, Nipun Kwatra, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC
//#####################################################################
//
// If an axis is set to be periodic then it is periodic with period m-1*grid.dx:U(m-1)=U(0), if repeats_at_last_node is true, otherwise its periodic with period m*grid.dx: U(m)=U(0).
//
//#####################################################################
#ifndef __BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC__
#define __BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC__

#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Vectors/VECTOR_3D.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class TV,class T2> // d=TV::m+2
class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC:public BOUNDARY<TV,T2>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef BOUNDARY<TV,T2> BASE;
    using BASE::Turn_Off_Constant_Extrapolation;using BASE::Set_Constant_Extrapolation;using BASE::Constant_Extrapolation;
    using BASE::constant_extrapolation;using BASE::Find_Ghost_Regions;using BASE::Fill_Single_Ghost_Region;

    VECTOR<bool,3> periodic,repeats_at_last_node;

    BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC(bool x_periodic=true,bool y_periodic=true,bool z_periodic=true,bool repeats_at_last_x_node=false,bool repeats_at_last_y_node=false,
            bool repeats_at_last_z_node=false):periodic(x_periodic,y_periodic,z_periodic),repeats_at_last_node(repeats_at_last_x_node,repeats_at_last_y_node,repeats_at_last_z_node)
    {
        Turn_Off_Constant_Extrapolation();
    }

    ~BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC() = default;

//#####################################################################
    // void Fill_Ghost_Cells_Helper(const GRID<TV>& grid,const ARRAY<VECTOR<T,3> ,VECTOR<int,1> >& u,ARRAY<VECTOR<T,3> ,VECTOR<int,1> >& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) const;
    // void Fill_Ghost_Cells_Helper(const GRID<TV>& grid,const ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& u,ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) const;
    // void Apply_Boundary_Condition_Helper(const GRID<TV>& grid,ARRAY<VECTOR<T,3> ,VECTOR<int,1> >& u,const T time) const;
    // void Apply_Boundary_Condition_Helper(const GRID<TV>& grid,ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& u,const T time) const;
    void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) const override;
    void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const T time) const override;
//#####################################################################
};
}
#endif

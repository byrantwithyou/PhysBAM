//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Frank Losasso, Avi Robinson-Mosher, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Math_Tools/maxabs.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> LEVELSET<VECTOR<T,2> >::
LEVELSET(GRID<TV>& grid_input,ARRAY<T,TV_INT>& phi_input,const int number_of_ghost_cells_input)
    :LEVELSET_BASE<TV>(grid_input,phi_input,number_of_ghost_cells_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> LEVELSET<VECTOR<T,2> >::
~LEVELSET()
{}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class T> VECTOR<T,1> LEVELSET<VECTOR<T,2> >::
Principal_Curvatures(const TV& X) const
{
    TV grad_phi=TV((Phi(TV(X.x+grid.dX.x,X.y))-Phi(TV(X.x-grid.dX.x,X.y)))/(2*grid.dX.x),
                                     (Phi(TV(X.x,X.y+grid.dX.y))-Phi(TV(X.x,X.y-grid.dX.y)))/(2*grid.dX.y));
    TV tangent=grad_phi.Perpendicular();T grad_phi_magnitude=tangent.Normalize();
    T curvature=TV::Dot_Product(tangent,Hessian(X)*tangent)/grad_phi_magnitude;
    return VECTOR<T,1>(curvature);
}
//#####################################################################
// Function Approximate_Length
//#####################################################################
// calculates the approximate perimeter using delta functions
template<class T> T LEVELSET<VECTOR<T,2> >::
Approximate_Length(const T interface_thickness,const T time) const
{
    int ghost_cells=number_of_ghost_cells;
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
    T interface_half_width=interface_thickness*grid.dX.Max()/2,one_over_two_dx=1/(2*grid.dX.x),one_over_two_dy=1/(2*grid.dX.y),length=0;
    for(int i=0;i<grid.counts.x;i++) for(int j=0;j<grid.counts.y;j++)
        length+=(LEVELSET_UTILITIES<T>::Delta(phi_ghost(i,j),interface_half_width)*sqrt(sqr((phi_ghost(i+1,j)-phi_ghost(i-1,j))*one_over_two_dx)+sqr((phi_ghost(i,j+1)-phi_ghost(i,j-1))*one_over_two_dy)));
    return length*grid.dX.x*grid.dX.y;
}
//#####################################################################
template class LEVELSET<VECTOR<float,2> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET<VECTOR<double,2> >;
#endif

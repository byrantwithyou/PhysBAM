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
template class LEVELSET<VECTOR<float,2> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET<VECTOR<double,2> >;
#endif

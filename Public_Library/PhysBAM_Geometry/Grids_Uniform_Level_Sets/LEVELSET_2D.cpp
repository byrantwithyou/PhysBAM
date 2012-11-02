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
// Function Compute_Normals
//#####################################################################
// note that sqrt(phix^2+phiy^2)=1 if it's a distance function
template<class T> void LEVELSET<VECTOR<T,2> >::
Compute_Normals(const T time)
{
    T one_over_two_dx=1/(2*grid.dX.x),one_over_two_dy=1/(2*grid.dX.y);
    int ghost_cells=number_of_ghost_cells;
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
        
    if(!normals) normals=new ARRAY<TV,TV_INT>(grid.Domain_Indices(ghost_cells-1));
    for(int i=normals->domain.min_corner.x;i<normals->domain.max_corner.x;i++) for(int j=normals->domain.min_corner.y;j<normals->domain.max_corner.y;j++){
        (*normals)(i,j)=TV((phi_ghost(i+1,j)-phi_ghost(i-1,j))*one_over_two_dx,(phi_ghost(i,j+1)-phi_ghost(i,j-1))*one_over_two_dy).Normalized();}
}
//#####################################################################
// Function Compute_Curvature
//#####################################################################
// kappa = - DIV(normal), negative for negative phi inside, positive for positive phi inside, sqrt(phix^2+phiy^2)=1 for a distance function
template<class T> void LEVELSET<VECTOR<T,2> >::
Compute_Curvature(const T time)
{
    int ghost_cells=number_of_ghost_cells;
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);

    if(!curvature) curvature=new ARRAY<T,TV_INT>(grid.Domain_Indices(2));
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid,2);iterator.Valid();iterator.Next())
        (*curvature)(iterator.Cell_Index())=Compute_Curvature(phi_ghost,iterator.Cell_Index());
}
//#####################################################################
// Function Compute_Curvature
//#####################################################################
template<class T> T LEVELSET<VECTOR<T,2> >::
Compute_Curvature(const ARRAY<T,TV_INT>& phi_input,const TV_INT& index) const
{
    T one_over_two_dx=(T).5*grid.one_over_dX.x,one_over_two_dy=(T).5*grid.one_over_dX.y;
    T one_over_dx_squared=sqr(grid.one_over_dX.x),one_over_dy_squared=sqr(grid.one_over_dX.y),one_over_four_dx_dy=(T).25*grid.one_over_dX.x*grid.one_over_dX.y;
    T max_curvature=1/grid.min_dX; // max resolution

    int i=index.x,j=index.y;
    T phix=(phi_input(i+1,j)-phi_input(i-1,j))*one_over_two_dx,phixx=(phi_input(i+1,j)-2*phi_input(i,j)+phi_input(i-1,j))*one_over_dx_squared,
        phiy=(phi_input(i,j+1)-phi_input(i,j-1))*one_over_two_dy,phiyy=(phi_input(i,j+1)-2*phi_input(i,j)+phi_input(i,j-1))*one_over_dy_squared,
        phixy=(phi_input(i+1,j+1)-phi_input(i+1,j-1)-phi_input(i-1,j+1)+phi_input(i-1,j-1))*one_over_four_dx_dy;
    T denominator=sqrt(sqr(phix)+sqr(phiy)),curvature;
    if(denominator >= small_number) curvature=-(sqr(phix)*phiyy-2*phix*phiy*phixy+sqr(phiy)*phixx)/cube(denominator);else curvature=LEVELSET_UTILITIES<T>::Sign(phi_input(i,j))*max_curvature;
    return minmag(curvature,LEVELSET_UTILITIES<T>::Sign(curvature)*max_curvature);
}
//#####################################################################
// Function Compute_Curvature
//#####################################################################
template<class T> T LEVELSET<VECTOR<T,2> >::
Compute_Curvature(const TV& location) const
{
    TV l2=(location-(T).5*grid.dX-grid.domain.min_corner)*grid.one_over_dX;
    TV_INT cell(floor(l2));

    TV w(l2-TV(cell));
    T k00=Compute_Curvature(phi,cell);
    cell.y++;
    T k01=Compute_Curvature(phi,cell);
    cell.x++;
    T k11=Compute_Curvature(phi,cell);
    cell.y--;
    T k10=Compute_Curvature(phi,cell);
    return LINEAR_INTERPOLATION<T,T>::Bilinear(k00,k10,k01,k11,w);
}
//#####################################################################
// Function Fast_Marching_Method
//#####################################################################
template<class T> void LEVELSET<VECTOR<T,2> >::
Fast_Marching_Method(const T time,const T stopping_distance,const ARRAY<TV_INT>* seed_indices,const bool add_seed_indices_for_ghost_cells,int process_sign)
{       
    Get_Signed_Distance_Using_FMM(phi,time,stopping_distance,seed_indices,add_seed_indices_for_ghost_cells,process_sign);
}
//#####################################################################
// Function Get_Signed_Distance_Using_FMM
//#####################################################################
template<class T> void LEVELSET<VECTOR<T,2> >::
Get_Signed_Distance_Using_FMM(ARRAY<T,TV_INT>& signed_distance,const T time,const T stopping_distance,const ARRAY<TV_INT>* seed_indices,const bool add_seed_indices_for_ghost_cells,int process_sign)
{       
    const int ghost_cells=2*number_of_ghost_cells+1;
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(ghost_cells),false);boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
    FAST_MARCHING_METHOD_UNIFORM<GRID<TV> > fmm(*this,ghost_cells);
    fmm.Fast_Marching_Method(phi_ghost,stopping_distance,seed_indices,add_seed_indices_for_ghost_cells,process_sign);
    ARRAY<T,TV_INT>::Get(signed_distance,phi_ghost);
    boundary->Apply_Boundary_Condition(grid,signed_distance,time);
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

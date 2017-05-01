//#####################################################################
// Copyright 2004-2007, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_PHI_WATER
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Incompressible/Boundaries/BOUNDARY_PHI_WATER.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BOUNDARY_PHI_WATER<TV>::
BOUNDARY_PHI_WATER(const TV_SIDES& constant_extrapolation)
    :sign(1),use_open_boundary_mode(false),open_boundary(6,true),V(0)
{
    Set_Constant_Extrapolation(constant_extrapolation);
    Set_Open_Boundary(false,false,false,false,false,false);
    Use_Extrapolation_Mode(false);
    Set_Tolerance();
}
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class TV> void BOUNDARY_PHI_WATER<TV>::
Fill_Ghost_Cells(const GRID<TV>& grid,const T_ARRAYS_BASE& u,T_ARRAYS_BASE& u_ghost,const T dt,const T time,const int number_of_ghost_cells) const
{
    assert(grid.Is_MAC_Grid() && V);T_ARRAYS_BASE::Put(u,u_ghost);
    RANGE<TV_INT> domain_indices=grid.Domain_Indices();
    VECTOR<RANGE<TV_INT>,2*TV::m> regions;Find_Ghost_Regions(grid,regions,number_of_ghost_cells); 
    for(int axis=0;axis<TV::m;axis++)for(int axis_side=0;axis_side<2;axis_side++){
        int side=2*axis+axis_side;
        Fill_Single_Ghost_Region(grid,u_ghost,side,regions(side));}
}
//#####################################################################
template class BOUNDARY_PHI_WATER<VECTOR<float,1> >;
template class BOUNDARY_PHI_WATER<VECTOR<float,2> >;
template class BOUNDARY_PHI_WATER<VECTOR<float,3> >;
template class BOUNDARY_PHI_WATER<VECTOR<double,1> >;
template class BOUNDARY_PHI_WATER<VECTOR<double,2> >;
template class BOUNDARY_PHI_WATER<VECTOR<double,3> >;
}

//#####################################################################
// Copyright 2005-2007, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_MPI
//#####################################################################
#ifndef __BOUNDARY_MPI__
#define __BOUNDARY_MPI__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Boundaries/BOUNDARY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_GRID_POLICY.h>
namespace PhysBAM{

template<class T_GRID> struct BOUNDARY_POLICY;

template<class TV,class T2=typename TV::SCALAR>
class BOUNDARY_MPI:public BOUNDARY<TV,T2>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;typedef VECTOR<bool,2> TV_BOOL2;typedef VECTOR<TV_BOOL2,GRID<TV>::dimension> TV_SIDES;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS;typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef typename MPI_GRID_POLICY<GRID<TV>>::MPI_GRID T_MPI_GRID;
    typedef ARRAY<T2,FACE_INDEX<TV::m> > T_FACE_ARRAYS_T2;
public:
    T_MPI_GRID* mpi_grid;
    BOUNDARY<TV,T2>& boundary;

    BOUNDARY_MPI(T_MPI_GRID* mpi_grid_input,BOUNDARY<TV,T2>& boundary_input);
    ~BOUNDARY_MPI();

    void Set_Constant_Extrapolation(const TV_SIDES& constant_extrapolation_input=TV_SIDES::Constant_Vector(TV_BOOL2::Constant_Vector(true)));
    bool Constant_Extrapolation(const int side) const PHYSBAM_OVERRIDE;
    void Set_Fixed_Boundary(const bool use_fixed_boundary_input=true,const T2 fixed_boundary_value_input=T2());
    void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells_input) PHYSBAM_OVERRIDE;
    void Fill_Ghost_Faces(const GRID<TV>& grid,const T_FACE_ARRAYS_T2& u,T_FACE_ARRAYS_T2& u_ghost,const T time,const int number_of_ghost_cells_input) PHYSBAM_OVERRIDE;
    void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const T time) PHYSBAM_OVERRIDE;
    void Apply_Boundary_Condition_Face(const GRID<TV>& grid,T_FACE_ARRAYS_T2& u,const T time) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Eran Guendelman, Avi Robinson-Mosher, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_2D  
//##################################################################### 
#ifndef __LEVELSET_2D__
#define __LEVELSET_2D__ 

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_BASE.h>
namespace PhysBAM{

template<class T>
class LEVELSET<VECTOR<T,2> >:public LEVELSET_BASE<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    typedef LEVELSET_BASE<TV> BASE;
    using BASE::grid;using BASE::phi;using BASE::normals;using BASE::curvature;using BASE::cell_range;
    using BASE::collision_body_list;using BASE::refine_fmm_initialization_with_iterative_solver;using BASE::fmm_initialization_iterations;using BASE::fmm_initialization_iterative_tolerance;
    using BASE::fmm_initialization_iterative_drift_fraction;using BASE::Hessian;using BASE::thread_queue;
    using BASE::levelset_callbacks;using BASE::small_number;using BASE::boundary;using BASE::max_time_step;using BASE::Compute_Curvature;
    using BASE::curvature_motion;using BASE::sigma;using BASE::interpolation;using BASE::Phi;using BASE::Extended_Phi;
    using BASE::curvature_interpolation;using BASE::normal_interpolation;using BASE::collision_aware_interpolation_minus;using BASE::number_of_ghost_cells;

    LEVELSET(GRID<TV>& grid_input,ARRAY<T,VECTOR<int,2> >& phi_input,const int number_of_ghost_cells_input=3);
    ~LEVELSET();

//#####################################################################
    VECTOR<T,1> Principal_Curvatures(const VECTOR<T,2>& X) const;
//#####################################################################
};
}
#endif

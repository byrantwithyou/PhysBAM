//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Andrew Selle, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_1D  
//##################################################################### 
#ifndef __LEVELSET_1D__
#define __LEVELSET_1D__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_BASE.h>
namespace PhysBAM{

template<class T>
class LEVELSET<VECTOR<T,1> >:public LEVELSET_BASE<VECTOR<T,1> >
{
    STATIC_ASSERT(IS_FLOATING_POINT<T>::value);
public:
    typedef VECTOR<T,1> TV;typedef VECTOR<int,1> TV_INT;
    typedef LEVELSET_BASE<TV> BASE;
    using BASE::grid;using BASE::phi;using BASE::normals;using BASE::curvature;using BASE::cell_range;
    using BASE::refine_fmm_initialization_with_iterative_solver;using BASE::fmm_initialization_iterations;using BASE::fmm_initialization_iterative_tolerance;
    using BASE::fmm_initialization_iterative_drift_fraction;using BASE::thread_queue;
    using BASE::levelset_callbacks;using BASE::small_number;using BASE::boundary;
    using BASE::max_time_step;using BASE::number_of_ghost_cells;using BASE::Compute_Curvature;
    using BASE::interpolation;using BASE::curvature_interpolation;using BASE::normal_interpolation;using BASE::Phi;using BASE::Extended_Phi;

    LEVELSET(GRID<TV>& grid_input,ARRAY<T,TV_INT>& phi_input,const int number_of_ghost_cells_input=3);
    ~LEVELSET();

//#####################################################################
    VECTOR<T,0> Principal_Curvatures(const TV& X) const;
public:
    void Fast_Marching_Method(const T time=0,const T stopping_distance=0,const ARRAY<TV_INT>* seed_indices=0,const bool add_seed_indices_for_ghost_cells=false,int process_sign=0);
    void Get_Signed_Distance_Using_FMM(ARRAY<T,TV_INT>& signed_distance,const T time=0,const T stopping_distance=0,const ARRAY<TV_INT>* seed_indices=0,
        const bool add_seed_indices_for_ghost_cells=false,int process_sign=0);
//#####################################################################
};
}
#endif

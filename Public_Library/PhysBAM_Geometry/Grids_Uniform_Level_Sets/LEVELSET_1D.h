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
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class LEVELSET_1D:public LEVELSET_UNIFORM<GRID<VECTOR<T_input,1> > >
{
    typedef T_input T;
    STATIC_ASSERT(IS_FLOATING_POINT<T>::value);
public:
    typedef VECTOR<T,1> TV;
    typedef LEVELSET_UNIFORM<GRID<TV> > BASE;
    using BASE::grid;using BASE::phi;using BASE::normals;using BASE::curvature;using BASE::cell_range;
    using BASE::refine_fmm_initialization_with_iterative_solver;using BASE::fmm_initialization_iterations;using BASE::fmm_initialization_iterative_tolerance;
    using BASE::fmm_initialization_iterative_drift_fraction;
    using BASE::levelset_callbacks;using BASE::small_number;using BASE::boundary;
    using BASE::max_time_step;using BASE::number_of_ghost_cells;
    using BASE::interpolation;using BASE::curvature_interpolation;using BASE::normal_interpolation;using BASE::Phi;using BASE::Extended_Phi;

    LEVELSET_1D(GRID<TV>& grid_input,ARRAY<T,VECTOR<int,1> >& phi_input,const int number_of_ghost_cells_input=3);
    ~LEVELSET_1D();

    T Compute_Curvature(const ARRAY<T,VECTOR<int,1> >& phi_input,const VECTOR<int,1>& index) const
    {return 0;}

    T Compute_Curvature(const VECTOR<int,1>& index) const
    {return Compute_Curvature(phi,index);}

    T Compute_Curvature(const VECTOR<T,1>& location) const
    {return 0;}

//#####################################################################
    VECTOR<T,0> Principal_Curvatures(const VECTOR<T,1>& X) const;
    void Compute_Normals(const T time=0);
    void Compute_Curvature(const T time=0);
    void Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists=true);
public:
    void Fast_Marching_Method(const T time=0,const T stopping_distance=0,const ARRAY<VECTOR<int,1> >* seed_indices=0,const bool add_seed_indices_for_ghost_cells=false,int process_sign=0);
    void Get_Signed_Distance_Using_FMM(ARRAY<T,VECTOR<int,1> >& signed_distance,const T time=0,const T stopping_distance=0,const ARRAY<VECTOR<int,1> >* seed_indices=0,
        const bool add_seed_indices_for_ghost_cells=false,int process_sign=0);
    void Fast_Marching_Method_Outside_Band(const T half_band_width,const T time=0,const T stopping_distance=0,int process_sign=0);
//#####################################################################
};
}
#endif

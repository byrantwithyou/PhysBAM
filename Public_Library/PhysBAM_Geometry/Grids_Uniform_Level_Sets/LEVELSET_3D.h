//#####################################################################
// Copyright 2002-2008, Zhaosheng Bao, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Frank Losasso, Neil Molino, Igor Neverov, Avi Robinson-Mosher, Jonathan Su, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_3D 
//##################################################################### 
#ifndef __LEVELSET_3D__
#define __LEVELSET_3D__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_BASE.h>
namespace PhysBAM{
template<class T> class TRIANGULATED_SURFACE;

template<class T>
class LEVELSET<VECTOR<T,3> >:public LEVELSET_BASE<VECTOR<T,3> >
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    typedef LEVELSET_BASE<TV> BASE;
    using BASE::grid;using BASE::phi;using BASE::normals;using BASE::curvature;using BASE::cell_range;
    using BASE::collision_body_list;using BASE::refine_fmm_initialization_with_iterative_solver;using BASE::fmm_initialization_iterations;using BASE::fmm_initialization_iterative_tolerance;
    using BASE::fmm_initialization_iterative_drift_fraction;using BASE::Hessian;
    using BASE::levelset_callbacks;using BASE::small_number;using BASE::boundary;using BASE::max_time_step;
    using BASE::curvature_motion;using BASE::sigma;using BASE::interpolation;using BASE::Compute_Curvature;
    using BASE::curvature_interpolation;using BASE::normal_interpolation;using BASE::secondary_interpolation;
    using BASE::collision_aware_interpolation_minus;using BASE::thread_queue;
    using BASE::number_of_ghost_cells;using BASE::Phi;using BASE::Extended_Phi;

    LEVELSET(GRID<TV>& grid_input,ARRAY<T,VECTOR<int,3> >& phi_input,const int number_of_ghost_cells_input=3);
    ~LEVELSET();

//#####################################################################
    VECTOR<T,2> Principal_Curvatures(const VECTOR<T,3>& X) const;
public:
    void Fast_Marching_Method(const T time=0,const T stopping_distance=0,const ARRAY<VECTOR<int,3> >* seed_indices=0,const bool add_seed_indices_for_ghost_cells=false,int process_sign=0);
    void Fast_Marching_Method(ARRAY<bool,VECTOR<int,3> >& seed_indices,const T time=0,const T stopping_distance=0,const bool add_seed_indices_for_ghost_cells=false,int process_sign=0);
    void Get_Signed_Distance_Using_FMM(ARRAY<T,VECTOR<int,3> >& signed_distance,const T time=0,const T stopping_distance=0,const ARRAY<VECTOR<int,3> >* seed_indices=0,
        const bool add_seed_indices_for_ghost_cells=false,int process_sign=0);
    void Get_Signed_Distance_Using_FMM(ARRAY<T,VECTOR<int,3> >& signed_distance,ARRAY<bool,VECTOR<int,3> >& seed_indices,const T time=0,const T stopping_distance=0,const bool add_seed_indices_for_ghost_cells=false,int process_sign=0);
    T Approximate_Surface_Area(const T interface_thickness=3,const T time=0) const;
    void Calculate_Triangulated_Surface_From_Marching_Tetrahedra(TRIANGULATED_SURFACE<T>& triangulated_surface,const bool include_ghost_values=false) const;
    void Calculate_Triangulated_Surface_From_Marching_Tetrahedra(const GRID<TV>& tet_grid,TRIANGULATED_SURFACE<T>& triangulated_surface) const;
private:
    int If_Zero_Crossing_Add_Particle_By_Index(TRIANGULATED_SURFACE<T>& triangulated_surface,const VECTOR<int,3>& index1,const VECTOR<int,3>& index2) const; 
    int If_Zero_Crossing_Add_Particle(TRIANGULATED_SURFACE<T>& triangulated_surface,const VECTOR<T,3>& x1,const VECTOR<T,3>& x2) const; 
    void Append_Triangles(TRIANGULATED_SURFACE<T>& triangulated_surface,const int e1,const int e2,const int e3,const int e4,const int e5,const int e6,const T phi1) const;
//#####################################################################
};   
}
#endif

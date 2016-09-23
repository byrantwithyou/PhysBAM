//#####################################################################
// Copyright 2009, Eftychios Sifakis,Aleka McAdams
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MULTIGRID_POISSON
//#####################################################################
#ifndef __MULTIGRID_POISSON__
#define __MULTIGRID_POISSON__

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Core/Utilities/NONCOPYABLE.h>
#include <Core/Vectors/VECTOR_2D.h>
#include <Core/Vectors/VECTOR_3D.h>
#include <Grid_Tools/Grids/GRID.h>

namespace PhysBAM{

template<class T,int d>
class MULTIGRID_POISSON:public NONCOPYABLE
{
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef VECTOR<int,d> TV_INT;
public:
    typedef enum {INTERIOR_CELL_TYPE=1,DIRICHLET_CELL_TYPE=2,NEUMANN_CELL_TYPE=3} CELL_TYPE;
private:

public:

    const T h;
    const T_INDEX n;
    const GRID<TV> grid;
    const GRID<TV> coarse_grid;

    ARRAY<unsigned char,TV_INT> cell_type;
    ARRAY<T,TV_INT> u;
    ARRAY<T,TV_INT> b;

    ARRAY<T,TV_INT> delta;
    ARRAY<T> one_over_diagonal_part; //defined only on boundary_indices
    ARRAY<T,TV_INT> diagonal_entries;

    int total_red_boundary_blocks,total_black_boundary_blocks;

    const int colors;
    const int& number_of_threads;

    // NOTE: in optimized version, boundary indices stored in flattened format
    // for unoptimized version, boundary indices stored in vector format
#ifndef MGPCG_UNOPTIMIZED

    ARRAY<int> boundary_block_start;
    ARRAY<int> boundary_block_end;
    ARRAY<int> boundary_indices;

    ARRAY<int> extended_boundary_indices;
#else
    ARRAY<int> boundary_block_start;
    ARRAY<int> boundary_block_end;
    ARRAY<T_INDEX> boundary_indices;

    ARRAY<T_INDEX> extended_boundary_indices;

#endif

    // for serial iterators
    RANGE<T_INDEX> padded_domain;
    RANGE<T_INDEX> unpadded_domain;
    RANGE<T_INDEX> padded_coarse_domain;
    RANGE<T_INDEX> unpadded_coarse_domain;

    enum WORKAROUND1 {boundary_block_size=4};
    enum WORKAROUND2 {boundary_block_padding=0};

    // bitmasks defined over 2x2x2 blocks 
    ARRAY<unsigned char,TV_INT> index_has_full_diagonal_coarse_bitmask;
    ARRAY<unsigned char,TV_INT> index_is_interior_coarse_bitmask;


//#####################################################################
public:
    MULTIGRID_POISSON(const T_INDEX& n,const T h,const int& number_of_threads_input);
    void Initialize();
    void Reinitialize();
    T Apply_System_Matrix(const T_INDEX& index);
    T Apply_System_Matrix(const T_INDEX& index,const ARRAY<T,TV_INT>& u_input) const;
    void Subtract_Multiple_Of_System_Matrix_And_Compute_Sum_And_Extrema(const ARRAY<T,TV_INT>& x,ARRAY<T,TV_INT>& y,double& sum,T& rmin,T& rmax) const;
    T Multiply_With_System_Matrix_And_Compute_Dot_Product(const ARRAY<T,TV_INT>&x,ARRAY<T,TV_INT>& y) const;

private:
    void Initialize_Boundary_Region();
    void Initialize_Boundary_Blocks();

    void Initialize_Interior_Bitmaps_And_Diagonal_Entries();

    void Build_System_Matrix();
    void Boundary_Relaxation(bool reverse_order,const int loops=1);
    T Interior_Relaxation(const int loops=1,const bool compute_dot_product=false);
    void Compute_Residuals_Boundary();

public:
    T Relaxation_Sweep(const bool reverse_order,const int boundary_pre_relaxations=2,const int interior_relaxations=1,const int boundary_post_relaxations=2,const bool compute_dot_product=false,const bool verbose=true);
    void Relax_And_Compute_Residuals(const int interior_relaxations,const int boundary_post_relaxations,const T nullspace_component);


    // Initialization helpers for debugging example
    void Initialize_Square_Domain();
    void Initialize_Square_Minus_Circle_Domain();
    void Initialize_Test_Domain();
    void Initialize_Test_Right_Hand_Side(ARRAY<T,TV_INT>& b_input);
    void Initialize_Test_Right_Hand_Side();
    void Initialize_Test_Initial_Guess(ARRAY<T,TV_INT>& u_input);
    void Initialize_Test_Initial_Guess();
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_UNIFORM
//#####################################################################
#ifndef __LAPLACE_UNIFORM__
#define __LAPLACE_UNIFORM__

#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Math_Tools/maxabs.h>
#include <Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <Grid_PDE/Interpolation/INTERPOLATION_UNIFORM.h>
#include <Grid_PDE/Poisson/LAPLACE.h>
#include <Grid_PDE/Poisson/LAPLACE_UNIFORM_MPI.h>
namespace PhysBAM{

template<class T> class PCG_SPARSE;
template<class T> class SPARSE_MATRIX_FLAT_MXN;
template<class TV> class LAPLACE_UNIFORM_MPI;

template<class TV>
class LAPLACE_UNIFORM:public LAPLACE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<int,TV_INT> T_ARRAYS_INT;
    typedef INTERPOLATION_UNIFORM<TV,T> T_INTERPOLATION_SCALAR;
public:
    typedef TV VECTOR_T;
    typedef GRID<TV> GRID_T;

    using LAPLACE<T>::tolerance;using LAPLACE<T>::number_of_regions;using LAPLACE<T>::solve_neumann_regions;using LAPLACE<T>::Find_Tolerance;

    GRID<TV> grid;
    ARRAY<T,TV_INT>& u;
    ARRAY<T,TV_INT> f; // f will be modified and reused as b in Ax=b for PCG
    PCG_SPARSE<T> pcg;
    T_ARRAYS_INT filled_region_colors;
    ARRAY<bool> filled_region_touches_dirichlet;
    ARRAY<bool,FACE_INDEX<TV::m> > psi_N;
    ARRAY<T,FACE_INDEX<TV::m> > psi_R;
    ARRAY<bool,TV_INT> psi_D;
    VECTOR<bool,TV::m> periodic_boundary;
    LAPLACE_UNIFORM_MPI<TV>* laplace_mpi;
    MPI_UNIFORM_GRID<TV>* mpi_grid;
    ARRAY<bool,TV_INT>* psi_D_save_for_sph;
    ARRAY<bool,FACE_INDEX<TV::m> >* psi_N_save_for_sph;
    bool enforce_compatibility;
    bool solve_single_cell_neumann_regions;
    bool use_psi_R;
public:

    LAPLACE_UNIFORM(const GRID<TV>& grid_input,ARRAY<T,TV_INT>& u_input,const bool initialize_grid,const bool enforce_compatibility_input);
    virtual ~LAPLACE_UNIFORM();

    bool All_Cell_Faces_Neumann(const TV_INT& cell_index) const
    {for(int axis=0;axis<TV::m;axis++)if(!psi_N.Component(axis)(cell_index) || !psi_N.Component(axis)(cell_index+TV_INT::Axis_Vector(axis))) return false;
    return true;}

    bool Any_Cell_Faces_Neumann(const TV_INT& cell_index) const
    {for(int axis=0;axis<TV::m;axis++)if(psi_N.Component(axis)(cell_index) || psi_N.Component(axis)(cell_index+TV_INT::Axis_Vector(axis))) return true;
    return false;}

    bool Any_Neighbor_Dirichlet(const TV_INT& cell_index) const
    {for(int axis=0;axis<TV::m;axis++){TV_INT offset=TV_INT::Axis_Vector(axis);
        if((!psi_N.Component(axis)(cell_index) && psi_D(cell_index-offset)) || (!psi_N.Component(axis)(cell_index+offset) && psi_D(cell_index+offset))) return true;}
    return false;}

//#####################################################################
    void Set_Neumann_Outer_Boundaries();
    void Set_Dirichlet_Outer_Boundaries();
    virtual void Initialize_Grid(const GRID<TV>& mac_grid_input);
    virtual void Solve(const T time=0,const bool solution_regions_already_computed=false);
    virtual void Find_A(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& A_array,ARRAY<ARRAY<T> >& b_array,const ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index);
    virtual void Find_A_Part_One(T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<int> >& row_counts);
    virtual void Find_A_Part_Two(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& A_array,ARRAY<ARRAY<T> >& b_array,T_ARRAYS_INT& cell_index_to_matrix_index);
    void Find_Solution_Regions();
    void Compute_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,T_ARRAYS_INT& cell_index_to_matrix_index);
    void Compute_Matrix_Indices(const RANGE<TV_INT>& domain,ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,T_ARRAYS_INT& cell_index_to_matrix_index);
    void Solve_Subregion(ARRAY<TV_INT>& matrix_index_to_cell_index,SPARSE_MATRIX_FLAT_MXN<T>& A,ARRAY<T>& b,const int color=0);
    void Build_Single_Solution_Region(ARRAY<bool,TV_INT>& solve);
    void Use_Psi_R();
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Jon Gretarsson, Eran Guendelman, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POISSON_UNIFORM
//#####################################################################
#ifndef __POISSON_UNIFORM__
#define __POISSON_UNIFORM__

#include <Tools/Grids_PDE_Linear/POISSON.h>
#include <Tools/Grids_Uniform_PDE_Linear/LAPLACE_UNIFORM.h>
namespace PhysBAM{

template<class T> class SPARSE_MATRIX_FLAT_MXN;

template<class TV>
class POISSON_UNIFORM:public POISSON<typename TV::SCALAR>,public LAPLACE_UNIFORM<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<int,TV_INT> T_ARRAYS_INT;
public:
    using LAPLACE_UNIFORM<TV>::grid;using POISSON<T>::use_variable_beta;using POISSON<T>::beta_given_on_faces;
    using POISSON<T>::use_weighted_divergence;using POISSON<T>::multiphase;
    using LAPLACE_UNIFORM<TV>::psi_N;using LAPLACE_UNIFORM<TV>::periodic_boundary;
    using LAPLACE_UNIFORM<TV>::filled_region_colors;using LAPLACE_UNIFORM<TV>::f;using LAPLACE_UNIFORM<TV>::u;using LAPLACE_UNIFORM<TV>::psi_D;
    using LAPLACE_UNIFORM<TV>::filled_region_touches_dirichlet;using LAPLACE_UNIFORM<TV>::solve_neumann_regions;

    ARRAY<T,FACE_INDEX<TV::m> > beta_face;
    ARRAY<T,TV_INT> variable_beta;
    ARRAY<T,FACE_INDEX<TV::m> > divergence_face_weights;
public:

    POISSON_UNIFORM(const GRID<TV>& grid_input,ARRAY<T,TV_INT>& u_input,const bool initialize_grid,const bool multiphase_input,const bool enforce_compatibility_input);
    virtual ~POISSON_UNIFORM();

    void Set_Variable_beta(const bool beta_given_on_faces_input=false)
    {use_variable_beta=true;beta_given_on_faces=beta_given_on_faces_input;
    if(!beta_given_on_faces) variable_beta.Resize(grid.Domain_Indices(1));else variable_beta.Clean_Memory();}

    void Use_Weighted_Divergence(bool use_face_velocity_weights_input=true)
    {use_weighted_divergence=use_face_velocity_weights_input;
    if(use_weighted_divergence) divergence_face_weights.Resize(grid,1);}

//#####################################################################
    void Initialize_Grid(const GRID<TV>& grid_input);
    void Find_Variable_beta();
    void Find_A_Part_Two(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& A_array,ARRAY<ARRAY<T> >& b_array,T_ARRAYS_INT& cell_index_to_matrix_index) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

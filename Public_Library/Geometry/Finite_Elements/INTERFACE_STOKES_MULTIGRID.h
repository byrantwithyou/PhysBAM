//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERFACE_STOKES_MULTIGRID
//#####################################################################
#ifndef __INTERFACE_STOKES_MULTIGRID__
#define __INTERFACE_STOKES_MULTIGRID__
#include <Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
<<<<<<< HEAD
#include <Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_COLOR.h>

=======
#include <Geometry/Finite_Elements/BOUNDARY_CONDITIONS_COLOR.h>
#include <Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_COLOR.h>
>>>>>>> aeb0b3f... Eliminate fluids dependencies in solids projects
namespace PhysBAM{

template<class TV>
class INTERFACE_STOKES_MULTIGRID:public INTERFACE_STOKES_SYSTEM_COLOR<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

public:
    enum WORKAROUND1{num_bc=3};
    typedef INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV> T_VECTOR;

    struct LEVEL
    {
        INTERFACE_STOKES_SYSTEM_COLOR<TV>* iss;
        ARRAY<SPARSE_MATRIX_FLAT_MXN<T> > pressure_poisson; // per color
        mutable T_VECTOR tmp0,tmp1,tmp2;
        ARRAY<int> interior_indices;
        ARRAY<int> boundary_indices;

        ARRAY<ARRAY<T,TV_INT> > phi_per_color;
        ARRAY<ARRAY<T,TV_INT> > phi_boundary;
        ARRAY<T,TV_INT> color_levelset_phi;
        ARRAY<int,TV_INT> color_levelset_color;
        SPARSE_MATRIX_FLAT_MXN<T> L,M;
        mutable ARRAY<KRYLOV_VECTOR_BASE<T>*> av;

        void Interior_Smoother(T_VECTOR& z,const T_VECTOR& x) const; // z should be initial guess
        void Boundary_Smoother(T_VECTOR& z,const T_VECTOR& x,int iterations) const; // z should be initial guess
        void Get_Change_Of_Variables_Matrix(SPARSE_MATRIX_FLAT_MXN<T>& M) const;
        void Exact_Solve(T_VECTOR& z,const T_VECTOR& rhs) const;
        void Initialize();
        LEVEL()
            :iss(0)
        {
        }
    };

    ARRAY<LEVEL> levels;
    ARRAY<TV_INT> p_restriction_stencil;
    VECTOR<ARRAY<TV_INT>,TV::m> u_restriction_stencil;
    int boundary_smoother_iterations;
    int number_of_ghost_cells;

    INTERFACE_STOKES_MULTIGRID(int num_levels,const GRID<TV>& grid_input,const ARRAY<T,TV_INT>& phi_value_input,
        const ARRAY<int,TV_INT>& phi_color_input,bool mac_phi,const ARRAY<ARRAY<T,TV_INT> >& phi_per_color_input,
        const VECTOR<ARRAY<T,TV_INT>,num_bc>& phi_boundary_input,int number_of_ghost_cells);
    ~INTERFACE_STOKES_MULTIGRID();

    void Construct_Level(int l);

    void Restriction(T_VECTOR& z,const T_VECTOR& x,int fine_level) const;
    void Prolongation(T_VECTOR& z,const T_VECTOR& x,int fine_level) const;

    void Apply_Preconditioner(T_VECTOR& z,const T_VECTOR& x,bool initial_guess) const;

    void Fill_Ghost(const GRID<TV>& grid,ARRAY<ARRAY<T,TV_INT> >& phi) const;
    void Fill_Ghost(const GRID<TV>& grid,ARRAY<T,TV_INT>& phi) const;
    void Coarsen_Levelset(const GRID<TV>& coarse_grid,const ARRAY<T,TV_INT>& fine_phi,ARRAY<T,TV_INT>& phi) const;
    void Coarsen_Levelset(const GRID<TV>& coarse_grid,const ARRAY<ARRAY<T,TV_INT> >& fine_phi,ARRAY<ARRAY<T,TV_INT> >& phi) const;
    void Fill_Color_Levelset(const GRID<TV>& grid,const ARRAY<ARRAY<T,TV_INT> >& cr_phis,const ARRAY<ARRAY<T,TV_INT> >& bc_phis,ARRAY<T,TV_INT>& color_phi,ARRAY<int,TV_INT>& colors) const;

    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

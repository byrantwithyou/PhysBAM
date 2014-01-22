//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERFACE_STOKES_MULTIGRID
//#####################################################################
#ifndef __INTERFACE_STOKES_MULTIGRID__
#define __INTERFACE_STOKES_MULTIGRID__
#include <Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_COLOR.h>
#include <suitesparse/umfpack.h>

namespace PhysBAM{

template<class TV>
class INTERFACE_STOKES_MULTIGRID:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

public:
    typedef INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV> T_VECTOR;

    struct LEVEL
    {
        INTERFACE_STOKES_SYSTEM_COLOR<TV>* iss;
        ARRAY<SPARSE_MATRIX_FLAT_MXN<T> > pressure_poisson; // per color
        T_VECTOR tmp0,tmp1,tmp2;
        ARRAY<int> interior_indices;
        ARRAY<int> boundary_indices;

        ARRAY<ARRAY<T,TV_INT> > phi_per_color;
        ARRAY<ARRAY<T,TV_INT> > phi_boundary;
        ARRAY<T,TV_INT> color_levelset_phi;
        ARRAY<int,TV_INT> color_levelset_color;
        
        void Interior_Smoother(T_VECTOR& z,const T_VECTOR& x) const; // z should be initial guess
        void Boundary_Smoother(T_VECTOR& z,const T_VECTOR& x,int iterations) const; // z should be initial guess
        void Get_Change_Of_Variables_Matrix(SPARSE_MATRIX_FLAT_MXN<T>& M) const;
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

    INTERFACE_STOKES_MULTIGRID(int num_levels,INTERFACE_STOKES_SYSTEM_COLOR<TV>* iss,const ARRAY<ARRAY<T,TV_INT> >& phi_per_color_input,const ARRAY<ARRAY<T,TV_INT> >& phi_boundary_input,int number_of_ghost_cells);
    ~INTERFACE_STOKES_MULTIGRID();

    void Construct_Level(int l);

    void Restriction(T_VECTOR& z,const T_VECTOR& x,int fine_level) const;
    void Prolongation(T_VECTOR& z,const T_VECTOR& x,int fine_level) const;
    void Exact_Solve(T_VECTOR& z,const T_VECTOR& rhs) const;

    void Apply_Preconditioner(T_VECTOR& z,const T_VECTOR& x,bool initial_guess);
    void Update();

    void Fill_Ghost(const GRID<TV>& grid,ARRAY<ARRAY<T,TV_INT> >& phi) const;
    void Fill_Ghost(const GRID<TV>& grid,ARRAY<T,TV_INT>& phi) const;
    void Coarsen_Levelset(const GRID<TV>& coarse_grid,const ARRAY<T,TV_INT>& fine_phi,ARRAY<T,TV_INT>& phi) const;
    void Coarsen_Levelset(const GRID<TV>& coarse_grid,const ARRAY<ARRAY<T,TV_INT> >& fine_phi,ARRAY<ARRAY<T,TV_INT> >& phi) const;
    void Fill_Color_Levelset(const GRID<TV>& grid,const ARRAY<ARRAY<T,TV_INT> >& cr_phis,const ARRAY<ARRAY<T,TV_INT> >& bc_phis,ARRAY<T,TV_INT>& color_phi,ARRAY<int,TV_INT>& colors) const;
    
//#####################################################################
};
}
#endif

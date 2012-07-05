//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERFACE_STOKES_SYSTEM_COLOR
//#####################################################################
#ifndef __INTERFACE_STOKES_SYSTEM_COLOR__
#define __INTERFACE_STOKES_SYSTEM_COLOR__
#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_VECTOR_COLOR.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Finite_Elements/BOUNDARY_CONDITIONS_COLOR.h>

namespace PhysBAM{

template<class TV> class GRID;
template<class T_GRID> class LEVELSET_UNIFORM;
template<class TV> class CELL_MANAGER_COLOR;
template<class TV> class CELL_DOMAIN_INTERFACE_COLOR;

template<class TV>
class INTERFACE_STOKES_SYSTEM_COLOR:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV> VECTOR_T;
    typedef KRYLOV_SYSTEM_BASE<T> BASE;

    VECTOR_T J; // Jacobi preconditioner 

public:

    //   #----# #----# #----#   #----#   #----# #----# #----# 
    //   | UU | | UV | | UW |   | UP |   | UN | | UT | | US |
    //   #----# #----# #----#   #----#   #----# #----# #----#
    //   #----# #----# #----#   #----#   #----# #----# #----# 
    //   | VU | | VV | | VW |   | VP |   | VN | | VT | | VS |
    //   #----# #----# #----#   #----#   #----# #----# #----#
    //   #----# #----# #----#   #----#   #----# #----# #----# 
    //   | WU | | WV | | WW |   | WP |   | WN | | WT | | WS |
    //   #----# #----# #----#   #----#   #----# #----# #----#
    //
    //   #----# #----# #----#   #---------------------------# 
    //   | PU | | PV | | PW |   |                           |
    //   #----# #----# #----#   |                           |
    //                          |                           |
    //   #----# #----# #----#   |                           |
    //   | NU | | NV | | NW |   |                           |
    //   #----# #----# #----#   |                           |
    //   #----# #----# #----#   |                           |
    //   | TV | | TV | | TW |   |                           |
    //   #----# #----# #----#   |                           |
    //   #----# #----# #----#   |                           |
    //   | SU | | SV | | SW |   |                           |
    //   #----# #----# #----#   #---------------------------# 
    
    VECTOR<VECTOR<ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >,TV::m>,TV::m> matrix_uu;
    VECTOR<ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >,TV::m> matrix_pu;
    VECTOR<ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >,TV::m> matrix_qu;

    VECTOR<ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >,TV::m> matrix_rhs_pu;
    VECTOR<ARRAY<VECTOR_ND<T> >,TV::m> rhs_surface;

    VECTOR<VECTOR_T,TV::m> null_u;
    VECTOR_T null_p;
    VECTOR<ARRAY<ARRAY<int> >,TV::m> inactive_u;
    ARRAY<ARRAY<int> > inactive_p;
    ARRAY<int> inactive_q;

    const GRID<TV>& grid;

    GRID<TV> phi_grid;
    ARRAY<T,TV_INT> phi_value;
    ARRAY<int,TV_INT> phi_color;
    
    bool run_self_tests;
    bool print_matrix;
    bool print_rhs;

    static int solve_id;

    VECTOR<CELL_MANAGER_COLOR<TV>*,TV::m> cm_u;
    CELL_MANAGER_COLOR<TV> *cm_p;
    CELL_DOMAIN_INTERFACE_COLOR<TV> *cdi;

    INTERFACE_STOKES_SYSTEM_COLOR(const GRID<TV>& grid_input,const ARRAY<T,TV_INT>& phi_value_input,const ARRAY<int,TV_INT>& phi_color_input);
    virtual ~INTERFACE_STOKES_SYSTEM_COLOR();

//#####################################################################
    void Set_Matrix(const ARRAY<T>& mu,bool wrap,BOUNDARY_CONDITIONS_COLOR<TV>* abc);
    void Set_RHS(VECTOR_T& rhs,const ARRAY<ARRAY<TV,TV_INT> >& f_volume,const ARRAY<ARRAY<T,FACE_INDEX<TV::m> > >& u);
    void Resize_Vector(KRYLOV_VECTOR_BASE<T>& x) const;
    void Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const;
    void Project(KRYLOV_VECTOR_BASE<T>& x) const;
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const;
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const;
private:
    void Set_Jacobi_Preconditioner();
//#####################################################################
};
}
#endif

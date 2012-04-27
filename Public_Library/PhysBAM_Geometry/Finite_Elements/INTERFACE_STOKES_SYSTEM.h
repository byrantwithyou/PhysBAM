//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERFACE_STOKES_SYSTEM
//#####################################################################
#ifndef __INTERFACE_STOKES_SYSTEM__
#define __INTERFACE_STOKES_SYSTEM__
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER.h>
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_VECTOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_UNIFORM.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

namespace PhysBAM{

template<class TV> class GRID;
template<class TV>
class INTERFACE_STOKES_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef INTERFACE_STOKES_SYSTEM_VECTOR<TV> VECTOR_T;
    typedef KRYLOV_SYSTEM_BASE<T> BASE;

    void Set_Jacobi_Preconditioner();

public:

    //   #----# #----# #----#   #----#   #----# #----# #----# 
    //   | UU | | UV | | UW |   | UP |   | UQ | |    | |    |
    //   #----# #----# #----#   #----#   #----# #----# #----#
    //   #----# #----# #----#   #----#   #----# #----# #----# 
    //   | VU | | VV | | VW |   | VP |   |    | | VR | |    |
    //   #----# #----# #----#   #----#   #----# #----# #----#
    //   #----# #----# #----#   #----#   #----# #----# #----# 
    //   | WU | | WV | | WW |   | WP |   |    | |    | | WS |
    //   #----# #----# #----#   #----#   #----# #----# #----#
    //
    //   #----# #----# #----#   #---------------------------# 
    //   | PU | | PV | | PW |   |                           |
    //   #----# #----# #----#   |                           |
    //                          |                           |
    //   #----# #----# #----#   |                           |
    //   | QU | |    | |    |   |                           |
    //   #----# #----# #----#   |                           |
    //   #----# #----# #----#   |                           |
    //   |    | | RV | |    |   |                           |
    //   #----# #----# #----#   |                           |
    //   #----# #----# #----#   |                           |
    //   |    | |    | | SW |   |                           |
    //   #----# #----# #----#   #---------------------------# 
    
    VECTOR<VECTOR<VECTOR<SPARSE_MATRIX_FLAT_MXN<T>,2>,TV::m>,TV::m> matrix_uu;
    VECTOR<VECTOR<SPARSE_MATRIX_FLAT_MXN<T>,2>,TV::m> matrix_pu;
    VECTOR<VECTOR<SPARSE_MATRIX_FLAT_MXN<T>,2>,TV::m> matrix_qu;
    VECTOR<VECTOR<SPARSE_MATRIX_FLAT_MXN<T>,2>,TV::m> matrix_rhs_qu;

    VECTOR<VECTOR_T,TV::m> null_u;
    VECTOR_T null_p;

    VECTOR_T J; // Jacobi preconditioner 

    const GRID<TV>& grid;
    const GRID<TV>& coarse_grid;
    const bool periodic_bc;
    
    LEVELSET_UNIFORM<GRID<TV> >* phi;
    GRID<TV> phi_grid;

    bool run_self_tests;
    bool print_matrix;
    bool print_rhs;

    static int solve_id;

    typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT object;
    VECTOR<CELL_MANAGER<TV>*,TV::m> cm_u;
    CELL_MANAGER<TV> *cm_p;
    CELL_DOMAIN_INTERFACE<TV> *cdi;

    INTERFACE_STOKES_SYSTEM(const GRID<TV>& grid_input,GRID<TV>& coarse_grid_input,ARRAY<T,TV_INT>& phi_input,bool periodic_bc_input=true);
    virtual ~INTERFACE_STOKES_SYSTEM();

//#####################################################################
    void Set_Matrix(const VECTOR<T,2>& mu);
    void Set_RHS(VECTOR_T& rhs,const VECTOR<ARRAY<TV,TV_INT>,2> f_body,const ARRAY<TV>& f_interface);
    void Resize_Vector(KRYLOV_VECTOR_BASE<T>& x) const;
    void Get_U_Part(const VECTOR_T& x,ARRAY<T,FACE_INDEX<TV::m> >& u) const;
    void Get_P_Part(const VECTOR_T& x,ARRAY<T,TV_INT>& p) const;
    void Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const;
    void Project(KRYLOV_VECTOR_BASE<T>& x) const;
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const;
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const;
//#####################################################################
};
}
#endif

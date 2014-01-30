//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERFACE_STOKES_SYSTEM_COLOR
//#####################################################################
#ifndef __INTERFACE_STOKES_SYSTEM_COLOR__
#define __INTERFACE_STOKES_SYSTEM_COLOR__
#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_VECTOR_COLOR.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <boost/function.hpp>

namespace PhysBAM{

template<class TV> class GRID;
template<class TV> class CELL_MANAGER_COLOR;
template<class TV> class CELL_DOMAIN_INTERFACE_COLOR;
template<class TV> struct VOLUME_FORCE_COLOR;

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

    VECTOR<ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >,TV::m> matrix_rhs_pu,matrix_inertial_rhs;
    VECTOR<ARRAY<ARRAY<T> >,TV::m> rhs_surface;
    ARRAY<T> q_rhs;
    VECTOR<VECTOR<ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >,TV::m>,TV::m> matrix_polymer_stress_rhs;

    ARRAY<VECTOR_T*> null_modes;
    VECTOR<ARRAY<ARRAY<int> >,TV::m> inactive_u;
    ARRAY<ARRAY<int> > inactive_p;
    ARRAY<int> inactive_q;

    GRID<TV> grid;

    GRID<TV> phi_grid;
    ARRAY<T,TV_INT> phi_value;
    ARRAY<int,TV_INT> phi_color;
    
    bool run_self_tests;
    bool print_matrix;
    bool print_rhs;
    bool use_p_null_mode;
    bool use_u_null_mode;
    bool use_polymer_stress;

    static int solve_id;

    VECTOR<CELL_MANAGER_COLOR<TV>*,TV::m> cm_u;
    CELL_MANAGER_COLOR<TV> *cm_p;
    CELL_MANAGER_COLOR<TV> *cm_ps;
    CELL_DOMAIN_INTERFACE_COLOR<TV> *cdi;

    INTERFACE_STOKES_SYSTEM_COLOR(const GRID<TV>& grid_input,const ARRAY<T,TV_INT>& phi_value_input,const ARRAY<int,TV_INT>& phi_color_input,bool mac_phi);
    virtual ~INTERFACE_STOKES_SYSTEM_COLOR();

//#####################################################################
    virtual void Set_Matrix(const ARRAY<T>& mu,bool use_discontinuous_velocity,
        boost::function<TV(const TV& X,int color0,int color1)> u_jump,
        boost::function<TV(const TV& X,int color0,int color1)> j_surface,
        ARRAY<T>* inertia,bool use_rhs);
    virtual void Set_RHS(VECTOR_T& rhs,boost::function<TV(const TV& X,int color)> body_force,const ARRAY<ARRAY<T,FACE_INDEX<TV::m> > >* u,bool analytic_velocity_correction);
    void Add_Polymer_Stress_RHS(VECTOR_T& rhs,const ARRAY<ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT> >& polymer_stress,T dt);
    void Resize_Vector(KRYLOV_VECTOR_BASE<T>& x) const;
    void Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const;
    void Project(KRYLOV_VECTOR_BASE<T>& x) const;
    void Clear_Unused_Entries(VECTOR_T& v) const;
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const;
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const;
    void Pack(const ARRAY<ARRAY<T,FACE_INDEX<TV::m> > >& u,VECTOR<ARRAY<ARRAY<T> >,TV::m>& v) const;
    void Pack(const ARRAY<ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT> >& polymer_stress,VECTOR<VECTOR<ARRAY<ARRAY<T> >,TV::m>,TV::m>& S,T scale) const;
    void Get_Sparse_Matrix(SPARSE_MATRIX_FLAT_MXN<T>& M) const;
private:
    void Set_Jacobi_Preconditioner();
//#####################################################################
};
}
#endif

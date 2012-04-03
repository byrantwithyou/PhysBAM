//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERFACE_FLUID_SYSTEM
//#####################################################################
#ifndef __INTERFACE_FLUID_SYSTEM__
#define __INTERFACE_FLUID_SYSTEM__
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>

namespace PhysBAM{

template<class TV> class GRID;
template<class TV>
class INTERFACE_FLUID_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef KRYLOV_VECTOR_WRAPPER<T,VECTOR_ND<T> > VECTOR_T;
    typedef KRYLOV_SYSTEM_BASE<T> BASE;

    SYSTEM_MATRIX_HELPER<T> *helper_rhs_q[TV::m],*helper_rhs_p[TV::m];

public:

    SPARSE_MATRIX_FLAT_MXN<T> matrix;
    const GRID<TV>& grid;
    const GRID<TV>& coarse_grid;
    const ARRAY<T,TV_INT>& phi;
    VECTOR_ND<T> null[TV::m],null_p;
    ARRAY<TV_INT> cell_map;
    ARRAY<bool> sign_map;
    bool run_self_tests;
    bool print_matrix;
    bool print_rhs;
    static int solve_id;
    int system_size;
    VECTOR_T rhs;
    VECTOR_T solution;
    typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT object;
    INTERVAL<int> index_range_u[TV::m],index_range_p,index_range_q[TV::m];

    INTERFACE_FLUID_SYSTEM(const GRID<TV>& grid_input,const GRID<TV>& coarse_grid_input,const ARRAY<T,TV_INT>& phi_input);
    virtual ~INTERFACE_FLUID_SYSTEM();

//#####################################################################
    void Set_Matrix(const VECTOR<T,2>& mu);
    void Set_RHS(VECTOR_T& rhs, const ARRAY<TV,TV_INT> f_body[2],const ARRAY<TV>& f_interface);
    void Get_U_Part(const VECTOR_ND<T>& x,ARRAY<T,TV_INT>& u_part,const int dir) const;
    void Get_P_Part(const VECTOR_ND<T>& x,ARRAY<T,TV_INT>& p_part) const;
    void Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const;
    void Project(KRYLOV_VECTOR_BASE<T>& x) const;
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const;
//#####################################################################
};
}
#endif

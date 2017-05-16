//#####################################################################
// Copyright 2016, Lin Huang, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace MPM_PROJECTION_SYSTEM
//#####################################################################
#ifndef __MPM_PROJECTION_SYSTEM__
#define __MPM_PROJECTION_SYSTEM__

#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
#include <Hybrid_Methods/Projection/MPM_PROJECTION_VECTOR.h>
#include <Hybrid_Methods/Examples_And_Drivers/PHASE_ID.h>

namespace PhysBAM{
//#####################################################################
// Class MPM_PROJECTION_SYSTEM
//#####################################################################
template<class TV>
class MPM_PROJECTION_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
public:
    typedef typename TV::SCALAR T;
    SPARSE_MATRIX_FLAT_MXN<T> A;
    mutable ARRAY<T> temp_vector;

    SPARSE_MATRIX_FLAT_MXN<T> gradient;
    SPARSE_MATRIX_FLAT_MXN<T> divergence;
    ARRAY<T> mass;
    ARRAY<FACE_INDEX<TV::m> > faces;
    ARRAY<PHASE_ID> phases;
    bool dc_present;
    MPM_PROJECTION_VECTOR<TV> null_u;
    
    MPM_PROJECTION_SYSTEM();
    virtual ~MPM_PROJECTION_SYSTEM();

    void Compute_Ones_Nullspace();
    void Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const override;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const override;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const override;
    void Project(KRYLOV_VECTOR_BASE<T>& x) const override;
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const override;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const override;
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const override;
    void Test() const;
};
//#####################################################################
}
#endif

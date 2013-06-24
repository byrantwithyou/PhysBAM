//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BW_BACKWARD_EULER_SYSTEM
//#####################################################################
#ifndef __BW_BACKWARD_EULER_SYSTEM__
#define __BW_BACKWARD_EULER_SYSTEM__

#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Tools/Matrices/MATRIX_POLICY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_MASS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
namespace PhysBAM{
template<class TV> class SOLIDS_EVOLUTION;
template<class TV> class BW_COLLISIONS;
template<class TV> class EXAMPLE_FORCES_AND_VELOCITIES;
//#####################################################################
// Class BW_BACKWARD_EULER_SYSTEM
//#####################################################################
template<class TV>
class BW_BACKWARD_EULER_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
public:
    typedef GENERALIZED_VELOCITY<TV> VECTOR_T;

    SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    BW_COLLISIONS<TV>& bw_collisions;
    EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities;
    T dt,time;

    ARRAY<TV> F_full,B_full;
    ARRAY<TWIST<TV> > rigid_F_full,rigid_B_full;

    BW_BACKWARD_EULER_SYSTEM(SOLID_BODY_COLLECTION<TV>& solid_body_collection,BW_COLLISIONS<TV>& bw_collisions_input,
        EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities,const T dt_input,const T time_input);

    virtual ~BW_BACKWARD_EULER_SYSTEM();

//#####################################################################
    void Force(const VECTOR_T& V,VECTOR_T& F) const;
    void Multiply(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& F) const PHYSBAM_OVERRIDE;
    void Project_Helper(KRYLOV_VECTOR_BASE<T>& V,const bool negate) const;
    void Project(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& V1,const KRYLOV_VECTOR_BASE<T>& V2) const PHYSBAM_OVERRIDE;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& R) const PHYSBAM_OVERRIDE;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const PHYSBAM_OVERRIDE;
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_VISCOSITY_UNIFORM_SYSTEM
//#####################################################################
#ifndef __LEVELSET_VISCOSITY_UNIFORM_SYSTEM__
#define __LEVELSET_VISCOSITY_UNIFORM_SYSTEM__

#include <Core/Vectors/VECTOR.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Incompressible/Incompressible_Flows/LEVELSET_FACE_POISSON_UNIFORM.h>
namespace PhysBAM{

template<class TV> class BOUNDARY_CONDITIONS_CALLBACKS;
template<class TV> class LEVELSET_INDEX_MAP_UNIFORM;

template<class TV>
class LEVELSET_VISCOSITY_UNIFORM_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
public:
    typedef KRYLOV_VECTOR_WRAPPER<T,ARRAY<T> > VECTOR_T;

    LEVELSET_FACE_POISSON_UNIFORM<TV> poisson;
    T scale;

    LEVELSET_VISCOSITY_UNIFORM_SYSTEM(LEVELSET_INDEX_MAP_UNIFORM<TV>& index_map_input);
    virtual ~LEVELSET_VISCOSITY_UNIFORM_SYSTEM();

//#####################################################################
    void Compute(int axis,T scale_input);
    void Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result,bool transpose=false) const override;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const override;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const override;
    void Project(KRYLOV_VECTOR_BASE<T>& x) const override;
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const override;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const override;
    void Add_Constant_Part(KRYLOV_VECTOR_BASE<T>& x) const;
//#####################################################################
};
}
#endif

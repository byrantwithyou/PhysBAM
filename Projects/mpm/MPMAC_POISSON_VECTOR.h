//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPMAC_POISSON_VECTOR
//#####################################################################
#ifndef __MPMAC_POISSON_VECTOR__
#define __MPMAC_POISSON_VECTOR__

#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
namespace PhysBAM{

template<class TV>
class MPMAC_POISSON_VECTOR:public KRYLOV_VECTOR_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    ARRAY<T,TV_INT> v;

    MPMAC_POISSON_VECTOR();
    virtual ~MPMAC_POISSON_VECTOR();

    const MPMAC_POISSON_VECTOR& operator=(const MPMAC_POISSON_VECTOR& bv);
    KRYLOV_VECTOR_BASE<T>& operator+=(const KRYLOV_VECTOR_BASE<T>& bv) PHYSBAM_OVERRIDE;
    KRYLOV_VECTOR_BASE<T>& operator-=(const KRYLOV_VECTOR_BASE<T>& bv) PHYSBAM_OVERRIDE;
    KRYLOV_VECTOR_BASE<T>& operator*=(const T a) PHYSBAM_OVERRIDE;
    void Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv) PHYSBAM_OVERRIDE;
    void Copy(const T c1,const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2) PHYSBAM_OVERRIDE;
    T Dot(const KRYLOV_VECTOR_BASE<T>& bv) const PHYSBAM_OVERRIDE;
    int Raw_Size() const PHYSBAM_OVERRIDE;
    T& Raw_Get(int i) PHYSBAM_OVERRIDE;
    KRYLOV_VECTOR_BASE<T>* Clone_Default() const PHYSBAM_OVERRIDE;
    void Resize(const KRYLOV_VECTOR_BASE<T>& v) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_KRYLOV_VECTOR
//#####################################################################
#ifndef __MPM_KRYLOV_VECTOR__
#define __MPM_KRYLOV_VECTOR__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
namespace PhysBAM{

template<class TV> class TWIST;
template<class TV>
class MPM_KRYLOV_VECTOR:public KRYLOV_VECTOR_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    ARRAY<TV,TV_INT> u;
    ARRAY<int>& valid_indices;

    MPM_KRYLOV_VECTOR(ARRAY<int>& valid_indices);
    virtual ~MPM_KRYLOV_VECTOR();

    const MPM_KRYLOV_VECTOR& operator=(const MPM_KRYLOV_VECTOR& bv);
    KRYLOV_VECTOR_BASE<T>& operator+=(const KRYLOV_VECTOR_BASE<T>& bv) PHYSBAM_OVERRIDE;
    KRYLOV_VECTOR_BASE<T>& operator-=(const KRYLOV_VECTOR_BASE<T>& bv) PHYSBAM_OVERRIDE;
    KRYLOV_VECTOR_BASE<T>& operator*=(const T a) PHYSBAM_OVERRIDE;
    void Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv) PHYSBAM_OVERRIDE;
    void Copy(const T c1,const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2) PHYSBAM_OVERRIDE;
    int Raw_Size() const PHYSBAM_OVERRIDE;
    T& Raw_Get(int i) PHYSBAM_OVERRIDE;
    KRYLOV_VECTOR_BASE<T>* Clone_Default() const PHYSBAM_OVERRIDE;
    void Resize(const KRYLOV_VECTOR_BASE<T>& v) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

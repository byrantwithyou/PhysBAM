//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_KRYLOV_VECTOR_RB
//#####################################################################
#ifndef __MPM_KRYLOV_VECTOR_RB__
#define __MPM_KRYLOV_VECTOR_RB__

#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Vectors/TWIST.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
namespace PhysBAM{

template<class TV> class TWIST;
template<class TV>
class MPM_KRYLOV_VECTOR_RB:public KRYLOV_VECTOR_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    ARRAY<TV,TV_INT> u;
    ARRAY<int>& valid_indices;
    ARRAY<TWIST<TV> > twists;

    MPM_KRYLOV_VECTOR_RB(ARRAY<int>& valid_indices);
    virtual ~MPM_KRYLOV_VECTOR_RB();

    const MPM_KRYLOV_VECTOR_RB& operator=(const MPM_KRYLOV_VECTOR_RB& bv);
    KRYLOV_VECTOR_BASE<T>& operator+=(const KRYLOV_VECTOR_BASE<T>& bv) override;
    KRYLOV_VECTOR_BASE<T>& operator-=(const KRYLOV_VECTOR_BASE<T>& bv) override;
    KRYLOV_VECTOR_BASE<T>& operator*=(const T a) override;
    void Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv) override;
    void Copy(const T c1,const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2) override;
    int Raw_Size() const override;
    T& Raw_Get(int i) override;
    KRYLOV_VECTOR_BASE<T>* Clone_Default() const override;
    void Resize(const KRYLOV_VECTOR_BASE<T>& v) override;
    void Get(ARRAY_VIEW<T> a) const override;
    void Set(ARRAY_VIEW<const T> a) override;
//#####################################################################
};
}
#endif

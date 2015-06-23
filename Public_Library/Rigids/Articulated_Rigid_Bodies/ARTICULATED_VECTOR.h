//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARTICULATED_VECTOR
//#####################################################################
#ifndef __ARTICULATED_VECTOR__
#define __ARTICULATED_VECTOR__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Rigids/Joints/JOINT_ID.h>
namespace PhysBAM{

template<class TV> class TWIST;
template<class TV>
class ARTICULATED_VECTOR:public KRYLOV_VECTOR_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
public:
    ARRAY<TWIST<TV>,JOINT_ID> v;

    ARTICULATED_VECTOR();
    virtual ~ARTICULATED_VECTOR();

    const ARTICULATED_VECTOR& operator=(const ARTICULATED_VECTOR& bv);
    KRYLOV_VECTOR_BASE<T>& operator+=(const KRYLOV_VECTOR_BASE<T>& bv) override;
    KRYLOV_VECTOR_BASE<T>& operator-=(const KRYLOV_VECTOR_BASE<T>& bv) override;
    KRYLOV_VECTOR_BASE<T>& operator*=(const T a) override;
    void Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv) override;
    void Copy(const T c1,const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2) override;
    int Raw_Size() const override;
    T& Raw_Get(int i) override;
    KRYLOV_VECTOR_BASE<T>* Clone_Default() const override;
    void Resize(const KRYLOV_VECTOR_BASE<T>& v) override;
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2009, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KRYLOV_VECTOR_WRAPPER
//#####################################################################
#ifndef __KRYLOV_VECTOR_WRAPPER__
#define __KRYLOV_VECTOR_WRAPPER__

#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
namespace PhysBAM{

template<class T,class TV>
class KRYLOV_VECTOR_WRAPPER:public KRYLOV_VECTOR_BASE<T>
{
public:
    TV v;
    bool deep_copy;

    KRYLOV_VECTOR_WRAPPER();
    KRYLOV_VECTOR_WRAPPER(TV vector);
    template<class VECTOR,class INDICES> KRYLOV_VECTOR_WRAPPER(VECTOR& vector,const INDICES& index);
    virtual ~KRYLOV_VECTOR_WRAPPER();

    KRYLOV_VECTOR_BASE<T>& operator+=(const KRYLOV_VECTOR_BASE<T>& bv) override;
    KRYLOV_VECTOR_BASE<T>& operator-=(const KRYLOV_VECTOR_BASE<T>& bv) override;
    KRYLOV_VECTOR_BASE<T>& operator*=(const T a) override;
    void Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv) override;
    void Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2) override;
    int Raw_Size() const override;
    T& Raw_Get(int i) override;
    KRYLOV_VECTOR_BASE<T>* Clone_Default() const override;
    void Resize(const KRYLOV_VECTOR_BASE<T>& v) override;
//#####################################################################
};
}
#endif

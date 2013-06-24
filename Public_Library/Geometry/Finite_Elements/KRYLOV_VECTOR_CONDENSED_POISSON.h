//#####################################################################
// Copyright 2009, Craig Schroeder,Russell Howes.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KRYLOV_VECTOR_CONDENSED_POISSON
//#####################################################################
#ifndef __KRYLOV_VECTOR_CONDENSED_POISSON__
#define __KRYLOV_VECTOR_CONDENSED_POISSON__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV>
class KRYLOV_VECTOR_CONDENSED_POISSON:public KRYLOV_VECTOR_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
public:
    ARRAY<T> v;

    KRYLOV_VECTOR_CONDENSED_POISSON();
    KRYLOV_VECTOR_CONDENSED_POISSON(TV vector);
    template<class VECTOR,class INDICES> KRYLOV_VECTOR_CONDENSED_POISSON(VECTOR& vector,const INDICES& index);
    virtual ~KRYLOV_VECTOR_CONDENSED_POISSON();

    KRYLOV_VECTOR_BASE<T>& operator+=(const KRYLOV_VECTOR_BASE<T>& bv) PHYSBAM_OVERRIDE;
    KRYLOV_VECTOR_BASE<T>& operator-=(const KRYLOV_VECTOR_BASE<T>& bv) PHYSBAM_OVERRIDE;
    KRYLOV_VECTOR_BASE<T>& operator*=(const T a) PHYSBAM_OVERRIDE;
    void Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv) PHYSBAM_OVERRIDE;
    void Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2) PHYSBAM_OVERRIDE;
    T Dot(const KRYLOV_VECTOR_BASE<T>& bv) const PHYSBAM_OVERRIDE;
    int Raw_Size() const PHYSBAM_OVERRIDE;
    T& Raw_Get(int i) PHYSBAM_OVERRIDE;
    KRYLOV_VECTOR_BASE<T>* Clone_Default() const PHYSBAM_OVERRIDE;
    void Resize(const KRYLOV_VECTOR_BASE<T>& w) PHYSBAM_OVERRIDE;
    T Max_Abs() const;
//#####################################################################
};
}
#endif

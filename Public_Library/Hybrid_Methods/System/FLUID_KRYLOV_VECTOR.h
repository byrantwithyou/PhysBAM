//#####################################################################
// Copyright 2015, Greg Klar, Andre Pradhana
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_KRYLOV_VECTOR
//#####################################################################
#ifndef __FLUID_KRYLOV_VECTOR__
#define __FLUID_KRYLOV_VECTOR__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
namespace PhysBAM{

template<class TV> class TWIST;
template<class TV>
class FLUID_KRYLOV_VECTOR:public KRYLOV_VECTOR_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    ARRAY<T,TV_INT> p;
    ARRAY<int>& valid_indices;

    FLUID_KRYLOV_VECTOR(ARRAY<int>& valid_indices);
    virtual ~FLUID_KRYLOV_VECTOR();

    const FLUID_KRYLOV_VECTOR& operator=(const FLUID_KRYLOV_VECTOR& bv);
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

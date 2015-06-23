//#####################################################################
// Copyright 2008.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PRESSURE_VELOCITY_VECTOR
//#####################################################################
#ifndef __PRESSURE_VELOCITY_VECTOR__
#define __PRESSURE_VELOCITY_VECTOR__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
namespace PhysBAM{
//#####################################################################
// Class PRESSURE_VELOCITY_VECTOR
//#####################################################################
template<class TV>
class PRESSURE_VELOCITY_VECTOR:public KRYLOV_VECTOR_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef KRYLOV_VECTOR_BASE<typename TV::SCALAR> BASE;
public:
    GENERALIZED_VELOCITY<TV>& solid_velocity;
    ARRAY<ARRAY<T> >& pressure;
    bool deep_copy;

    PRESSURE_VELOCITY_VECTOR(GENERALIZED_VELOCITY<TV>& solid_velocity_input,ARRAY<ARRAY<T> >& pressure_input);
    virtual ~PRESSURE_VELOCITY_VECTOR();

    PRESSURE_VELOCITY_VECTOR& operator=(const PRESSURE_VELOCITY_VECTOR& v);
    BASE& operator+=(const BASE& bv) override;
    BASE& operator-=(const BASE& bv) override;
    BASE& operator*=(const T a) override;
    void Copy(const T c,const BASE& bv) override;
    void Copy(const T c1,const BASE& bv1,const BASE& bv2) override;
    int Raw_Size() const override;
    T& Raw_Get(int i) override;
    KRYLOV_VECTOR_BASE<T>* Clone_Default() const override;
    void Resize(const KRYLOV_VECTOR_BASE<T>& v) override;
};
}
#endif

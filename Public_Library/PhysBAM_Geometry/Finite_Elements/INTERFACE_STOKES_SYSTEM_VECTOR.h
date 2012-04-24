//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERFACE_STOKES_SYSTEM_VECTOR
//#####################################################################
#ifndef __INTERFACE_STOKES_SYSTEM_VECTOR__
#define __INTERFACE_STOKES_SYSTEM_VECTOR__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
namespace PhysBAM{

template<class TV>
class INTERFACE_STOKES_SYSTEM_VECTOR:public KRYLOV_VECTOR_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef KRYLOV_VECTOR_BASE<T> BASE;
public:
    ARRAY<T> u[TV::m][2];
    ARRAY<T> q[TV::m][2];
    ARRAY<T> p[2];

    INTERFACE_STOKES_SYSTEM_VECTOR();
    ~INTERFACE_STOKES_SYSTEM_VECTOR();

    INTERFACE_STOKES_SYSTEM_VECTOR& operator=(const INTERFACE_STOKES_SYSTEM_VECTOR& v);
    BASE& operator+=(const BASE& bv) PHYSBAM_OVERRIDE;
    BASE& operator-=(const BASE& bv) PHYSBAM_OVERRIDE;
    BASE& operator*=(const T a) PHYSBAM_OVERRIDE;
    void Copy(const T c,const BASE& bv) PHYSBAM_OVERRIDE;
    void Copy(const T c1,const BASE& bv1,const BASE& bv2) PHYSBAM_OVERRIDE;
    void Print() const;
    int Raw_Size() const PHYSBAM_OVERRIDE;
    T& Raw_Get(int i) PHYSBAM_OVERRIDE;
    KRYLOV_VECTOR_BASE<T>* Clone_Default() const PHYSBAM_OVERRIDE;
    void Resize(const KRYLOV_VECTOR_BASE<T>& v) PHYSBAM_OVERRIDE;
};
}
#endif

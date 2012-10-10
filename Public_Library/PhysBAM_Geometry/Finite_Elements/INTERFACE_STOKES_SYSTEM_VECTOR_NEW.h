//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERFACE_STOKES_SYSTEM_VECTOR_NEW
//#####################################################################
#ifndef __INTERFACE_STOKES_SYSTEM_VECTOR_NEW__
#define __INTERFACE_STOKES_SYSTEM_VECTOR_NEW__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV>
class INTERFACE_STOKES_SYSTEM_VECTOR_NEW:public KRYLOV_VECTOR_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef KRYLOV_VECTOR_BASE<T> BASE;
public:
    VECTOR<VECTOR<ARRAY<T>,2>,TV::m> u;
    VECTOR<ARRAY<T>,2> p;
    ARRAY<T> q;

    INTERFACE_STOKES_SYSTEM_VECTOR_NEW();
    ~INTERFACE_STOKES_SYSTEM_VECTOR_NEW();

    INTERFACE_STOKES_SYSTEM_VECTOR_NEW& operator=(const INTERFACE_STOKES_SYSTEM_VECTOR_NEW& v);
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

    T Dot(const INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>& v) const;
    T Magnitude_Squared() const;
    T Magnitude() const;
    T Max_Abs() const;
    void Normalize();
    void Scale(const INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>& v,const INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>& c);
    void Scale(const INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>& c);
};
}
#endif

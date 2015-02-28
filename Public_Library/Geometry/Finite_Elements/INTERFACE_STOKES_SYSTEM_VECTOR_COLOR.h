//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERFACE_STOKES_SYSTEM_VECTOR_COLOR
//#####################################################################
#ifndef __INTERFACE_STOKES_SYSTEM_VECTOR_COLOR__
#define __INTERFACE_STOKES_SYSTEM_VECTOR_COLOR__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV>
class INTERFACE_STOKES_SYSTEM_VECTOR_COLOR:public KRYLOV_VECTOR_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef KRYLOV_VECTOR_BASE<T> BASE;
public:
    VECTOR<ARRAY<ARRAY<T> >,TV::m> u;
    ARRAY<ARRAY<T> > p;
    ARRAY<T> q;

    int colors;

    INTERFACE_STOKES_SYSTEM_VECTOR_COLOR();
    ~INTERFACE_STOKES_SYSTEM_VECTOR_COLOR();

    INTERFACE_STOKES_SYSTEM_VECTOR_COLOR& operator=(const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR& v);
    BASE& operator+=(const BASE& bv) PHYSBAM_OVERRIDE;
    BASE& operator-=(const BASE& bv) PHYSBAM_OVERRIDE;
    BASE& operator*=(const T a) PHYSBAM_OVERRIDE;
    void Copy(const T c1,const BASE& bv1) PHYSBAM_OVERRIDE;
    void Copy(const T c1,const BASE& bv1,const BASE& bv2) PHYSBAM_OVERRIDE;
    void Print() const;
    int Raw_Size() const PHYSBAM_OVERRIDE;
    T& Raw_Get(int i) PHYSBAM_OVERRIDE;
    KRYLOV_VECTOR_BASE<T>* Clone_Default() const PHYSBAM_OVERRIDE;
    void Resize(const KRYLOV_VECTOR_BASE<T>& v) PHYSBAM_OVERRIDE;

    T Max_Abs() const;
    void Scale(const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>& v,const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>& s);
};
}
#endif

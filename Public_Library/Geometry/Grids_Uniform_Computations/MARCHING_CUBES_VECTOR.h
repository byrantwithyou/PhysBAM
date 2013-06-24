//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################
// Class MARCHING_CUBES_VECTOR
//#####################################################################
#ifndef __MARCHING_CUBES_VECTOR__
#define __MARCHING_CUBES_VECTOR__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>

namespace PhysBAM{
//#####################################################################
// Class MARCHING_CUBES_VECTOR
//#####################################################################
template<class TV>
class MARCHING_CUBES_VECTOR:public KRYLOV_VECTOR_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef KRYLOV_VECTOR_BASE<T> BASE;

public:

    ARRAY<TV> x; 

    MARCHING_CUBES_VECTOR(){}
    virtual ~MARCHING_CUBES_VECTOR(){}

    BASE& operator+=(const BASE& bv) PHYSBAM_OVERRIDE
    {x+=debug_cast<const MARCHING_CUBES_VECTOR&>(bv).x;return *this;}

    BASE& operator-=(const BASE& bv) PHYSBAM_OVERRIDE
    {x-=debug_cast<const MARCHING_CUBES_VECTOR&>(bv).x;return *this;}

    BASE& operator*=(const T a) PHYSBAM_OVERRIDE
    {x*=a;return *this;}

    void Copy(const T c1,const BASE& bv1) PHYSBAM_OVERRIDE
    {x=debug_cast<const MARCHING_CUBES_VECTOR&>(bv1).x*c1;}

    void Copy(const T c1,const BASE& bv1,const BASE& bv2) PHYSBAM_OVERRIDE
    {x=debug_cast<const MARCHING_CUBES_VECTOR&>(bv1).x*c1+debug_cast<const MARCHING_CUBES_VECTOR&>(bv2).x;}

    T Dot(const KRYLOV_VECTOR_BASE<T>& bv) const PHYSBAM_OVERRIDE
    {return x.Dot(debug_cast<const MARCHING_CUBES_VECTOR&>(bv).x);}

    int Raw_Size() const PHYSBAM_OVERRIDE
    {return x.m*TV::m;}

    T& Raw_Get(int i) PHYSBAM_OVERRIDE
    {return x(i/TV::m)(i%TV::m);}

    KRYLOV_VECTOR_BASE<T>* Clone_Default() const PHYSBAM_OVERRIDE
    {MARCHING_CUBES_VECTOR* V=new MARCHING_CUBES_VECTOR;V->x.Resize(x.m);return V;}

    void Resize(const KRYLOV_VECTOR_BASE<T>& bv) PHYSBAM_OVERRIDE
    {x.Resize(debug_cast<const MARCHING_CUBES_VECTOR&>(bv).x.m);}
};
}
#endif

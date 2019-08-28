//#####################################################################
// Copyright 2009, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KRYLOV_VECTOR_BASE
//#####################################################################
#ifndef __KRYLOV_VECTOR_BASE__
#define __KRYLOV_VECTOR_BASE__

#include <Core/Arrays/ARRAY_VIEW.h>
#include <Core/Utilities/DEBUG_CAST.h>
#include <iostream>
namespace PhysBAM{

template<class T>
class KRYLOV_VECTOR_BASE
{
public:
    KRYLOV_VECTOR_BASE();
    virtual ~KRYLOV_VECTOR_BASE();

    KRYLOV_VECTOR_BASE& operator= (const KRYLOV_VECTOR_BASE& bv)
    {Copy((T)1,bv);return *this;}

    virtual KRYLOV_VECTOR_BASE& operator+=(const KRYLOV_VECTOR_BASE& bv)=0;
    virtual KRYLOV_VECTOR_BASE& operator-=(const KRYLOV_VECTOR_BASE& bv)=0;
    virtual KRYLOV_VECTOR_BASE& operator*=(const T a)=0;
    virtual void Copy(const T c,const KRYLOV_VECTOR_BASE& bv)=0;
    virtual void Copy(const T c1,const KRYLOV_VECTOR_BASE& bv1,const KRYLOV_VECTOR_BASE& bv2)=0;
    virtual int Raw_Size() const=0;
    virtual T& Raw_Get(int i)=0;
    const T& Raw_Get(int i) const;
    virtual void Get(ARRAY_VIEW<T> a) const=0;
    virtual void Set(ARRAY_VIEW<const T> a)=0;
    virtual KRYLOV_VECTOR_BASE* Clone_Default() const=0;
    virtual void Resize(const KRYLOV_VECTOR_BASE& v)=0;
//#####################################################################
};
template<class T> std::ostream& operator<<(std::ostream& output,const KRYLOV_VECTOR_BASE<T>& x);
}
#endif

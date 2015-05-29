//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OLDROYD_CONSTITUTIVE_MODEL
//#####################################################################
#ifndef __OLDROYD_CONSTITUTIVE_MODEL__
#define __OLDROYD_CONSTITUTIVE_MODEL__

#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{

template<class TV>
class OLDROYD_CONSTITUTIVE_MODEL:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:

    OLDROYD_CONSTITUTIVE_MODEL(){}
    virtual ~OLDROYD_CONSTITUTIVE_MODEL(){}

//#####################################################################
    virtual void Resize(int n)=0;
    virtual void Precompute(const MATRIX<T,TV::m>& F,const SYMMETRIC_MATRIX<T,TV::m>& S,int p)=0;
    virtual T Energy_Density(int p) const=0;
    virtual void Gradient(MATRIX<T,TV::m>& dF,SYMMETRIC_MATRIX<T,TV::m>& dS,int p) const=0;
    virtual void Hessian(const MATRIX<T,TV::m>& F,const SYMMETRIC_MATRIX<T,TV::m>& S,
        MATRIX<T,TV::m>& dF,SYMMETRIC_MATRIX<T,TV::m>& dS,int p) const=0;
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OLDROYD_NEO_HOOKEAN
//#####################################################################
#ifndef __OLDROYD_NEO_HOOKEAN__
#define __OLDROYD_NEO_HOOKEAN__

#include <Core/Arrays/ARRAY.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Hybrid_Methods/Forces/OLDROYD_CONSTITUTIVE_MODEL.h>
namespace PhysBAM{

template<class TV>
class OLDROYD_NEO_HOOKEAN:public OLDROYD_CONSTITUTIVE_MODEL<TV>
{
    typedef typename TV::SCALAR T;
public:

    T mu,lambda;

    ARRAY<T> psi,b,c;
    ARRAY<MATRIX<T,TV::m> > P,H;
    ARRAY<SYMMETRIC_MATRIX<T,TV::m> > Q;

    OLDROYD_NEO_HOOKEAN();
    virtual ~OLDROYD_NEO_HOOKEAN();

//#####################################################################
    void Resize(int n) override;
    void Precompute(const MATRIX<T,TV::m>& F,const SYMMETRIC_MATRIX<T,TV::m>& S,int p) override;
    T Energy_Density(int p) const override;
    void Gradient(MATRIX<T,TV::m>& dF,SYMMETRIC_MATRIX<T,TV::m>& dS,int p) const override;
    void Hessian(const MATRIX<T,TV::m>& F,const SYMMETRIC_MATRIX<T,TV::m>& S,
        MATRIX<T,TV::m>& dF,SYMMETRIC_MATRIX<T,TV::m>& dS,int p) const override;
//#####################################################################
};
}
#endif

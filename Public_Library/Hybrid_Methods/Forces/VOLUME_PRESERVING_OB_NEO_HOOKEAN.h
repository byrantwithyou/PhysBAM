//#####################################################################
// Copyright 2015
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VOLUME_PRESERVING_OB_NEO_HOOKEAN
//#####################################################################
#ifndef __VOLUME_PRESERVING_OB_NEO_HOOKEAN__
#define __VOLUME_PRESERVING_OB_NEO_HOOKEAN__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/IDENTITY_MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Hybrid_Methods/Forces/OLDROYD_CONSTITUTIVE_MODEL.h>
namespace PhysBAM{

template<class TV>
class VOLUME_PRESERVING_OB_NEO_HOOKEAN:public OLDROYD_CONSTITUTIVE_MODEL<TV>
{
    typedef typename TV::SCALAR T;
public:

    T mu,lambda;

    ARRAY<T> psi,b,c;
    ARRAY<MATRIX<T,TV::m> > P,H,G,F_mat;
    ARRAY<SYMMETRIC_MATRIX<T,TV::m> > Q,N,M,S_mat;

    VOLUME_PRESERVING_OB_NEO_HOOKEAN();
    virtual ~VOLUME_PRESERVING_OB_NEO_HOOKEAN();

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

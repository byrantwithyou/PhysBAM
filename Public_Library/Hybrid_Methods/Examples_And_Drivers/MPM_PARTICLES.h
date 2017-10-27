//#####################################################################
// Copyright 2004-2009, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_PARTICLES
//#####################################################################
#ifndef __MPM_PARTICLES__
#define __MPM_PARTICLES__

#include <Core/Arrays/ARRAY.h>
#include <Core/Matrices/MATRIX.h>
#include <Tools/Particles/PARTICLES.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Hybrid_Methods/Examples_And_Drivers/PHASE_ID.h>

namespace PhysBAM{

template<class TV>
class MPM_PARTICLES:public CLONEABLE<MPM_PARTICLES<TV>,DEFORMABLE_PARTICLES<TV> >
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<MPM_PARTICLES<TV>,DEFORMABLE_PARTICLES<TV> > BASE;
public:
    using BASE::Add_Array;using BASE::Remove_Array;

    bool store_Fp,store_B,store_C,store_S,store_lame,store_lame0,store_plastic_def,store_phase,store_vort;
    ARRAY_VIEW<T> volume;
    ARRAY_VIEW<MATRIX<T,TV::m> > F,Fp,B,C;
    ARRAY_VIEW<SYMMETRIC_MATRIX<T,TV::m> > S;
    ARRAY_VIEW<bool> valid;
    ARRAY_VIEW<T> mu,lambda,mu0,lambda0;
    ARRAY_VIEW<PHASE_ID> phase;
    ARRAY_VIEW<T> vort;

    MPM_PARTICLES();
    virtual ~MPM_PARTICLES();
    void Store_Fp(bool store);
    void Store_B(bool store);
    void Store_C(bool store);
    void Store_S(bool store);
    void Store_Lame(bool store);
    void Store_Lame0(bool store);
    void Store_Phase(bool store);
    void Store_Vort(bool store);

//#####################################################################
};
}
#endif

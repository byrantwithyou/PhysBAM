//#####################################################################
// Copyright 2004-2009, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_PARTICLES
//#####################################################################
#ifndef __MPM_PARTICLES__
#define __MPM_PARTICLES__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Particles/PARTICLES.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>

namespace PhysBAM{

const ATTRIBUTE_ID ATTRIBUTE_ID_F(40);
const ATTRIBUTE_ID ATTRIBUTE_ID_VOLUME(41);
const ATTRIBUTE_ID ATTRIBUTE_ID_B(42);
const ATTRIBUTE_ID ATTRIBUTE_ID_VALID(43);
const ATTRIBUTE_ID ATTRIBUTE_ID_S(44);
const ATTRIBUTE_ID ATTRIBUTE_ID_C(45);
const ATTRIBUTE_ID ATTRIBUTE_ID_ONE_OVER_LAMBDA(46);

template<class TV>
class MPM_PARTICLES:public CLONEABLE<MPM_PARTICLES<TV>,DEFORMABLE_PARTICLES<TV> >
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<MPM_PARTICLES<TV>,DEFORMABLE_PARTICLES<TV> > BASE;
public:
    using BASE::Add_Array;using BASE::Remove_Array;

    bool store_B,store_S,store_C,store_one_over_lambda;
    ARRAY_VIEW<T> volume;
    ARRAY_VIEW<MATRIX<T,TV::m> > F,B,C;
    ARRAY_VIEW<SYMMETRIC_MATRIX<T,TV::m> > S;
    ARRAY_VIEW<bool> valid;
    ARRAY_VIEW<T> one_over_lambda;

    MPM_PARTICLES();
    virtual ~MPM_PARTICLES();
    void Store_B(bool store=true);
    void Store_S(bool store=true);
    void Store_C(bool store=true);
    void Store_One_Over_Lambda(bool store=true);

//#####################################################################
};
}
#endif

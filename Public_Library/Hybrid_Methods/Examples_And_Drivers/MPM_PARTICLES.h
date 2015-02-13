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

namespace PhysBAM{

const ATTRIBUTE_ID ATTRIBUTE_ID_F(40);
const ATTRIBUTE_ID ATTRIBUTE_ID_VOLUME(41);
const ATTRIBUTE_ID ATTRIBUTE_ID_B(42);

template<class TV>
class MPM_PARTICLES:public CLONEABLE<MPM_PARTICLES<TV>,PARTICLES<TV> >
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<MPM_PARTICLES<TV>,PARTICLES<TV> > BASE;
public:
    using BASE::Add_Array;

    ARRAY_VIEW<TV> X,V;
    ARRAY_VIEW<T> mass;
    ARRAY_VIEW<T> volume;
    ARRAY_VIEW<MATRIX<T,TV::m> > F,B;

    MPM_PARTICLES();
    virtual ~MPM_PARTICLES();
//#####################################################################
};
}
#endif

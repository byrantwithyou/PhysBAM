//#####################################################################
// Copyright 2007, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_EIGENSYSTEM_BASE
//##################################################################### 
#ifndef __EULER_EIGENSYSTEM_BASE__
#define __EULER_EIGENSYSTEM_BASE__   

#include <Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
#include <Compressible/Euler_Equations/EULER.h>
namespace PhysBAM{

template<class TV>
class EULER_EIGENSYSTEM_BASE:public EIGENSYSTEM<typename TV::SCALAR,VECTOR<typename TV::SCALAR,TV::m+2> >
{
    typedef typename TV::SCALAR T;typedef VECTOR<T,TV::m+2> TV_DIMENSION;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<TV_DIMENSION,TV_INT> T_ARRAYS_DIMENSION_SCALAR;
    typedef EIGENSYSTEM<T,TV_DIMENSION> EIGENSYSTEM_BASE;
public:
    using EIGENSYSTEM_BASE::Flux;using EIGENSYSTEM_BASE::Eigenvalues;using EIGENSYSTEM_BASE::Eigenvectors;using EIGENSYSTEM_BASE::slice_index;
    EOS<T>* eos;

    EULER_EIGENSYSTEM_BASE(EOS<T>* eos)
        :eos(eos)
    {}
};
}
#endif


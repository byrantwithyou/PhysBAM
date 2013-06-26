//#####################################################################
// Copyright 2007, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_EIGENSYSTEM
//##################################################################### 
#ifndef __EULER_EIGENSYSTEM__
#define __EULER_EIGENSYSTEM__   

#include <Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
#include <Compressible/Euler_Equations/EULER.h>
namespace PhysBAM{

template<class T_GRID>
class EULER_EIGENSYSTEM:public EIGENSYSTEM<typename T_GRID::SCALAR,VECTOR<typename T_GRID::SCALAR,T_GRID::dimension+2> >,public EULER<T_GRID>
{
    typedef typename T_GRID::SCALAR T;typedef VECTOR<T,T_GRID::dimension+2> TV_DIMENSION;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef EULER<T_GRID> EULER_BASE;
    typedef EIGENSYSTEM<T,TV_DIMENSION> EIGENSYSTEM_BASE;
public:
    using EULER_BASE::eos;using EIGENSYSTEM_BASE::Flux;using EIGENSYSTEM_BASE::Eigenvalues;using EIGENSYSTEM_BASE::Eigenvectors;using EIGENSYSTEM_BASE::slice_index;

    EULER_EIGENSYSTEM()
    {}
};
}
#endif


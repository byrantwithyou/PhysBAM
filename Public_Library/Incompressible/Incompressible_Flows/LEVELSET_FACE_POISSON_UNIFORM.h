//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_FACE_POISSON_UNIFORM
//#####################################################################
#ifndef __LEVELSET_FACE_POISSON_UNIFORM__
#define __LEVELSET_FACE_POISSON_UNIFORM__

#include <Core/Arrays/ARRAY.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
namespace PhysBAM{

template<class TV> class BOUNDARY_CONDITIONS_CALLBACKS;
template<class TV> class LEVELSET_INDEX_MAP_UNIFORM;

template<class TV>
class LEVELSET_FACE_POISSON_UNIFORM
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
public:

    LEVELSET_INDEX_MAP_UNIFORM<TV>& index_map;
    SPARSE_MATRIX_FLAT_MXN<T> P;
    ARRAY<T> b;
    T theta_threshold;

    LEVELSET_FACE_POISSON_UNIFORM(LEVELSET_INDEX_MAP_UNIFORM<TV>& index_map_input);
    LEVELSET_FACE_POISSON_UNIFORM(const LEVELSET_FACE_POISSON_UNIFORM&) = delete;
    void operator=(const LEVELSET_FACE_POISSON_UNIFORM&) = delete;
    virtual ~LEVELSET_FACE_POISSON_UNIFORM();

//#####################################################################
    void Compute(int axis);
//#####################################################################
};
}
#endif

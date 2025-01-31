//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_FLUID_GRADIENT_CUT
//#####################################################################
#ifndef __MATRIX_FLUID_GRADIENT_CUT__
#define __MATRIX_FLUID_GRADIENT_CUT__
#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Dynamics/Coupled_Evolution/MATRIX_FLUID_GRADIENT_BASE.h>

namespace PhysBAM{

template<class TV> class COLLISION_AWARE_INDEX_MAP;
template<class TV> class GRID;

template<class TV>
class MATRIX_FLUID_GRADIENT_CUT:public MATRIX_FLUID_GRADIENT_BASE<TV>
{
    enum WORKAROUND {d=TV::m};
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,d> TV_INT;
    typedef MATRIX_FLUID_GRADIENT_BASE<TV> BASE;

public:
    using BASE::gradient;using BASE::index_map;using BASE::ghost_gradient;

    // TODO: the only reason index_map isn't const is because Array_View doesn't work for const quite yet
    MATRIX_FLUID_GRADIENT_CUT(COLLISION_AWARE_INDEX_MAP<TV>& index_map_input);
    virtual ~MATRIX_FLUID_GRADIENT_CUT();

//#####################################################################
    void Compute(const ARRAY<bool,FACE_INDEX<d> >& psi_N_domain_boundary) override;
//#####################################################################
};
}
#endif

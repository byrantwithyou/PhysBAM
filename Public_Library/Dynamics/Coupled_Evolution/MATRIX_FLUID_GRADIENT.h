//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_FLUID_GRADIENT
//#####################################################################
#ifndef __MATRIX_FLUID_GRADIENT__
#define __MATRIX_FLUID_GRADIENT__
#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Dynamics/Coupled_Evolution/MATRIX_FLUID_GRADIENT_BASE.h>

namespace PhysBAM{

template<class TV> class COLLISION_AWARE_INDEX_MAP;
template<class TV> class GRID;

template<class TV>
class MATRIX_FLUID_GRADIENT:public MATRIX_FLUID_GRADIENT_BASE<TV>
{
    enum WORKAROUND {d=TV::m};
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,d> TV_INT;
    typedef MATRIX_FLUID_GRADIENT_BASE<TV> BASE;

public:
    using BASE::gradient;using BASE::index_map;using BASE::ghost_gradient;using BASE::interface_gradient;
    VECTOR<bool,d> periodic_boundary;

    // TODO: the only reason index_map isn't const is because Array_View doesn't work for const quite yet
    MATRIX_FLUID_GRADIENT(COLLISION_AWARE_INDEX_MAP<TV>& index_map_input);
    virtual ~MATRIX_FLUID_GRADIENT();

//#####################################################################
    void Compute(const ARRAY<bool,FACE_INDEX<d> >& psi_N_domain_boundary);
    void Add_Cell(int face_lookup,int axis,const TV_INT& cell_index,T weight);
    void Add_Interface(int face_lookup,const FACE_INDEX<d>& face_index,T weight);
    void Add_Interface(int face_lookup,const SIDED_FACE_INDEX<d>& face_index,T weight);
//#####################################################################
};
}
#endif

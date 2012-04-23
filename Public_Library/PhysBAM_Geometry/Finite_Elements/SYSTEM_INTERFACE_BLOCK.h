//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYSTEM_INTERFACE_BLOCK 
//#####################################################################

#ifndef __SYSTEM_INTERFACE_BLOCK__
#define __SYSTEM_INTERFACE_BLOCK__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Symbolics/STATIC_POLYNOMIAL.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER.h>

namespace PhysBAM{

template<class TV,int static_degree>
class SYSTEM_INTERFACE_BLOCK
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

public:

    struct OVERLAP_POLYNOMIAL
    {
        TV_INT index_offset;
        int subcell; // flags indicating fine cells
        STATIC_POLYNOMIAL<T,TV::m,static_degree> polynomial;
        ARRAY<int,TV_INT> flat_diff_index;
    };

    T scale;
    bool ignore_orientation;
    CELL_MANAGER<TV> *cm;
    CELL_DOMAIN_INTERFACE<TV> *cdi;
    ARRAY<OVERLAP_POLYNOMIAL> overlap_polynomials;
    ARRAY<int> flat_diff;
    MATRIX_MXN<T> data[2]; // inside and outside

    inline void Add_Entry(int interface_element,int flat_diff_index,int inside,T value)
    {
        data[inside](interface_element,flat_diff_index)+=value*scale;
        cm->Set_Active(cdi->Get_Flat_Base(interface_element)+flat_diff(flat_diff_index),inside);
    }

    template<int d>
    void Initialize(const BASIS_STENCIL_UNIFORM<TV,d>& s,CELL_MANAGER<TV>& cm_input,CELL_DOMAIN_INTERFACE<TV>& cdi,
        T scale_input,bool ignore_orientation_input);
};
}
#endif

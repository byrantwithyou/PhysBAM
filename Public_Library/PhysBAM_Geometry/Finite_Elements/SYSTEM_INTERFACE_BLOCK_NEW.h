//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYSTEM_INTERFACE_BLOCK_NEW 
//#####################################################################

#ifndef __SYSTEM_INTERFACE_BLOCK_NEW__
#define __SYSTEM_INTERFACE_BLOCK_NEW__

#include <PhysBAM_Tools/Symbolics/STATIC_POLYNOMIAL.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_INTERFACE_BLOCK_HELPER_NEW.h>

namespace PhysBAM{

template<class TV,int static_degree>
class SYSTEM_INTERFACE_BLOCK_NEW:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    CELL_DOMAIN_INTERFACE_NEW<TV>* cdi;
    SYSTEM_INTERFACE_BLOCK_HELPER_NEW<TV> *helper;

public:

    struct OVERLAP_POLYNOMIAL
    {
        int flat_index_offset;
        int flat_index_diff_ref;
        int subcell; // flags indicating fine cells
        STATIC_POLYNOMIAL<T,TV::m,static_degree> polynomial;
    };

    
    T scale;
    int axis;
    ARRAY<OVERLAP_POLYNOMIAL> overlap_polynomials;

    template<int d>
    void Initialize(SYSTEM_INTERFACE_BLOCK_HELPER_NEW<TV>& helper_input,const BASIS_STENCIL_UNIFORM<TV,d>& s,int axis_input,T scale_input);

    void Add_Entry(int interface_dof,int orientation,int flat_index_diff_ref,int inside,T value)
    {helper->data[inside](orientation*cdi->interface_dofs+interface_dof,flat_index_diff_ref)+=value*scale;}

    int Flat_Diff(int i)
    {return helper->flat_diff(i);}
};
}
#endif

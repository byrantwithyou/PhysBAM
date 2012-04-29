//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYSTEM_INTERFACE_BLOCK 
//#####################################################################

#ifndef __SYSTEM_INTERFACE_BLOCK__
#define __SYSTEM_INTERFACE_BLOCK__

#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_INTERFACE_BLOCK_HELPER.h>

namespace PhysBAM{

template<class TV,int static_degree>
class SYSTEM_INTERFACE_BLOCK:NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    CELL_DOMAIN_INTERFACE<TV>* cdi;
    SYSTEM_INTERFACE_BLOCK_HELPER<TV> *helper;

public:

    struct OVERLAP_POLYNOMIAL
    {
        int flat_index_offset;
        ARRAY<int,TV_INT> flat_index_diff;
        int subcell; // flags indicating fine cells
        STATIC_POLYNOMIAL<T,TV::m,static_degree> polynomial;
    };

    T scale;
    bool ignore_orientation;
    ARRAY<OVERLAP_POLYNOMIAL> overlap_polynomials;

    template<int d>
    void Initialize(SYSTEM_INTERFACE_BLOCK_HELPER<TV>& helper_input,const BASIS_STENCIL_UNIFORM<TV,d>& s,T scale_input,bool ignore_orientation_input)
    {
        scale=scale_input;
        ignore_orientation=ignore_orientation_input;
        helper=&helper_input;
        cdi=helper->cdi;
        
        overlap_polynomials.Resize(s.diced.m);
        for(int i=0;i<overlap_polynomials.m;i++){
            OVERLAP_POLYNOMIAL& op=overlap_polynomials(i);
            const typename BASIS_STENCIL_UNIFORM<TV,d>::DICED& diced=s.diced(i);
            op.flat_index_offset=cdi->Flatten_Diff(diced.index_offset);
            op.flat_index_diff.Resize(cdi->coarse_range);
            for(RANGE_ITERATOR<TV::m> it(cdi->coarse_range);it.Valid();it.Next())
                op.flat_index_diff(it.index)=helper->flat_diff.Binary_Search(cdi->Flatten_Diff(it.index)+op.flat_index_offset);
            op.subcell=diced.subcell;
            op.polynomial=diced.polynomial;}
    }

    inline void Add_Entry(int interface_element,int flat_index_diff,int inside,T value)
    {helper->data[inside](interface_element,flat_index_diff)+=value*scale;}
};
}
#endif

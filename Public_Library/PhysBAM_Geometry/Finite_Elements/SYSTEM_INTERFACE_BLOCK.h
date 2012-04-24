//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYSTEM_INTERFACE_BLOCK 
//#####################################################################

#ifndef __SYSTEM_INTERFACE_BLOCK__
#define __SYSTEM_INTERFACE_BLOCK__

#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_INTERFACE_BLOCK_HELPER.h>

namespace PhysBAM{

template<class TV,int static_degree>
class SYSTEM_INTERFACE_BLOCK
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    CELL_DOMAIN_INTERFACE<TV>* cdi;
    SYSTEM_INTERFACE_BLOCK_HELPER<TV> *helper;

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
    ARRAY<OVERLAP_POLYNOMIAL> overlap_polynomials;

    template<int d>
    void Initialize(SYSTEM_INTERFACE_BLOCK_HELPER<TV> helper_input,const BASIS_STENCIL_UNIFORM<TV,d>& s,T scale_input,bool ignore_orientation_input)
    {
        scale=scale_input;
        ignore_orientation=ignore_orientation_input;
        helper=&helper_input;
        cdi=helper->cdi;
        
        overlap_polynomials.Resize(s.diced.m);
        for(int i=0;i<overlap_polynomials.m;i++){
            OVERLAP_POLYNOMIAL& op=overlap_polynomials(i);
            const typename BASIS_STENCIL_UNIFORM<TV,d>::DICED& diced=s.diced(i);
            op.index_offset=diced.index_offset;
            op.subcell=diced.subcell;
            op.polynomial=diced.polynomial;
            op.flat_diff_index.Resize(cdi->coarse_range);
            for(RANGE_ITERATOR<TV::m> it(cdi->coarse_range);it.Valid();it.Next())
                op.flat_diff_index(it.index)=helper->flat_diff.Binary_Search(cdi->Flatten_Diff(it.index+op.index_offset));}
    }

    inline void Add_Entry(int interface_element,int flat_diff_index,int inside,T value)
    {helper->data[inside](interface_element,flat_diff_index)+=value*scale;}
};
}
#endif

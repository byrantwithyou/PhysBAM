//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYSTEM_VOLUME_BLOCK
//#####################################################################

#ifndef __SYSTEM_VOLUME_BLOCK__
#define __SYSTEM_VOLUME_BLOCK__

#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK_HELPER.h>

namespace PhysBAM{

template<class TV,int static_degree>
class SYSTEM_VOLUME_BLOCK
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    CELL_DOMAIN_INTERFACE<TV>* cdi;
    SYSTEM_VOLUME_BLOCK_HELPER<TV> *helper;

public:
    
    struct OPEN_ENTRY
    {
        TV_INT index0, index1;
        T x;
        int flat_diff_index;
        
        bool operator< (const OPEN_ENTRY& me) const
        {
            assert(flat_diff_index==me.flat_diff_index);
            if(index0!=me.index0) return LEXICOGRAPHIC_COMPARE()(index0,me.index0);
            return LEXICOGRAPHIC_COMPARE()(index1,me.index1);
        }
        
        void Merge(const OPEN_ENTRY& me){x+=me.x;}
    };
    
    struct OVERLAP_POLYNOMIAL
    {
        TV_INT index_offset0,index_offset1;
        int subcell; // flags indicating fine cells
        STATIC_POLYNOMIAL<T,TV::m,static_degree> polynomial;
        int flat_diff_index;
    };

    VECTOR<T,2> scale;
    ARRAY<OVERLAP_POLYNOMIAL> overlap_polynomials;
    ARRAY<OPEN_ENTRY> open_entries,open_subcell_entries[1<<TV::m];

    template<int d0,int d1>
    void Initialize(SYSTEM_VOLUME_BLOCK_HELPER<TV>& helper_input,const BASIS_STENCIL_UNIFORM<TV,d0>& s0,
        const BASIS_STENCIL_UNIFORM<TV,d1>& s1,const VECTOR<T,2>& scale_input)
    {
        scale=scale_input;
        helper=&helper_input;
        cdi=helper->cdi;
        
        for(int i=0;i<s0.diced.m;i++)
            for(int j=0;j<s1.diced.m;j++){
                int overlap=s0.diced(i).subcell&s1.diced(j).subcell;
                if(overlap){
                    const typename BASIS_STENCIL_UNIFORM<TV,d0>::DICED& diced0=s0.diced(i);
                    const typename BASIS_STENCIL_UNIFORM<TV,d1>::DICED& diced1=s1.diced(j);
                    OVERLAP_POLYNOMIAL op={diced0.index_offset,diced1.index_offset,overlap};
                    op.polynomial=diced0.polynomial*diced1.polynomial;
                    op.flat_diff_index=helper->flat_diff.Binary_Search(cdi->Flatten_Diff(op.index_offset1-op.index_offset0));
                    overlap_polynomials.Append(op);}}
    }

    inline void Add_Entry(const TV_INT& index,int flat_diff_index,int inside,T value)
    {helper->data[inside](cdi->Flatten(index),flat_diff_index)+=value*scale(inside);}
    inline void Add_Open_Entry(const TV_INT& cell,int inside,OPEN_ENTRY& oe)
    {Add_Entry(oe.index0+cell,oe.flat_diff_index,inside,oe.x);}

    void Add_Open_Entries(const TV_INT& cell,int inside)
    {for(int j=0;j<open_entries.m;j++) Add_Open_Entry(cell,inside,open_entries(j));}

    void Add_Open_Subcell_Entries(const TV_INT& cell,int block,int inside)
    {for(int j=0;j<open_entries.m;j++) Add_Open_Entry(cell,block,open_subcell_entries[block](j));}
};
}
#endif


//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYSTEM_VOLUME_BLOCK
//#####################################################################
#ifndef __SYSTEM_VOLUME_BLOCK__
#define __SYSTEM_VOLUME_BLOCK__

#include <PhysBAM_Tools/Symbolics/STATIC_POLYNOMIAL.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK_HELPER.h>

namespace PhysBAM{

template<class TV,int static_degree>
class SYSTEM_VOLUME_BLOCK:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;

    CELL_DOMAIN_INTERFACE<TV>* cdi;
    SYSTEM_VOLUME_BLOCK_HELPER<TV> *helper;

public:
    
    struct OPEN_ENTRY
    {
        int flat_index_offset;
        int flat_index_diff;
        T x;
        
        bool operator< (const OPEN_ENTRY& me) const
        {
            if(flat_index_offset!=me.flat_index_offset) return flat_index_offset<me.flat_index_offset; 
            return flat_index_diff<me.flat_index_diff; 
        }
        
        void Merge(const OPEN_ENTRY& me){x+=me.x;}
    };
    
    struct OVERLAP_POLYNOMIAL
    {
        int flat_index_offset;
        int flat_index_diff;
        int subcell; // flags indicating fine cells
        STATIC_POLYNOMIAL<T,TV::m,static_degree> polynomial;
    };

    VECTOR<T,2> scale;
    ARRAY<OVERLAP_POLYNOMIAL> overlap_polynomials;
    ARRAY<OPEN_ENTRY> open_entries,open_subcell_entries[1<<TV::m];

    void Add_Entry(int flat_index,int flat_index_diff,int inside,T value)
    {helper->data[inside](flat_index,flat_index_diff)+=value*scale(inside);}

    void Add_Open_Entry(int flat_index,int inside,OPEN_ENTRY& oe)
    {Add_Entry(flat_index+oe.flat_index_offset,oe.flat_index_diff,inside,oe.x);}

    template<int d0,int d1>
    void Initialize(SYSTEM_VOLUME_BLOCK_HELPER<TV>& helper_input,const BASIS_STENCIL_UNIFORM<TV,d0>& s0,
        const BASIS_STENCIL_UNIFORM<TV,d1>& s1,const VECTOR<T,2>& scale_input);
    void Add_Open_Entries(int flat_index,int inside);
    void Add_Open_Subcell_Entries(int flat_index,int block,int inside);
};
}
#endif


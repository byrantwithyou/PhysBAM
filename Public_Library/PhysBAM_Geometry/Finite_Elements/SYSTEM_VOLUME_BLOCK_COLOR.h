//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYSTEM_VOLUME_BLOCK_COLOR
//#####################################################################
#ifndef __SYSTEM_VOLUME_BLOCK_COLOR__
#define __SYSTEM_VOLUME_BLOCK_COLOR__

#include <PhysBAM_Tools/Symbolics/STATIC_POLYNOMIAL.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK_HELPER_COLOR.h>

namespace PhysBAM{

template<class TV,int static_degree>
class SYSTEM_VOLUME_BLOCK_COLOR:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;

    SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV> *helper;

public:

    struct OPEN_ENTRY
    {
        int flat_index_offset;
        int flat_index_diff_ref;
        T x;
        
        bool operator< (const OPEN_ENTRY& oe) const
        {
            if(flat_index_offset!=oe.flat_index_offset) return flat_index_offset<oe.flat_index_offset; 
            return flat_index_diff_ref<oe.flat_index_diff_ref; 
        }
        
        void Merge(const OPEN_ENTRY& oe){x+=oe.x;}
    };
    
    struct OVERLAP_POLYNOMIAL
    {
        int flat_index_offset;
        int flat_index_diff_ref;
        int subcell; // flags indicating fine cells
        STATIC_POLYNOMIAL<T,TV::m,static_degree> polynomial;
    };

    ARRAY<T> scale;
    ARRAY<OVERLAP_POLYNOMIAL> overlap_polynomials;
    ARRAY<OPEN_ENTRY> open_entries,open_subcell_entries[1<<TV::m];

    void Add_Entry(int flat_index,int flat_index_diff_ref,int color,T value)
    {helper->data(color)(flat_index,flat_index_diff_ref)+=value*scale(color);}

    void Add_Open_Entry(int flat_index,int color,OPEN_ENTRY& oe)
    {Add_Entry(flat_index+oe.flat_index_offset,oe.flat_index_diff_ref,color,oe.x);}

    template<int d0,int d1>
    void Initialize(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>& helper_input,const BASIS_STENCIL_UNIFORM<TV,d0>& s0,
        const BASIS_STENCIL_UNIFORM<TV,d1>& s1,const ARRAY<T>& scale_input);
    void Add_Open_Entries(int flat_index,int color);
    void Add_Open_Subcell_Entries(int flat_index,int subcell,int color);
};
}
#endif


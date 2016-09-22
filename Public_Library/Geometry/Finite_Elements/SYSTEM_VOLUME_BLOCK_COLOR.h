//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYSTEM_VOLUME_BLOCK_COLOR
//#####################################################################
#ifndef __SYSTEM_VOLUME_BLOCK_COLOR__
#define __SYSTEM_VOLUME_BLOCK_COLOR__

#include <Core/Utilities/NONCOPYABLE.h>
#include <Tools/Symbolics/STATIC_POLYNOMIAL.h>
#include <Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK_HELPER_COLOR.h>

namespace PhysBAM{

template<class TV,int static_degree>
class SYSTEM_VOLUME_BLOCK_COLOR:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
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

    void Add_Entry(int flat_index,int flat_index_diff_ref,int color,T value);
    void Add_Open_Entry(int flat_index,int color,OPEN_ENTRY& oe);
    template<int d0,class ...Args>
    void Initialize(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>& helper_input,const ARRAY<T>& scale_input,
        const BASIS_STENCIL_UNIFORM<TV,d0>& s0,Args&& ...args);
    template<int d1,int pd,class ...Args>
    void Initialize_Helper(int subcell,int flat_index_offset,const TV_INT& index_offset,
        const STATIC_POLYNOMIAL<T,TV::m,pd>& poly,ARRAY<int>& diffs,const BASIS_STENCIL_UNIFORM<TV,d1>& s1,
        Args&& ...args);
    template<int pd>
    void Initialize_Helper(int subcell,int flat_index_offset,const TV_INT& index_offset,
        const STATIC_POLYNOMIAL<T,TV::m,pd>& poly,ARRAY<int>& diffs);
    void Add_Open_Entries(int flat_index,int color);
    void Add_Open_Subcell_Entries(int flat_index,int subcell,int color);
};
}
#endif


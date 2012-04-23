//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYSTEM_VOLUME_BLOCK
//#####################################################################

#ifndef __SYSTEM_VOLUME_BLOCK__
#define __SYSTEM_VOLUME_BLOCK__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Symbolics/STATIC_POLYNOMIAL.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER.h>

namespace PhysBAM{

template<class TV,int static_degree>
class SYSTEM_VOLUME_BLOCK
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

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
    CELL_MANAGER<TV> *cm0,*cm1;
    CELL_DOMAIN_INTERFACE<TV> *cdi;
    ARRAY<OVERLAP_POLYNOMIAL> overlap_polynomials;
    ARRAY<OPEN_ENTRY> open_entries,open_subcell_entries[1<<TV::m];
    ARRAY<int> flat_diff;
    MATRIX_MXN<T> data[2];  // inside and outside

    template<int d0,int d1>
    void Initialize(const BASIS_STENCIL_UNIFORM<TV,d0>& s0,const BASIS_STENCIL_UNIFORM<TV,d1>& s1,
        CELL_MANAGER<TV>& cm0_input,CELL_MANAGER<TV>& cm1_input,CELL_DOMAIN_INTERFACE<TV> &cdi,const VECTOR<T,2>& scale_input);

    inline void Add_Entry(const TV_INT& index,int flat_diff_index,int inside,T value)
    {data[inside](cdi->Flatten(index),flat_diff_index)+=value*scale(inside);}
    inline void Add_Open_Entry(const TV_INT& cell,int inside,OPEN_ENTRY& oe)
    {Add_Entry(oe.index0+cell,oe.flat_diff_index,inside,oe.x);}

    void Add_Open_Entries(const TV_INT& cell,int inside);
    void Add_Open_Subcell_Entries(const TV_INT& cell,int block,int inside);
    void Mark_Active_Cells(T tol=0);
};
}
#endif

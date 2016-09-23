//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYSTEM_VOLUME_BLOCK_HELPER_COLOR
//#####################################################################

#ifndef __SYSTEM_VOLUME_BLOCK_HELPER_COLOR__
#define __SYSTEM_VOLUME_BLOCK_HELPER_COLOR__

#include <Core/Arrays/ARRAY.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include <Core/Utilities/NONCOPYABLE.h>

namespace PhysBAM{

template<class TV> class CELL_DOMAIN_INTERFACE_COLOR;
template<class TV> class CELL_MANAGER_COLOR;
template<class TV,int d> struct BASIS_STENCIL_UNIFORM;
template<class T> class SPARSE_MATRIX_FLAT_MXN;

template<class TV>
class SYSTEM_VOLUME_BLOCK_HELPER_COLOR:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

public:
    ARRAY<MATRIX_MXN<T> > data;
    CELL_DOMAIN_INTERFACE_COLOR<TV>* cdi;
    ARRAY<CELL_MANAGER_COLOR<TV>*> cm; // size = # basis functions
    ARRAY<ARRAY<int> > flat_diff;

    template<int d0,class ...Args> 
    void Initialize(CELL_DOMAIN_INTERFACE_COLOR<TV> &cdi_input,const BASIS_STENCIL_UNIFORM<TV,d0>& s0,CELL_MANAGER_COLOR<TV>& cm0,Args&& ...args);
    template<int d1,class ...Args> 
    void Initialize_Helper(int subcell,ARRAY<int>& diffs,const TV_INT& index_offset,const BASIS_STENCIL_UNIFORM<TV,d1>& s1,
        CELL_MANAGER_COLOR<TV>& cm1,Args&& ...args);
    void Initialize_Helper(int subcell,ARRAY<int>& diffs,const TV_INT& index_offset);
    template<int d1,class ...Args> 
    void Store_Cell_Manager(const BASIS_STENCIL_UNIFORM<TV,d1>& s1,CELL_MANAGER_COLOR<TV>& cm1,Args&& ...args);
    void Store_Cell_Manager();
    void Mark_Active_Cells(T tol=0);
    void Build_Matrix(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& matrix);
    template<class ...Args>
    void Build_Matrix_With_Contract(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& matrix,Args&& ...args);
    template<class ...Args>
    void Build_Matrix_With_Contract(int c,SPARSE_MATRIX_FLAT_MXN<T>& matrix,int index,
        const ARRAY<ARRAY<T> >& vec,Args&& ...args);
    template<class ...Args>
    void Setup_Contract_Instructions(ARRAY<ARRAY<T> >& contract_row,ARRAY<ARRAY<VECTOR<int,2> > >& instructions,
        ARRAY<ARRAY<const ARRAY<T>*> >& vecs,const ARRAY<ARRAY<int> >& flat_diff,ARRAY<ARRAY<int> >& flatter_diff,
        ARRAY<CELL_MANAGER_COLOR<TV>*>& pruned_cm,ARRAY<int>& contract_index,int index,const ARRAY<ARRAY<T> >& vec,
        Args&& ...args);
    void Setup_Contract_Instructions(ARRAY<ARRAY<T> >& contract_row,ARRAY<ARRAY<VECTOR<int,2> > >& instructions,
        ARRAY<ARRAY<const ARRAY<T>*> >& vecs,const ARRAY<ARRAY<int> >& flat_diff,ARRAY<ARRAY<int> >& flatter_diff,
        ARRAY<CELL_MANAGER_COLOR<TV>*>& pruned_cm,ARRAY<int>& contract_index);
};
}
#endif

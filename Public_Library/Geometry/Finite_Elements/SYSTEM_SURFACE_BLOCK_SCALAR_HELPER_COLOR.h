//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR
//#####################################################################

#ifndef __SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR__
#define __SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR__

#include <Core/Arrays/ARRAY.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include <Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>

namespace PhysBAM{

template<class TV> class CELL_MANAGER_COLOR;
template<class TV,int d> struct BASIS_STENCIL_UNIFORM;
template<class T> class SPARSE_MATRIX_FLAT_MXN;

template<class TV>
class SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

public:
    
    ARRAY<MATRIX_MXN<T> > data;
    ARRAY<ARRAY<T> > rhs_data;
    CELL_DOMAIN_INTERFACE_COLOR<TV> *cdi;
    CELL_MANAGER_COLOR<TV> *cm;
    ARRAY<int> flat_diff;

    SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR() = default;
    SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR(const SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR&) = delete;
    void operator=(const SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR&) = delete;
    
    template<int d> 
    void Initialize(CELL_DOMAIN_INTERFACE_COLOR<TV> &cdi_input,const BASIS_STENCIL_UNIFORM<TV,d>& s,CELL_MANAGER_COLOR<TV>& cm_input);
    void Mark_Active_Cells(T tol=0);
    void Build_Matrix(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& matrix,ARRAY<T>& constraint_rhs);
    void Resize();
};
}
#endif

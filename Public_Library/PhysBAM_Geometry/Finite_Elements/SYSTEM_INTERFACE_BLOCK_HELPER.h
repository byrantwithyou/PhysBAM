//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYSTEM_INTERFACE_BLOCK_HELPER
//#####################################################################

#ifndef __SYSTEM_INTERFACE_BLOCK_HELPER__
#define __SYSTEM_INTERFACE_BLOCK_HELPER__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>

namespace PhysBAM{

template<class TV> class CELL_DOMAIN_INTERFACE;
template<class TV> class CELL_MANAGER;
template<class TV,int d> struct BASIS_STENCIL_UNIFORM;
template<class T> class SPARSE_MATRIX_FLAT_MXN;

template<class TV>
class SYSTEM_INTERFACE_BLOCK_HELPER:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

public:
    
    VECTOR<MATRIX_MXN<T>,2> data;
    CELL_DOMAIN_INTERFACE<TV> *cdi;        
    CELL_MANAGER<TV> *cm;
    ARRAY<int> flat_diff;
    
    template<int d> 
    void Initialize(const BASIS_STENCIL_UNIFORM<TV,d>& s,CELL_MANAGER<TV>& cm_input,CELL_DOMAIN_INTERFACE<TV> &cdi_input);
    void Mark_Active_Cells(T tol=0);
    void Build_Matrix(VECTOR<SPARSE_MATRIX_FLAT_MXN<T>,2>& matrix);
};
}
#endif

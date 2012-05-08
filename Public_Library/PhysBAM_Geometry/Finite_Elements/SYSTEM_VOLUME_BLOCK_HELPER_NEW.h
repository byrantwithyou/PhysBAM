//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYSTEM_VOLUME_BLOCK_HELPER_NEW
//#####################################################################

#ifndef __SYSTEM_VOLUME_BLOCK_HELPER_NEW__
#define __SYSTEM_VOLUME_BLOCK_HELPER_NEW__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>

namespace PhysBAM{

template<class TV> class CELL_DOMAIN_INTERFACE_NEW;
template<class TV> class CELL_MANAGER_NEW;
template<class TV,int d> class BASIS_STENCIL_UNIFORM;
template<class T> class SPARSE_MATRIX_FLAT_MXN;

template<class TV>
class SYSTEM_VOLUME_BLOCK_HELPER_NEW:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

public:
    VECTOR<MATRIX_MXN<T>,2> data;
    CELL_DOMAIN_INTERFACE_NEW<TV>* cdi;
    CELL_MANAGER_NEW<TV> *cm0,*cm1;
    ARRAY<int> flat_diff;

    template<int d0,int d1> 
    void Initialize(const BASIS_STENCIL_UNIFORM<TV,d0>& s0,const BASIS_STENCIL_UNIFORM<TV,d1>& s1,
        CELL_MANAGER_NEW<TV>& cm0_input,CELL_MANAGER_NEW<TV>& cm1_input,CELL_DOMAIN_INTERFACE_NEW<TV> &cdi_input);
    void Mark_Active_Cells(T tol=0);
    void Build_Matrix(VECTOR<SPARSE_MATRIX_FLAT_MXN<T>,2>& matrix);
};
}
#endif

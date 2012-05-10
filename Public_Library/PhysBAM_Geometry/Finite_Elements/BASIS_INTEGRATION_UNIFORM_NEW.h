//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BASIS_INTEGRATION_UNIFORM_NEW
//#####################################################################
#ifndef __BASIS_INTEGRATION_UNIFORM_NEW__
#define __BASIS_INTEGRATION_UNIFORM_NEW__

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/STATIC_TENSOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_INTERFACE_BLOCK_NEW.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK_NEW.h>

namespace PhysBAM{

template<class TV> class GRID;
template<class TV> class CELL_MAPPING;
template<class T,int rank,int d> class STATIC_TENSOR;

template<class TV,int static_degree>
class BASIS_INTEGRATION_UNIFORM_NEW:public NONCOPYABLE
{
public:
    typedef SYSTEM_VOLUME_BLOCK_NEW<TV,static_degree> VOLUME_BLOCK;
    typedef SYSTEM_INTERFACE_BLOCK_NEW<TV,static_degree> INTERFACE_BLOCK;
    typedef typename VOLUME_BLOCK::OPEN_ENTRY OPEN_ENTRY;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE T_FACE;
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    const GRID<TV>& grid;
    const GRID<TV>& phi_grid;
    CELL_DOMAIN_INTERFACE_NEW<TV>& cdi;
    enum WORKAROUND {subcell_elements=(TV::m==3?5:2)};

    const ARRAY<T,TV_INT>& phi;
    STATIC_TENSOR<bool,TV::m,static_degree+1> volume_monomials_needed,surface_monomials_needed;

    ARRAY<VOLUME_BLOCK*> volume_blocks;
    ARRAY<INTERFACE_BLOCK*> interface_blocks;

    BASIS_INTEGRATION_UNIFORM_NEW(const GRID<TV>& grid_input,const GRID<TV>& phi_grid_input,
        const ARRAY<T,TV_INT>& phi_input,CELL_DOMAIN_INTERFACE_NEW<TV>& cdi_input);
    ~BASIS_INTEGRATION_UNIFORM_NEW();

    void Compute_Entries();
    void Compute_Open_Entries();
    void Add_Uncut_Cell(const TV_INT& cell,int enclose_inside);
    void Add_Uncut_Fine_Cell(const TV_INT& cell,int block,int enclose_inside);
    void Add_Cut_Fine_Cell(const TV_INT& cell,int block,ARRAY<T_FACE>& interface,ARRAY<T_FACE>& sides,
        int direction,bool enclose_inside,int cut_cell_index,MATRIX<T,TV::m>& base_orientation);
    template<int d0,int d1>
    void Add_Volume_Block(SYSTEM_VOLUME_BLOCK_HELPER_NEW<TV>& helper,const BASIS_STENCIL_UNIFORM<TV,d0>& s0,
        const BASIS_STENCIL_UNIFORM<TV,d1>& s1,const VECTOR<T,2>& scale);
    template<int d>
    void Add_Interface_Block(SYSTEM_INTERFACE_BLOCK_HELPER_NEW<TV>& helper,const BASIS_STENCIL_UNIFORM<TV,d>& s,
        T scale,bool ignore_orientation);
};
}
#endif

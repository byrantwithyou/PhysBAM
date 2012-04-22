//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BASIS_INTEGRATION_UNIFORM
//#####################################################################
#ifndef __BASIS_INTEGRATION_UNIFORM__
#define __BASIS_INTEGRATION_UNIFORM__

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Vectors/STATIC_TENSOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_INTERFACE_BLOCK.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>

namespace PhysBAM{

template<class TV> class GRID;
template<class TV> class CELL_MAPPING;
template<class T,int rank,int d> class STATIC_TENSOR;

template<class TV,int static_degree>
class BASIS_INTEGRATION_UNIFORM:public NONCOPYABLE
{
public:
    typedef SYSTEM_VOLUME_BLOCK<TV,static_degree> VOLUME_BLOCK;
    typedef SYSTEM_INTERFACE_BLOCK<TV,static_degree> INTERFACE_BLOCK;
    typedef typename VOLUME_BLOCK::OPEN_ENTRY OPEN_ENTRY;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE T_FACE;
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    const GRID<TV>& grid;
    const GRID<TV>& phi_grid;
    CELL_DOMAIN_INTERFACE<TV>& cdi;
    enum WORKAROUND {subcell_elements=(TV::m==3?5:2)};

    const ARRAY<T,TV_INT>& phi;
    STATIC_TENSOR<bool,TV::m,static_degree+1> volume_monomials_needed,surface_monomials_needed;
    int coarse_factor;
    RANGE<TV_INT> double_coarse_range,coarse_range;

    ARRAY<VOLUME_BLOCK*> volume_blocks;
    ARRAY<INTERFACE_BLOCK*> interface_blocks;

    BASIS_INTEGRATION_UNIFORM(const GRID<TV>& grid_input,const GRID<TV>& phi_grid_input,const ARRAY<T,TV_INT>& phi_input,CELL_DOMAIN_INTERFACE<TV>& cdi_input);
    ~BASIS_INTEGRATION_UNIFORM();

    void Compute();
    void Compute_Open_Entries();
    void Cut_Elements(ARRAY<ARRAY<PAIR<T_FACE,int> >,TV_INT>& cut_elements,const ARRAY<T_FACE>& elements,
        const RANGE<TV_INT>& range,const RANGE<TV>& domain,int dir,int e);
    void Add_Uncut_Coarse_Cell(const TV_INT& coarse_cell,int inside);
    void Add_Uncut_Fine_Cell(const TV_INT& cell,int block,int inside);
    template<int d0,int d1>
    int Add_Block(const BASIS_STENCIL_UNIFORM<TV,d0>& s0,const BASIS_STENCIL_UNIFORM<TV,d1>& s1,
        CELL_MANAGER<TV>& cm0,CELL_MANAGER<TV>& cm1,const VECTOR<T,2>& scale);
    template<int d>
    int Add_Block(const BASIS_STENCIL_UNIFORM<TV,d>& s,CELL_MANAGER<TV>& cm,
        T scale,bool ignore_orientation);
    void Add_Cut_Subcell(const ARRAY<PAIR<T_FACE,int> >& side_elements,const ARRAY<PAIR<T_FACE,int> >& interface_elements,
        const TV_INT& cell,const TV_INT& subcell_cell,int dir,bool enclose_inside,int block,int element_base);
};
}
#endif

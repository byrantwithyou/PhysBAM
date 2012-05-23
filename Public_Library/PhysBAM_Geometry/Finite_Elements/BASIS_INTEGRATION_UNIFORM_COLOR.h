//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BASIS_INTEGRATION_UNIFORM_COLOR
//#####################################################################
#ifndef __BASIS_INTEGRATION_UNIFORM_COLOR__
#define __BASIS_INTEGRATION_UNIFORM_COLOR__

#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/STATIC_TENSOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_SURFACE_BLOCK_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK_COLOR.h>

namespace PhysBAM{

template<class TV> class GRID;
template<class TV> class CELL_MAPPING;
template<class T,int rank,int d> class STATIC_TENSOR;

template<class TV,int static_degree>
class BASIS_INTEGRATION_UNIFORM_COLOR:public NONCOPYABLE
{
public:
    typedef SYSTEM_VOLUME_BLOCK_COLOR<TV,static_degree> VOLUME_BLOCK;
    typedef SYSTEM_SURFACE_BLOCK_COLOR<TV,static_degree> SURFACE_BLOCK;
    typedef typename VOLUME_BLOCK::OPEN_ENTRY OPEN_ENTRY;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE T_FACE;
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    const GRID<TV>& grid;
    const GRID<TV>& phi_grid;
    const ARRAY<T,TV_INT>& phi_value;
    const ARRAY<int,TV_INT>& phi_color;

    CELL_DOMAIN_INTERFACE_COLOR<TV>& cdi;

    STATIC_TENSOR<bool,TV::m,static_degree+1> volume_monomials_needed;
    STATIC_TENSOR<bool,TV::m,static_degree+1> surface_monomials_needed;

    ARRAY<VOLUME_BLOCK*> volume_blocks;
    ARRAY<SURFACE_BLOCK*> surface_blocks;

    BASIS_INTEGRATION_UNIFORM_COLOR(const GRID<TV>& grid_input,const GRID<TV>& phi_grid_input,
        const ARRAY<T,TV_INT>& phi_value_input,const ARRAY<int,TV_INT>& phi_color_input,CELL_DOMAIN_INTERFACE_COLOR<TV>& cdi_input);
    ~BASIS_INTEGRATION_UNIFORM_COLOR();

    void Compute_Entries(VECTOR<ARRAY<VECTOR_ND<T> >,TV::m>* f_surface);
    void Compute_Open_Entries();
    void Add_Uncut_Cell(const TV_INT& cell,int color);
    void Add_Uncut_Fine_Cell(const TV_INT& cell,int block,int color);
    void Add_Cut_Fine_Cell(const TV_INT& cell,int block,const TV& block_offset,ARRAY<TRIPLE<T_FACE,int,int> >& surface,ARRAY<PAIR<T_FACE,int> >& sides,
        const ARRAY<MATRIX<T,TV::m> >& base_orientation,VECTOR<ARRAY<VECTOR_ND<T> >,TV::m>* f_surface,const ARRAY<int>& map_constraints,
        const HASHTABLE<VECTOR<int,2>,int>& ht_color_pairs);

    template<int d0,int d1>
    void Add_Volume_Block(SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>& helper,const BASIS_STENCIL_UNIFORM<TV,d0>& s0,
        const BASIS_STENCIL_UNIFORM<TV,d1>& s1,const ARRAY<T>& scale);
    template<int d>
    void Add_Surface_Block(SYSTEM_SURFACE_BLOCK_HELPER_COLOR<TV>& helper,const BASIS_STENCIL_UNIFORM<TV,d>& s,
        ANALYTIC_BOUNDARY_CONDITIONS_COLOR<TV>* abc,int axis,T scale);
};
}
#endif

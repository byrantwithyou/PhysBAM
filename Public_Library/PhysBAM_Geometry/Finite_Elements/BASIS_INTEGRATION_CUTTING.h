//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BASIS_INTEGRATION_CUTTING
//#####################################################################
#ifndef __BASIS_INTEGRATION_CUTTING__
#define __BASIS_INTEGRATION_CUTTING__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/STATIC_TENSOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>
namespace PhysBAM{

template<class TV> class GRID;
template<class TV> class CELL_MAPPING;
template<class T,int rank,int d> class STATIC_TENSOR;

template<class TV,int static_degree>
class BASIS_INTEGRATION_CUTTING:public NONCOPYABLE
{
public:
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE T_FACE;
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    const GRID<TV>& grid;
    const GRID<TV>& phi_grid;
    enum WORKAROUND {subcell_elements=(TV::m==3?5:2)};

    enum BOUNDARY_CONDITION {undefined, periodic, dirichlet, neumann};

    RANGE<TV_INT> boundary_conditions;
    const ARRAY<T,TV_INT>& phi;
    STATIC_TENSOR<bool,TV::m,static_degree+2> monomials_needed;
    int coarse_factor;
    RANGE<TV_INT> double_coarse_range,coarse_range;

    struct OVERLAP_VOLUME_POLYNOMIALS
    {
        TV_INT index_offset0, index_offset1;
        int subcell; // flags
        STATIC_POLYNOMIAL<T,TV::m,static_degree> polynomial;
    };

    struct VOLUME_MATRIX_ENTRY
    {
        TV_INT index0, index1;
        T x;

        bool operator< (const VOLUME_MATRIX_ENTRY& me) const
        {
            if(index0!=me.index0) return LEXICOGRAPHIC_COMPARE()(index0,me.index0);
            return LEXICOGRAPHIC_COMPARE()(index1,me.index1);
        }

        void Merge(const VOLUME_MATRIX_ENTRY& me)
        {
            x+=me.x;
        }
    };

    struct VOLUME_BLOCK
    {
        CELL_MAPPING<TV>* cm0,*cm1;
        SYSTEM_MATRIX_HELPER<T>* helper;
        ARRAY<OVERLAP_VOLUME_POLYNOMIALS> overlap;
        ARRAY<VOLUME_MATRIX_ENTRY> open_entries,open_subcell_entries[1<<TV::m];
        T scale;
    };

    ARRAY<VOLUME_BLOCK*> volume_blocks;

    struct OVERLAP_INTERFACE_POLYNOMIALS
    {
        TV_INT index_offset;
        int subcell; // flags indicating fine cells
        STATIC_POLYNOMIAL<T,TV::m,static_degree> polynomial;
    };

    struct INTERFACE_BLOCK
    {
        CELL_MAPPING<TV>* cm;
        SYSTEM_MATRIX_HELPER<T>* helper;
        ARRAY<OVERLAP_INTERFACE_POLYNOMIALS> overlap;
        T scale;
    };

    ARRAY<INTERFACE_BLOCK*> interface_blocks;

    BASIS_INTEGRATION_CUTTING(const RANGE<TV_INT> boundary_conditions_input,const GRID<TV>& grid_input,const GRID<TV>& phi_grid_input,
        const ARRAY<T,TV_INT>& phi_input);
    ~BASIS_INTEGRATION_CUTTING();

    void Compute();
    void Compute_Open_Entries();
    void Cut_Elements(ARRAY<ARRAY<PAIR<T_FACE,int> >,TV_INT>& cut_elements,const ARRAY<T_FACE>& elements,
        const RANGE<TV_INT>& range,const RANGE<TV>& domain,int dir,int e);
    void Add_Uncut_Coarse_Cell(const TV_INT& coarse_cell,int inside);
    void Add_Uncut_Fine_Cell(const TV_INT& cell,int block,int inside);
    void Apply_Matrix_Entry(VOLUME_BLOCK* vb,const TV_INT& cell,bool inside,VOLUME_MATRIX_ENTRY me);
    template<int d0,int d1>
    int Add_Block(SYSTEM_MATRIX_HELPER<T>& helper,const BASIS_STENCIL_UNIFORM<TV,d0>& s0,const BASIS_STENCIL_UNIFORM<TV,d1>& s1,
        CELL_MAPPING<TV>& cm0,CELL_MAPPING<TV>& cm1,T scale);
    template<int d>
    int Add_Block(SYSTEM_MATRIX_HELPER<T>& helper,const BASIS_STENCIL_UNIFORM<TV,d>& s,CELL_MAPPING<TV>& cm,T scale);
    void Add_Cut_Subcell(const ARRAY<PAIR<T_FACE,int> >& side_elements,const ARRAY<PAIR<T_FACE,int> >& interface_elements,
        const TV_INT& cell,const TV_INT& subcell_cell,int dir,bool enclose_inside,int block,int element_base);
};
}
#endif

//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BASIS_INTEGRATION_UNIFORM
//#####################################################################
#ifndef __BASIS_INTEGRATION_UNIFORM__
#define __BASIS_INTEGRATION_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Symbolics/MULTIVARIATE_POLYNOMIAL.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
namespace PhysBAM{

template<class TV> class GRID;
template<class TV> class BASIS_STENCIL_UNIFORM;
template<class TV> class CELL_MAPPING;
template<class T> class SYSTEM_MATRIX_HELPER;

template<class TV>
class BASIS_INTEGRATION_UNIFORM:public NONCOPYABLE
{
public:
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE T_FACE;
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    const GRID<TV>& grid;
    const GRID<TV>& phi_grid;

    enum BOUNDARY_CONDITION {undefined, periodic, dirichlet, neumann};

    RANGE<TV_INT> boundary_conditions;
    const BASIS_STENCIL_UNIFORM<TV>& s0,s1;
    CELL_MAPPING<TV>& cm0,&cm1;
    const ARRAY<T,TV_INT>& phi;

    struct MATRIX_ENTRY
    {
        TV_INT index0, index1;
        T x;

        bool operator< (const MATRIX_ENTRY& me) const
        {
            if(index0!=me.index0) return LEXICOGRAPHIC_COMPARE()(index0,me.index0);
            return LEXICOGRAPHIC_COMPARE()(index1,me.index1);
        }

        void Merge(const MATRIX_ENTRY& me)
        {
            x+=me.x;
        }
    };
    ARRAY<MATRIX_ENTRY> open_entries; // stencil with no boundary conditions

    struct OVERLAP_POLYNOMIALS
    {
        TV_INT index_offset0, index_offset1;
        RANGE<TV_INT> range; // Subset of [-1,1)
        MULTIVARIATE_POLYNOMIAL<TV> polynomial;
    };
    ARRAY<OVERLAP_POLYNOMIALS> overlap_polynomials;

    BASIS_INTEGRATION_UNIFORM(const RANGE<TV_INT> boundary_conditions_input,const GRID<TV>& grid_input,const GRID<TV>& phi_grid_input,
        const BASIS_STENCIL_UNIFORM<TV>& s0_input,const BASIS_STENCIL_UNIFORM<TV>& s1_input,CELL_MAPPING<TV>& cm0_input,
        CELL_MAPPING<TV>& cm1_input,const ARRAY<T,TV_INT>& phi_input);
    ~BASIS_INTEGRATION_UNIFORM();

    void Compute_Matrix(SYSTEM_MATRIX_HELPER<T>& helper);
    void Add_Uncut_Stencil(SYSTEM_MATRIX_HELPER<T>& helper,const TV_INT& cell,bool inside);
    void Add_Cut_Stencil(SYSTEM_MATRIX_HELPER<T>& helper,const ARRAY<T_FACE>& elements,const TV_INT& cell,int dir,bool inside,const TV_INT& sub_cell,const OVERLAP_POLYNOMIALS& op);
    void Cut_Elements(ARRAY<ARRAY<T_FACE>,TV_INT>& cut_elements,const ARRAY<T_FACE>& elements,const RANGE<TV_INT>& range,const RANGE<TV>& domain,int dir);
    void Compute_Overlap_Polynomials();
    void Compute_Open_Entries();
};
}
#endif

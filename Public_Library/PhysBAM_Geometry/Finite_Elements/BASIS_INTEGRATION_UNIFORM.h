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
namespace PhysBAM{

template<class TV> class GRID;
template<class TV> class BASIS_STENCIL_UNIFORM;
template<class T> class SYSTEM_MATRIX_HELPER;

template<class TV>
class BASIS_INTEGRATION_UNIFORM:public NONCOPYABLE
{
public:
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    const GRID<TV>& grid;

    enum BOUNDARY_CONDITION {undefined, periodic, dirichlet, neumann};

    RANGE<TV_INT> boundary_conditions;

    struct MATRIX_ENTRY
    {
        TV_INT index0, index1;
        T x;
    };

    struct OVERLAP_POLYNOMIALS
    {
        TV_INT index_offset0, index_offset1;
        RANGE<TV_INT> range; // Subset of [-1,1)
        MULTIVARIATE_POLYNOMIAL<TV> polynomial;
    };

    BASIS_INTEGRATION_UNIFORM(const GRID<TV>& grid_input);
    ~BASIS_INTEGRATION_UNIFORM();

    void Compute_Matrix(SYSTEM_MATRIX_HELPER<T>& helper,const BASIS_STENCIL_UNIFORM<TV>& s0, const BASIS_STENCIL_UNIFORM<TV>& s1, const ARRAY<int,TV_INT>& index_map0, const ARRAY<int,TV_INT>& index_map1);
};
}
#endif

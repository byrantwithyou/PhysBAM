//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BASIS_INTEGRATION_BOUNDARY_UNIFORM
//#####################################################################
#ifndef __BASIS_INTEGRATION_BOUNDARY_UNIFORM__
#define __BASIS_INTEGRATION_BOUNDARY_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Symbolics/MULTIVARIATE_POLYNOMIAL.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{

template<class TV> class GRID;
template<class TV> class BASIS_STENCIL_BOUNDARY_UNIFORM;
template<class TV> class BASIS_STENCIL_UNIFORM;
template<class TV> class CELL_MAPPING;
template<class T> class SYSTEM_MATRIX_HELPER;

template<class TV>
class BASIS_INTEGRATION_BOUNDARY_UNIFORM:public NONCOPYABLE
{
public:
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE T_FACE;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_OBJECT;
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    const GRID<TV>& grid;

    const BASIS_STENCIL_UNIFORM<TV>& stencil;
    const BASIS_STENCIL_BOUNDARY_UNIFORM<TV>& boundary;
    CELL_MAPPING<TV>& cm;

    BASIS_INTEGRATION_BOUNDARY_UNIFORM(const GRID<TV>& grid_input,const BASIS_STENCIL_UNIFORM<TV>& stencil_input,const BASIS_STENCIL_BOUNDARY_UNIFORM<TV>& boundary_input,CELL_MAPPING<TV>& cm_input);
    ~BASIS_INTEGRATION_BOUNDARY_UNIFORM();

    void Compute_Matrix(SYSTEM_MATRIX_HELPER<T>& helper);
};
}
#endif

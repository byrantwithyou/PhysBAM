//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __CUT_CELL_PROJECTION__
#define __CUT_CELL_PROJECTION__
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/GRID.h>
namespace PhysBAM
{
// bc_type: 0=periodic, 1=free, 2=slip
template<class TV>
struct CUT_CELL_PROJECTION
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    const ARRAY<T,TV_INT>* object_phi=0;
    const ARRAY<T,TV_INT>* surface_phi=0;
    ARRAY<T,TV_INT>* pressure=0;
    ARRAY<bool,FACE_INDEX<TV::m> >* valid_u=0;
    ARRAY<bool,TV_INT>* valid_p=0;
    VECTOR<int,2*TV::m> bc_type;
    bool use_preconditioner=true;
    bool print_matrix=false;
    bool print_residual=false;
    bool test_system=false;
    T solver_tolerance=std::numeric_limits<T>::epsilon();
    int solver_iterations=INT_MAX;
    std::function<TV(const TV& X)> bc_v=0;
    std::function<T(const TV& X)> bc_p=0;
    bool use_warm_start=false;
    
    void Cut_Cell_Projection(const GRID<TV>& grid,int ghost,
        ARRAY<T,FACE_INDEX<TV::m> >& u,T density,T dt);
};
}
#endif

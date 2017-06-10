//#####################################################################
// Copyright 2005, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FAST_MARCHING_METHOD_UNIFORM  
//#####################################################################
#ifndef __FAST_MARCHING_METHOD_UNIFORM__
#define __FAST_MARCHING_METHOD_UNIFORM__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <functional>
namespace PhysBAM{
template<class TV>
class FAST_MARCHING_METHOD_UNIFORM
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    std::function<bool(const FACE_INDEX<TV::m>& face)> Neighbor_Visible=0;
    bool add_seed_indices_for_ghost_cells=false;
    const ARRAY<TV_INT>* seed_indices=0;
    int process_sign=0;
    bool correct_interface_phi=true;
    struct ESTIMATION_PARAMETERS
    {
        T iterative_tolerance=(T)1e-2;
        bool refine_with_iterative_solver=true;
        int iterations=10;
        T iterative_drift_fraction=(T).1;
    } estimation_parameters;

    FAST_MARCHING_METHOD_UNIFORM()=default;
    ~FAST_MARCHING_METHOD_UNIFORM()=default;

//#####################################################################
    void Fast_Marching_Method(const GRID<TV>& grid,int ghost,
        ARRAY<T,TV_INT>& phi,const T stopping_distance=0);
//#####################################################################
};
}
#endif

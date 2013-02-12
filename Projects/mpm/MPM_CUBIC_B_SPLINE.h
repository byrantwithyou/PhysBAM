//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_CUBIC_B_SPLINE
//#####################################################################
#ifndef __MPM_CUBIC_B_SPLINE__
#define __MPM_CUBIC_B_SPLINE__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include "MPM_GRID_BASIS_FUNCTION.h"
namespace PhysBAM{

template<class TV,int order>
class MPM_CUBIC_B_SPLINE:public MPM_GRID_BASIS_FUNCTION<TV,order>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef MPM_GRID_BASIS_FUNCTION<TV,order> BASE;
public:
    using BASE::IN;

    MPM_CUBIC_B_SPLINE();
    ~MPM_CUBIC_B_SPLINE();
    
    void Build_Weights_And_Grad_Weight_Over_Weights(const TV& X,const GRID<TV>& grid,TV_INT& influence_corner,VECTOR<TV,IN>& weight,VECTOR<TV,IN>& grad_weight_over_weight) PHYSBAM_OVERRIDE;
    void Build_Weights_And_Grad_Weights_Exact(const TV& X,const GRID<TV>& grid,TV_INT& influence_corner,ARRAY<T>& weight,ARRAY<TV>& grad_weight) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

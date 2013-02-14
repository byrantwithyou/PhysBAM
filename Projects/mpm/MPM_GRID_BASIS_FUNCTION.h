//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_GRID_BASIS_FUNCTION
//#####################################################################
#ifndef __MPM_GRID_BASIS_FUNCTION__
#define __MPM_GRID_BASIS_FUNCTION__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
namespace PhysBAM{

template<class TV,int order>
class MPM_GRID_BASIS_FUNCTION
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    enum WORKAROUND {IN=order+1};
    MPM_GRID_BASIS_FUNCTION(){}
    ~MPM_GRID_BASIS_FUNCTION(){}
    
    virtual void Build_Weights_And_Grad_Weight_Over_Weights(const TV&X,const GRID<TV>& grid,TV_INT& influence_corner,VECTOR<TV,IN>& weight,VECTOR<TV,IN>& grad_weight_over_weight)=0;
    virtual void Build_Weights_And_Grad_Weights_Exact(const TV& X,const GRID<TV>& grid,TV_INT& influence_corner,ARRAY<T,TV_INT>& weight,ARRAY<TV,TV_INT>& grad_weight)=0;
//#####################################################################
};
}
#endif

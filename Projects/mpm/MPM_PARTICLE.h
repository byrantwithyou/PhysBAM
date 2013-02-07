//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_PARTICLE
//#####################################################################
#ifndef __MPM_PARTICLE__
#define __MPM_PARTICLE__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
namespace PhysBAM{

template<class TV,int IN> // IN means Influence N, which is 4 nodes for cubic B-spline, 2 nodes for bi/tri-linear
class MPM_PARTICLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

public:
    
    T mu;
    T lambda;
    TV X;
    TV V;
    T mass;

    T volume0;

    TV_INT influence_corner;
    VECTOR<TV,IN> weights;
    VECTOR<TV,IN> grad_weights;

    MPM_PARTICLE();

    MPM_PARTICLE(const T muI,const T lambdaI,const TV XI,const TV VI,const T massI,const GRID<TV>& grid);

    void Evaluate_Weights_And_Grad_Weights(const GRID<TV>& grid);

    ~MPM_PARTICLE();

//#####################################################################
};
}
#endif

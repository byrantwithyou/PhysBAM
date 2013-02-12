//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/pow.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include "MPM_PARTICLES.h"
#include "MPM_CONSTITUTIVE_MODEL.h"
#include "MPM_CUBIC_B_SPLINE.h"
#include "MPM_SIMULATION.h"
using namespace PhysBAM;

int main(int argc,char *argv[])
{
    static const int dimension=2;
    typedef float T;
    typedef VECTOR<T,dimension> TV;
    typedef VECTOR<int,dimension> TV_INT;

    TV_INT grid_counts(6,7);
    RANGE<TV> grid_box(TV(0,0),TV(5,6));
    GRID<TV> grid(grid_counts,grid_box);

    RANGE<TV_INT> range(TV_INT(),TV_INT()+4);
    for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
        LOG::cout<<it.index<<std::endl;}


    // {
    //     TV X(2.75,3.5);
    //     TV_INT influence_corner;
    //     VECTOR<TV,4> weight;
    //     VECTOR<TV,4> grad_weight_over_weight;
    //     MPM_CUBIC_B_SPLINE<TV,3> Bspline;
    //     Bspline.Build_Weights_And_Grad_Weight_Over_Weights(X,grid,influence_corner,weight,grad_weight_over_weight);
    //     LOG::cout<<influence_corner<<std::endl;
    //     LOG::cout<<weight<<std::endl;
    //     LOG::cout<<grad_weight_over_weight<<std::endl;
    //     ARRAY<T> weights; weights.Resize(16);
    //     ARRAY<TV> grad_weights;grad_weights.Resize(16);
    //     Bspline.Build_Weights_And_Grad_Weights_Exact(X,grid,influence_corner,weights,grad_weights);
    //     LOG::cout<<weights<<std::endl;
    //     LOG::cout<<grad_weights<<std::endl;

    // }
    // {
    //     TV X(3.5,3);
    //     TV_INT influence_corner;
    //     VECTOR<TV,4> weight;
    //     VECTOR<TV,4> grad_weight_over_weight;
    //     MPM_CUBIC_B_SPLINE<TV,3> Bspline;
    //     Bspline.Build_Weights_And_Grad_Weight_Over_Weights(X,grid,influence_corner,weight,grad_weight_over_weight);
    //     LOG::cout<<influence_corner<<std::endl;
    //     LOG::cout<<weight<<std::endl;
    //     LOG::cout<<grad_weight_over_weight<<std::endl;
    //     ARRAY<T> weights; weights.Resize(16);
    //     ARRAY<TV> grad_weights;grad_weights.Resize(16);
    //     Bspline.Build_Weights_And_Grad_Weights_Exact(X,grid,influence_corner,weights,grad_weights);
    //     for(int i=0;i<16;i++) LOG::cout<<i<<"   "<<weights(i)<<std::endl;
    //     for(int i=0;i<16;i++) LOG::cout<<i<<"   "<<grad_weights(i)<<std::endl;
    // }

    // MPM_CONSTITUTIVE_MODEL<TV> test;
    // test.Derivative_Test();

    return 0;
}

//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include "MPM_PARTICLE.h"
namespace PhysBAM{

//#####################################################################
// Constructor
//#####################################################################
template<class TV,int IN> MPM_PARTICLE<TV,IN>::
MPM_PARTICLE()
{}

//#####################################################################
// Constructor
//#####################################################################
template<class TV,int IN> MPM_PARTICLE<TV,IN>::
MPM_PARTICLE(const T muI,const T lambdaI,const TV XI,const TV VI,const T massI,const GRID<TV>& grid)
{
    mu=muI;lambda=lambdaI;X=XI;V=VI;mass=massI;
    Evaluate_Weights_And_Grad_Weights(grid);
    //TODO: compute volume0
}

//#####################################################################
// Destructor
//#####################################################################
template<class TV,int IN> MPM_PARTICLE<TV,IN>::
~MPM_PARTICLE()
{}

//#####################################################################
// Function Evaluate_Weights_And_Grad_Weights
//#####################################################################
template<class TV,int IN> MPM_PARTICLE<TV,IN>::
Evaluate_Weights_And_Grad_Weights(const GRID<TV>& grid)
{
    TV_INT my_corner=grid.Cell(X,0);
    if(IN==2) influence_corner=my_corner;
    else if(IN==4) influence_corner=my_corner-1;

    //TODO: build weights

    //TODO: build grad_weights
}

//#####################################################################
template class MPM_PARTICLE<VECTOR<float,1>,2>;
template class MPM_PARTICLE<VECTOR<float,2>,2>;
template class MPM_PARTICLE<VECTOR<float,3>,2>;
template class MPM_PARTICLE<VECTOR<double,1>,2>;
template class MPM_PARTICLE<VECTOR<double,2>,2>;
template class MPM_PARTICLE<VECTOR<double,3>,2>;
template class MPM_PARTICLE<VECTOR<float,1>,4>;
template class MPM_PARTICLE<VECTOR<float,2>,4>;
template class MPM_PARTICLE<VECTOR<float,3>,4>;
template class MPM_PARTICLE<VECTOR<double,1>,4>;
template class MPM_PARTICLE<VECTOR<double,2>,4>;
template class MPM_PARTICLE<VECTOR<double,3>,4>;
}

//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Math_Tools/pow.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
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

    //TODO compute particle initial volume
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
template<class TV,int IN> void MPM_PARTICLE<TV,IN>::
Evaluate_Weights_And_Grad_Weights(const GRID<TV>& grid)
{
    TV_INT my_corner=grid.Cell(X,0);
    if(IN==4){
        influence_corner=my_corner-1;
        for(int n=0;n<IN;n++){
            for(int d=0;d<TV::m;d++){
                T x=grid.One_Over_DX()(d)*(X(d)-(influence_corner(d)+n)*grid.DX()(d));
                T abs_x=abs(x);
                T x_square=abs_x*abs_x;
                T abs_x_cube=x_square*abs_x;
                if(abs_x<1){
                    weights(n)(d)=(T)0.5*abs_x_cube-x_square+(T)0.66667;
                    if(x>0) grad_weights(n)(d)=(T)1.5*x_square-(T)2.0*x;
                    else grad_weights(n)(d)=-(T)1.5*x_square-(T)2.0*x;
                    grad_weights(n)(d)*=grid.One_Over_DX()(d);}
                else if(abs_x<2){
                    weights(n)(d)=-(T)0.16667*abs_x_cube+x_square-(T)2.0*abs_x+(T)1.3333;
                    if(x>0) grad_weights(n)(d)=-(T)0.5*x_square+(T)2.0*x-(T)2.0;
                    else grad_weights(n)(d)=(T)0.5*x_square+(T)2.0*x+(T)2.0;
                    grad_weights(n)(d)*=grid.One_Over_DX()(d);}
                else weights(n)(d)=0;}}}
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

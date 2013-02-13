//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include "Fanfu_Utilities/FLATTEN_INDEX.h"
#include "MPM_CUBIC_B_SPLINE.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int order> MPM_CUBIC_B_SPLINE<TV,order>::
MPM_CUBIC_B_SPLINE()
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int order> MPM_CUBIC_B_SPLINE<TV,order>::
~MPM_CUBIC_B_SPLINE()
{}
//#####################################################################
// Function Build_Weights_And_Grad_Weight_Over_Weights
//#####################################################################
template<class TV,int order> void MPM_CUBIC_B_SPLINE<TV,order>::
Build_Weights_And_Grad_Weight_Over_Weights(const TV& X,const GRID<TV>& grid,TV_INT& influence_corner,VECTOR<TV,IN>& weight,VECTOR<TV,IN>& grad_weight_over_weight)
{  
    static T eps=1e-5;
    TV_INT my_corner=grid.Cell(X,0);
    influence_corner=my_corner-1;
    for(int n=0;n<IN;n++){
        for(int d=0;d<TV::m;d++){
            T x=grid.one_over_dX(d)*(X(d)-(influence_corner(d)+n)*grid.dX(d)-grid.domain.min_corner(d));
            T abs_x=abs(x),x_square=abs_x*abs_x,abs_x_cube=x_square*abs_x;
            if(abs_x<1){
                weight(n)(d)=(T)0.5*abs_x_cube-x_square+(T)0.66666666666666667;
                if(x>0) grad_weight_over_weight(n)(d)=(T)1.5*x_square-(T)2.0*x;
                else grad_weight_over_weight(n)(d)=-(T)1.5*x_square-(T)2.0*x;
                grad_weight_over_weight(n)(d)*=grid.one_over_dX(d);
                grad_weight_over_weight(n)(d)/=weight(n)(d);}
            else if(abs_x<2){
                weight(n)(d)=-(T)0.16666666666666666667*abs_x_cube+x_square-(T)2.0*abs_x+(T)1.333333333333333333;
                if(abs(weight(n)(d))<eps){
                    weight(n)(d)=0;
                    grad_weight_over_weight(n)(d)=0;
                    continue;}
                if(x>0) grad_weight_over_weight(n)(d)=-(T)0.5*x_square+(T)2.0*x-(T)2.0;
                else grad_weight_over_weight(n)(d)=(T)0.5*x_square+(T)2.0*x+(T)2.0;
                grad_weight_over_weight(n)(d)*=grid.one_over_dX(d);
                grad_weight_over_weight(n)(d)/=weight(n)(d);}
            else{
                weight(n)(d)=0;
                grad_weight_over_weight(n)(d)=0;}}}
}
//#####################################################################
// Function Build_Weights_And_Grad_Weights_Exact
//#####################################################################
template<class TV,int order> void MPM_CUBIC_B_SPLINE<TV,order>::
Build_Weights_And_Grad_Weights_Exact(const TV& X,const GRID<TV>& grid,TV_INT& influence_corner,ARRAY<T>& weight,ARRAY<TV>& grad_weight)
{  
    VECTOR<TV,IN> weight_short;
    VECTOR<TV,IN> grad_weight_over_weight_short;
    Build_Weights_And_Grad_Weight_Over_Weights(X,grid,influence_corner,weight_short,grad_weight_over_weight_short);
    // LOG::cout<<X<<std::endl;
    // LOG::cout<<influence_corner<<std::endl;
    // LOG::cout<<weight_short<<std::endl;

    TV_INT local_grid_count;
    local_grid_count.Fill(IN);
    RANGE<TV_INT> range(TV_INT(),TV_INT()+IN);    
    weight.Fill(T(1));
    for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
        int ind=Flatten_Index(it.index,local_grid_count);
        for(int d=0;d<TV::m;d++) weight(ind)*=weight_short(it.index(d))(d);
        TV gw;
        for(int d=0;d<TV::m;d++){
            TV lg;
            lg(d)=grad_weight_over_weight_short(it.index(d))(d);
            gw+=lg;}
        grad_weight(ind)=gw;
        grad_weight(ind)*=weight(ind);}
}
//#####################################################################
template class MPM_CUBIC_B_SPLINE<VECTOR<float,2>,3>;
template class MPM_CUBIC_B_SPLINE<VECTOR<float,3>,3>;
template class MPM_CUBIC_B_SPLINE<VECTOR<double,2>,3>;
template class MPM_CUBIC_B_SPLINE<VECTOR<double,3>,3>;
}

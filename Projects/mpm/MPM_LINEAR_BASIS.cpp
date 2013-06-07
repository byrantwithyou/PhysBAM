//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include "MPM_LINEAR_BASIS.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int order> MPM_LINEAR_BASIS<TV,order>::
MPM_LINEAR_BASIS()
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int order> MPM_LINEAR_BASIS<TV,order>::
~MPM_LINEAR_BASIS()
{}
//#####################################################################
// Function Build_Weights
//#####################################################################
template<class TV,int order> void MPM_LINEAR_BASIS<TV,order>::
Build_Weights(const TV& X,const GRID<TV>& grid,TV_INT& influence_corner,VECTOR<TV,IN>& weight)
{  
    static T eps=1e-5;
    TV_INT my_corner=grid.Cell(X,0);
    influence_corner=my_corner;
    for(int n=0;n<IN;n++){
        for(int d=0;d<TV::m;d++){
            T x=grid.one_over_dX(d)*(X(d)-(influence_corner(d)+n)*grid.dX(d)-grid.domain.min_corner(d));
            T abs_x=abs(x);
            if(abs_x<1){
                weight(n)(d)=(T)1.0-abs_x;
                if(abs(weight(n)(d))<eps){
                    weight(n)(d)=0;
                    continue;}
            else{
                weight(n)(d)=0;}}}}
}
//#####################################################################
// Function Build_Weights_Exact
//#####################################################################
template<class TV,int order> void MPM_LINEAR_BASIS<TV,order>::
Build_Weights_Exact(const TV& X,const GRID<TV>& grid,TV_INT& influence_corner,ARRAY<T,TV_INT>& weight)
{
    VECTOR<TV,IN> weight_short;
    Build_Weights(X,grid,influence_corner,weight_short);
    RANGE<TV_INT> range(TV_INT(),TV_INT()+IN);    
    weight.Fill(T(1));
    for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
        for(int d=0;d<TV::m;d++) weight(it.index)*=weight_short(it.index(d))(d);}
}
//#####################################################################
// Function Build_Weights_And_Grad_Weight_Over_Weights
//#####################################################################
template<class TV,int order> void MPM_LINEAR_BASIS<TV,order>::
Build_Weights_And_Grad_Weight_Over_Weights(const TV& X,const GRID<TV>& grid,TV_INT& influence_corner,VECTOR<TV,IN>& weight,VECTOR<TV,IN>& grad_weight_over_weight)
{  
    static T eps=1e-5;
    TV_INT my_corner=grid.Cell(X,0);
    influence_corner=my_corner;
    for(int n=0;n<IN;n++){
        for(int d=0;d<TV::m;d++){
            T x=grid.one_over_dX(d)*(X(d)-(influence_corner(d)+n)*grid.dX(d)-grid.domain.min_corner(d));
            T abs_x=abs(x);
            if(abs_x<1){
                weight(n)(d)=(T)1.0-abs_x;
                if(abs(weight(n)(d))<eps){
                    weight(n)(d)=0;
                    grad_weight_over_weight(n)(d)=0;
                    continue;}
                if(x>0) grad_weight_over_weight(n)(d)=-(T)1.0;
                else grad_weight_over_weight(n)(d)=(T)1.0;
                grad_weight_over_weight(n)(d)*=grid.one_over_dX(d);
                grad_weight_over_weight(n)(d)/=weight(n)(d);}
            else{
                weight(n)(d)=0;
                grad_weight_over_weight(n)(d)=0;}}}
}
//#####################################################################
// Function Build_Weights_And_Grad_Weights_Exact
//#####################################################################
template<class TV,int order> void MPM_LINEAR_BASIS<TV,order>::
Build_Weights_And_Grad_Weights_Exact(const TV& X,const GRID<TV>& grid,TV_INT& influence_corner,ARRAY<T,TV_INT>& weight,ARRAY<TV,TV_INT>& grad_weight)
{
    VECTOR<TV,IN> weight_short;
    VECTOR<TV,IN> grad_weight_over_weight_short;
    Build_Weights_And_Grad_Weight_Over_Weights(X,grid,influence_corner,weight_short,grad_weight_over_weight_short);
    RANGE<TV_INT> range(TV_INT(),TV_INT()+IN);    
    weight.Fill(T(1));
    for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
        for(int d=0;d<TV::m;d++) weight(it.index)*=weight_short(it.index(d))(d);
        TV gw;
        for(int d=0;d<TV::m;d++){
            TV lg;
            lg(d)=grad_weight_over_weight_short(it.index(d))(d);
            gw+=lg;}
        grad_weight(it.index)=gw;
        grad_weight(it.index)*=weight(it.index);}
}
//#####################################################################
template class MPM_LINEAR_BASIS<VECTOR<float,2>,1>;
template class MPM_LINEAR_BASIS<VECTOR<float,3>,1>;
template class MPM_LINEAR_BASIS<VECTOR<double,2>,1>;
template class MPM_LINEAR_BASIS<VECTOR<double,3>,1>;
}

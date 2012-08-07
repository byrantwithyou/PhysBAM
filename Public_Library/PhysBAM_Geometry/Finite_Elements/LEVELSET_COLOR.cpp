//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Finite_Elements/LEVELSET_COLOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/REINITIALIZATION.h>
using namespace PhysBAM;
template<class T> static T Phi_And_Color_Helper(VECTOR<T,1>& phi,VECTOR<int,1>& color,const VECTOR<T,0>&,int& c)
{
    c=color(0);
    return phi(0);
}
template<class T,int d,int e> static T Phi_And_Color_Helper(VECTOR<T,e>& phi,VECTOR<int,e>& color,const VECTOR<T,d>& frac,int& c)
{
    VECTOR<T,e/2> phi_half;
    VECTOR<int,e/2> color_half;
    T a=frac(0);
    for(int i=0;i<e/2;i++){
        int c0=color(2*i),c1=color(2*i+1);
        if(c0==c1){
            phi_half(i)=phi(2*i)*(1-a)+phi(2*i+1)*a;
            color_half(i)=c0;}
        else{
            T ph=phi(2*i)*(1-a)-phi(2*i+1)*a;
            if(ph>=0){
                phi_half(i)=ph;
                color_half(i)=c0;}
            else{
                phi_half(i)=-ph;
                color_half(i)=c1;}}}
    return Phi_And_Color_Helper(phi_half,color_half,frac.Remove_Index(0),c);
}
//#####################################################################
// Function Phi_And_Color
//#####################################################################
template<class TV> typename TV::SCALAR LEVELSET_COLOR<TV>::
Phi_And_Color(const TV& X,int& c) const
{
    const VECTOR<TV_INT,1<<TV::m>& bits=GRID<TV>::Binary_Counts(TV_INT());
    VECTOR<T,1<<TV::m> phi_vector;
    VECTOR<int,1<<TV::m> color_vector;
    TV grid_space_X=(X-grid.domain.min_corner)*grid.one_over_dX-grid.MAC_offset;
    TV_INT clamped_index=clamp(TV_INT(floor(grid_space_X)),TV_INT(),grid.counts-2);
    TV interp_frac_raw=grid_space_X-TV(clamped_index),interp_frac=clamp(interp_frac_raw,TV(),TV()+1);
    for(int i=0;i<1<<TV::m;i++){
        TV_INT ind(clamped_index+bits(i));
        phi_vector(i)=phi(ind);
        color_vector(i)=color(ind);}

    return Phi_And_Color_Helper(phi_vector,color_vector,interp_frac,c);
}
//#####################################################################
// Function Get_Raw_Levelset_For_Color
//#####################################################################
template<class TV> void LEVELSET_COLOR<TV>::
Get_Raw_Levelset_For_Color(ARRAY<T,TV_INT>& color_phi,int c,int ghost) const
{
    color_phi.Resize(grid.Domain_Indices(ghost));
    for(RANGE_ITERATOR<TV::m> it(phi.domain);it.Valid();it.Next()){
        int k=color(it.index);
        T p=phi(it.index);
        if(k!=c) p=-p;
        color_phi(it.index)=p;}
}
//#####################################################################
// Function Get_Levelset_For_Color
//#####################################################################
template<class TV> void LEVELSET_COLOR<TV>::
Get_Levelset_For_Color(ARRAY<T,TV_INT>& color_phi,int c,int ghost) const
{
    Get_Raw_Levelset_For_Color(color_phi,c,ghost);
    FAST_LEVELSET<GRID<TV> > fl(grid,color_phi,ghost);
    Reinitialize(&fl,10*ghost,(T)0,(T)ghost,(T)0,(T).5,3,3);
}
template class LEVELSET_COLOR<VECTOR<float,2> >;
template class LEVELSET_COLOR<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_COLOR<VECTOR<double,2> >;
template class LEVELSET_COLOR<VECTOR<double,3> >;
#endif

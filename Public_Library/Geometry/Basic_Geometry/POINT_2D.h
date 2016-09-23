//#####################################################################
// Copyright 2007, Nipun Kwatra, Jonathan Su
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POINT_2D  
//##################################################################### 
#ifndef __POINT_2D__
#define __POINT_2D__

#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Log/LOG.h>
#include <Core/Utilities/TYPE_UTILITIES.h>
#include <Core/Vectors/VECTOR_2D.h>
#include <Tools/Polynomials/QUADRATIC.h>
namespace PhysBAM{

template<class T>
class POINT_2D:public VECTOR<T,2>
{
    typedef VECTOR<T,2> TV;
public:
    POINT_2D()
    {}

    POINT_2D(const TV& x_input)
        :TV(x_input)
    {}

    bool Edge_Edge_Interaction(const POINT_2D<T>& point,const T interaction_distance,T& distance,VECTOR<T,2>& normal) const
    {normal=*this-point;distance=normal.Magnitude();return distance<=interaction_distance;}

    void Edge_Edge_Interaction_Data(const POINT_2D<T>& point,const TV& v1,const TV& v2,const T& distance,TV& normal,const T small_number) const
    {if(distance > small_number) normal/=distance;
    else normal=(v1-v2).Normalized();}

    template<class T_ARRAY>
    bool Edge_Edge_Interaction(const POINT_2D<T>& point,const INDIRECT_ARRAY<T_ARRAY,VECTOR<int,2>&> V_edges,const T interaction_distance,
        T& distance,VECTOR<T,2>& normal,VECTOR<T,3>& weights,bool allow_negative_weights,const T small_number,const bool exit_early=false) const
    {if(!Edge_Edge_Interaction(point,interaction_distance,distance,normal)) return false;
    if(!exit_early) Edge_Edge_Interaction_Data(point,V_edges(0),V_edges(1),distance,normal,small_number);
    return true;}

    template<class T_ARRAY>
    bool Edge_Edge_Collision(const POINT_2D<T>& point,const INDIRECT_ARRAY<T_ARRAY,VECTOR<int,3>&> V_edges,const T dt,const T collision_thickness,T& collision_time,
        VECTOR<T,2>& normal,VECTOR<T,3>& weights,bool allow_negative_weights,const T small_number=0,const bool exit_early=false) const
    {
        VECTOR<T,2> x2_minus_x1=point-*this,v2_minus_v1=V_edges(1)-V_edges(0);
        QUADRATIC<double> quadratic(v2_minus_v1.Magnitude_Squared(),2*VECTOR<T,2>::Dot_Product(x2_minus_x1,v2_minus_v1),x2_minus_x1.Magnitude_Squared()-sqr(collision_thickness));
        quadratic.Compute_Roots_In_Interval(0,dt);
        if(quadratic.roots==0)return false;
        else if(quadratic.roots==-1){
            LOG::cout<<"VERY SINGULAR ON QUADRATIC SOLVE"<<std::endl;
            collision_time=0;}
        else if(quadratic.roots==1)collision_time=(T)quadratic.root1;
        else collision_time=min((T)quadratic.root1,(T)quadratic.root2);
        POINT_2D<T> new_point(*this+collision_time*V_edges(0));
        T distance;
        return new_point.Edge_Edge_Interaction(POINT_2D<T>(point+collision_time*V_edges(1)),V_edges,collision_thickness,distance,normal,weights,true,small_number,exit_early);
    }

    template<class T_ARRAY>
    bool Edge_Edge_Collision(const POINT_2D<T>& point,const INDIRECT_ARRAY<T_ARRAY,VECTOR<int,2>&> V_edges,const T dt,const T collision_thickness,T& collision_time,
        VECTOR<T,2>& normal,VECTOR<T,3>& weights,const T small_number=0,const bool exit_early=false) const
    {PHYSBAM_NOT_IMPLEMENTED();}
};
}
#endif

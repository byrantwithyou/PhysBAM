//#####################################################################
// Copyright 2004-2009, Ronald Fedkiw, Frank Losasso, Avi Robinson-Mosher, Andy Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROJECTION_DYNAMICS
//#####################################################################
#ifndef __PROJECTION_DYNAMICS__
#define __PROJECTION_DYNAMICS__

#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Data_Structures/TRIPLE.h>
namespace PhysBAM{

template<class T>
class PROJECTION_DYNAMICS
{
public:
    bool flame;

    ARRAY<T> densities;
    // flame speed constants:  gives the reaction speed information between any two phases
    // flame_speed=(.x+.y*curvature)*normal; jump=.z*flame_speed
    // .z is the jump constant (-density_1*(1/density_1-1/density_2))
    // if .z is zero getting the jump is shortcircuited
    // flame_speed_constants is assumed to be symmetric with respect to x and z and skew-symmetric with respect to y
    ARRAY<TRIPLE<T,T,T> ,VECTOR<int,2> > flame_speed_constants;

public:
    PROJECTION_DYNAMICS(const bool flame_input)  
        :flame(flame_input)
    {}
    PROJECTION_DYNAMICS(const PROJECTION_DYNAMICS&) = delete;
    void operator=(const PROJECTION_DYNAMICS&) = delete;

    virtual ~PROJECTION_DYNAMICS()
    {}

    void Set_Densities(const ARRAY<T>& densities_input)
    {densities=densities_input;}
    
    void Set_Flame_Speed_Constants(const ARRAY<T>& densities_input,const ARRAY<T,VECTOR<int,2> >& normal_flame_speeds,const ARRAY<T,VECTOR<int,2> >& curvature_flame_speeds)
    {flame_speed_constants.Resize(VECTOR<int,2>()+densities_input.m,no_init);flame_speed_constants.Fill({});
    densities=densities_input;
    for(int i=0;i<densities_input.m;i++)for(int j=0;j<densities_input.m;j++)if(i!=j&&(normal_flame_speeds(i,j)!=0||curvature_flame_speeds(i,j)!=0)){
        flame_speed_constants(i,j).x=normal_flame_speeds(i,j);assert((normal_flame_speeds(i,j)-normal_flame_speeds(j,i))<(T)1e-6); // must be symmetric
        flame_speed_constants(i,j).y=curvature_flame_speeds(i,j);assert((curvature_flame_speeds(i,j)+curvature_flame_speeds(j,i))<(T)1e-6); // must be skew symmetric
        if(densities(i)>densities(j)) flame_speed_constants(i,j).z=-densities(i)*(1/densities(i)-1/densities(j)); // make sure that the higher density is the reaction source
        else flame_speed_constants(i,j).z=-densities(j)*(1/densities(j)-1/densities(i));}}

//#####################################################################
};
}
#endif


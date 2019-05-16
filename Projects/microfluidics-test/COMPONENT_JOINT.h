//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __COMPONENT_JOINT__
#define __COMPONENT_JOINT__

#include <map>
#include "CANONICAL_BLOCK.h"
#include "CANONICAL_COMPONENT.h"

namespace PhysBAM{

extern double comp_tol;

template<class TV> struct BLOCK_MESHING_ITERATOR;

template<class T>
struct JOINT_KEY
{
    int num_dofs;
    T width;
    ARRAY<T> angles;
    bool operator<(const JOINT_KEY& p) const
    {
        if(num_dofs<p.num_dofs) return true;
        if(p.num_dofs<num_dofs) return false;
        if(width<p.width-comp_tol) return true;
        if(p.width<width-comp_tol) return false;
        if(angles.m!=p.angles.m) return angles.m<p.angles.m;
        for(int i=0;i<angles.m;i++)
        {
            if(angles(i)<p.angles(i)-comp_tol) return true;
            if(p.angles(i)<angles(i)-comp_tol) return false;
        }
        return false;
    }
};

template<class T>
struct COMPONENT_JOINT
{
    typedef VECTOR<T,2> TV;
    
    T target_length;
    std::map<JOINT_KEY<T>,PAIR<CANONICAL_COMPONENT<T>*,ARRAY<T> > > canonical_joints;
    int num_j2=0,num_j3_avg=0,num_j3_small=0,num_j4=0;

    PAIR<CANONICAL_COMPONENT<T>*,ARRAY<T> > Make_Component(int d,T width,const ARRAY<T>& angles);
    PAIR<CANONICAL_COMPONENT<T>*,ARRAY<T> > Make_Joint_2(int d,T width,const ARRAY<T>& angles);
    PAIR<CANONICAL_COMPONENT<T>*,ARRAY<T> > Make_Joint_3(int d,T width,const ARRAY<T>& angles);
    PAIR<CANONICAL_COMPONENT<T>*,ARRAY<T> > Make_Joint_3_Small(int d,T width,const ARRAY<T>& angles);
    PAIR<CANONICAL_COMPONENT<T>*,ARRAY<T> > Make_Joint_3_Average(int d,T width,const ARRAY<T>& angles);
    PAIR<CANONICAL_COMPONENT<T>*,ARRAY<T> > Make_Joint_4_Right_Angle(int d,T width);
    std::tuple<TV,T,T> Elbow_Pit(T angle,T width) const;
    TV Elbow_Pit_Oriented(T angle,T width) const;
    VECTOR<TV,2> Extrude(const TV& v0,const TV& v1,const TV& n) const;
    PAIR<ARRAY<TV>,ARRAY<TV> > Arc(const TV& c,T angle,T len_arm,T ext0,T ext1) const;
    ARRAY<TV> Polyline(const ARRAY<TV>& points,T dx) const;
    void Joint_Connection(int offset,BLOCK_MESHING_ITERATOR<TV>& it,
        CANONICAL_BLOCK<T>* cb,CC_IRREGULAR_CONNECTION* ic0,CC_IRREGULAR_CONNECTION* ic1,
        ARRAY<CC_BLOCK_CONNECTION,CON_ID>& con,CON_ID prev) const;
};

}

#endif

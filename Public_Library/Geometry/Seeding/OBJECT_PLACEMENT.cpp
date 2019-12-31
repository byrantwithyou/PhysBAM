//#####################################################################
// Copyright 2019, Craig Schroeder, Yunxin Sun.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace OBJECT_PLACEMENT
//#####################################################################
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Seeding/OBJECT_PLACEMENT.h>

namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> OBJECT_PLACEMENT<TV>::
OBJECT_PLACEMENT(RANDOM_NUMBERS<T>& random,IMPLICIT_OBJECT<TV>* io)
    :random(random),io(io)
{
    io->Update_Box();
    seed_box=io->Box();
}
//#####################################################################
// Function Seed_Object
//#####################################################################
template<class TV> bool OBJECT_PLACEMENT<TV>::
Seed_Object(SEEDING& obj)
{
    for(int i=0;i<max_intersect_tries;i++)
    {
        if(!Find_Seed_Location(obj.frame,obj.scale,obj.index)) return false;
        if(!Test_Location(obj.frame,obj.scale,obj.index)) continue;

        const auto& o=objects(obj.index);
        obj.twist.angular=random.template Get_Vector_In_Unit_Sphere<typename TV::SPIN>()*o.max_angulare_vel;
        random.Fill_Uniform(obj.twist.linear,o.velocity_range);
        return true;
    }
    return false;
}
//#####################################################################
// Function Test_Location
//#####################################################################
template<class TV> bool OBJECT_PLACEMENT<TV>::
Test_Location(const FRAME<TV>& frame,T scale,int index)
{
    const auto& o=objects(index);
    ORIENTED_BOX<TV> orient_box(o.box*scale,frame);
    RANGE<TV> rot_box=orient_box.Axis_Aligned_Bounding_Box();
    for(auto r:boxes)
        if(rot_box.Intersection(r,min_sep/2))
            return false;
    for(auto p:pts)
        if(orient_box.Thickened(min_sep).Lazy_Inside(p))
            return false;
    boxes.Append(rot_box);
    return true;
}
//#####################################################################
// Function Find_Seed_Location
//#####################################################################
template<class TV> bool OBJECT_PLACEMENT<TV>::
Find_Seed_Location(FRAME<TV>& frame,T& scale,int& index)
{
    for(int j=0;j<max_seed_tries;j++)
    {
        index=random.Get_Uniform_Integer(0,objects.m-1);
        const auto& o=objects(index);
        random.Fill_Uniform(scale,o.scale);
        Random_Fill_Uniform(random,frame.r);
        RANGE<TV> sc_box=o.box*scale;
        RANGE<TV> rot_box=ORIENTED_BOX<TV>(sc_box,frame.r).Axis_Aligned_Bounding_Box();
        RANGE<TV> new_box(seed_box.min_corner-rot_box.min_corner,seed_box.max_corner-rot_box.max_corner);
        if(new_box.Empty()) continue;
        Random_Fill_Uniform(random,frame.t,new_box);
        VECTOR<TV,1<<TV::m> corners;
        sc_box.Corners(corners);
        bool in=true;
        for(auto x:corners)
            if(!io->Lazy_Inside(frame*x)){
                in=false;
                break;}
        if(in) return true;
    }
    return false;
}
//#####################################################################
// Function Seed_Objects
//#####################################################################
template<class TV> void OBJECT_PLACEMENT<TV>::
Seed_Objects(ARRAY<SEEDING>& objs,int max_objects)
{
    SEEDING obj;
    for(int i=0;i<max_objects && Seed_Object(obj);i++)
        objs.Append(obj);
}
template class OBJECT_PLACEMENT<VECTOR<float,1> >;
template class OBJECT_PLACEMENT<VECTOR<float,2> >;
template class OBJECT_PLACEMENT<VECTOR<float,3> >;
template class OBJECT_PLACEMENT<VECTOR<double,1> >;
template class OBJECT_PLACEMENT<VECTOR<double,2> >;
template class OBJECT_PLACEMENT<VECTOR<double,3> >;
}

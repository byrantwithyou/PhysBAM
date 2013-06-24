//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_STANDARD_TESTS
//#####################################################################
#ifndef __RIGIDS_STANDARD_TESTS__
#define __RIGIDS_STANDARD_TESTS__

#include <Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Rigids/Joints/JOINT_ID.h>
namespace PhysBAM{

template<class TV> class RIGID_BODY;
template<class TV> class RIGID_BODY_COLLECTION;

template<class TV>
class RIGIDS_STANDARD_TESTS
{
    typedef typename TV::SCALAR T;
public:
    STREAM_TYPE stream_type;
    std::string data_directory;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;

    RIGIDS_STANDARD_TESTS(STREAM_TYPE stream_type,const std::string& data_directory,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input);

//#####################################################################
    RIGID_BODY<TV>& Add_Rigid_Body(const std::string& name,const T scaling_factor,const T friction,const bool read_implicit=true,const bool always_read_object=false);
    RIGID_BODY<TV>& Add_Ground(const T friction=(T).3,const T height=0,const T coefficient_of_restitution=(T).5,const T scale=(T)1);
    RIGID_BODY<TV>& Add_Analytic_Box(const VECTOR<T,1>& scaling_factor);
    RIGID_BODY<TV>& Add_Analytic_Box(const VECTOR<T,2>& scaling_factor,int segments_per_side=1);
    RIGID_BODY<TV>& Add_Analytic_Box(const VECTOR<T,3>& scaling_factor);
    RIGID_BODY<TV>& Add_Analytic_Torus(const T inner_radius,const T outer_radius,int inner_resolution=8,int outer_resolution=16);
    RIGID_BODY<TV>& Add_Analytic_Cylinder(const T height,const T radius,int resolution_radius=16,int resolution_height=4);
    RIGID_BODY<TV>& Add_Analytic_Shell(const T height,const T outer_radius,const T inner_radius,int resolution=40);
    RIGID_BODY<TV>& Add_Analytic_Bowl(const T height,const T outer_radius,const T inner_radius,int res_radial=32,int res_vertical=8);
    RIGID_BODY<TV>& Add_Analytic_Sphere(const T radius,const T density,int levels=4);
    RIGID_BODY<TV>& Add_Analytic_Smooth_Gear(const TV& dimensions,int cogs,int levels=4);
    void Make_Lathe_Chain(const FRAME<TV>& frame,const T scale=1,const T friction=(T).6,const T cor=(T).3);
    void Set_Joint_Frames(JOINT_ID id,const TV& location);
    JOINT_ID Connect_With_Point_Joint(RIGID_BODY<TV>& parent,RIGID_BODY<TV>& child,const TV& location);
//#####################################################################
};
}
#endif

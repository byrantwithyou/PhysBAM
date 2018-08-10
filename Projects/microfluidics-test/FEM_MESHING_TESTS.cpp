//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include "FEM_MESHING_TESTS.h"
#include "FLUID_LAYOUT_FEM.h"

namespace PhysBAM{
//#####################################################################
// Function Test_Degree2_Joint
//#####################################################################
template<typename TV>
void Test_Degree2_Joint(JOINT_TYPE jt,typename TV::SCALAR a0,typename TV::SCALAR a1,typename TV::SCALAR da)
{
    typedef typename TV::SCALAR T;

    for(T rad=a0;rad<a1;rad+=da){
        PARSE_DATA_FEM<TV> pd;
        pd.half_width=4;
        pd.unit_length=0.5;
        pd.pts.Append({TV(-10,1),dirichlet_v,default_joint});
        pd.pts.Append({TV(),nobc,jt});
        TV p=TV(10,0);
        p={cos(rad)*p(0)-sin(rad)*p(1),sin(rad)*p(0)+cos(rad)*p(1)};
        pd.pts.Append({p,traction,default_joint});
        pd.pipes.Append({0,1});
        pd.pipes.Append({1,2});
        pd.joints.Set(0,{0});
        pd.joints.Set(1,{0,1});
        pd.joints.Set(2,{1});

        FLUID_LAYOUT_FEM<TV> fl;
        fl.Dump_Input(pd);
        Flush_Frame<TV>("init");
        fl.Compute(pd);
        fl.Print_Statistics();
        fl.Dump_Layout();
        Flush_Frame<TV>("blocks");}
}
//#####################################################################
// Function Test_Degree2_Circle
//#####################################################################
template<typename TV>
void Test_Degree2_Circle(JOINT_TYPE jt,typename TV::SCALAR h0,typename TV::SCALAR h1,typename TV::SCALAR dh)
{
    typedef typename TV::SCALAR T;

    for(T h=h0;h<h1;h+=dh){
        PARSE_DATA_FEM<TV> pd;
        pd.half_width=4;
        pd.unit_length=h;
        pd.pts.Append({TV(-4,-4),nobc,jt});
        pd.pts.Append({TV(4,-4),nobc,jt});
        pd.pts.Append({TV(4,4),nobc,jt});
        pd.pts.Append({TV(-4,4),nobc,jt});
        pd.pipes.Append({0,1});
        pd.pipes.Append({1,2});
        pd.pipes.Append({3,2});
        pd.pipes.Append({3,0});
        pd.joints.Set(0,{0,3});
        pd.joints.Set(1,{0,1});
        pd.joints.Set(2,{1,2});
        pd.joints.Set(3,{3,2});

        FLUID_LAYOUT_FEM<TV> fl;
        fl.Dump_Input(pd);
        Flush_Frame<TV>("init");
        fl.Compute(pd);
        fl.Print_Statistics();
        fl.Dump_Layout();
        Flush_Frame<TV>("blocks");}
}
//#####################################################################
// Function Test_Degree3_Joint
//#####################################################################
template<typename TV>
void Test_Degree3_Joint(JOINT_TYPE jt,typename TV::SCALAR h,int n,int seed)
{
    typedef typename TV::SCALAR T;

    RANDOM_NUMBERS<T> random(seed);
    for(int i=0;i<n;i++){
        PARSE_DATA_FEM<TV> pd;
        pd.half_width=4;
        pd.unit_length=h;
        pd.pts.Append({TV(),nobc,jt});
        for(int k=0;k<3;k++){
            T a=random.Get_Uniform_Number(0,2*pi);
            TV coord=TV(cos(a),sin(a))*5;
            pd.pts.Append({coord,traction,jt});
            pd.pipes.Append({0,k+1});
            pd.joints.Set(k+1,{k});}
        pd.joints.Set(0,{0,1,2});
        FLUID_LAYOUT_FEM<TV> fl;
        fl.Dump_Input(pd);
        Flush_Frame<TV>("init");
        fl.Compute(pd);
        fl.Print_Statistics();
        fl.Dump_Layout();
        Flush_Frame<TV>("blocks");}
}

template void Test_Degree2_Joint<VECTOR<double,2> >(JOINT_TYPE,double,double,double);
template void Test_Degree2_Circle<VECTOR<double,2> >(JOINT_TYPE,double,double,double);
template void Test_Degree3_Joint<VECTOR<double,2> >(JOINT_TYPE,double,int,int);
}

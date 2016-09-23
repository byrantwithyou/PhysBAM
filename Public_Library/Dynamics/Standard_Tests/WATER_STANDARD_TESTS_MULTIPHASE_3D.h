//#####################################################################
// Copyright 2005, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class WATER_STANDARD_TESTS_MULTIPHASE_3D
//#####################################################################
// Test descriptions:
//   11. Two drops colliding
//   11. 4 phase splash
// Also supports a variety of standard resolutions in powers of 2.
//#####################################################################
#ifndef __WATER_STANDARD_TESTS_MULTIPHASE_3D__
#define __WATER_STANDARD_TESTS_MULTIPHASE_3D__

#include <Core/Random_Numbers/NOISE.h>
#include <Dynamics/Standard_Tests/WATER_STANDARD_TESTS_3D.h>
#include <Dynamics/Standard_Tests/WATER_STANDARD_TESTS_MULTIPHASE.h>
namespace PhysBAM{

template<class TV>
class WATER_STANDARD_TESTS_MULTIPHASE_3D:public WATER_STANDARD_TESTS_MULTIPHASE<TV,WATER_STANDARD_TESTS_3D<TV> >
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,3> TV_INT;
public:
    typedef WATER_STANDARD_TESTS_MULTIPHASE<TV,WATER_STANDARD_TESTS_3D<TV> > BASE;
    using BASE::rigid_body_collection;using BASE::fluids_parameters;using BASE::grid;using BASE::example;using BASE::sphere;using BASE::test_number;using BASE::world_to_source;
    using BASE::source_velocity;using BASE::source_region;using BASE::sources;

    WATER_STANDARD_TESTS_MULTIPHASE_3D(SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>& example_input,FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters_input,FLUID_COLLECTION<TV>& fluid_collection,
        RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
        :WATER_STANDARD_TESTS_MULTIPHASE<TV,WATER_STANDARD_TESTS_3D<TV> >(example_input,fluids_parameters_input,fluid_collection,rigid_body_collection_input)
    {
    }

void Initialize(const int test_number_input,const int resolution,const int restart_frame)
{
    BASE::Initialize(test_number_input,resolution,restart_frame);
    test_number=test_number_input;
    LOG::cout<<"Running Multiphase Standard Test Number "<<test_number<<" at resolution "<<resolution<<std::endl;
    int cells=1*resolution;
    if(test_number==11){
        fluids_parameters.domain_walls[2][2]=true;
        example.first_frame=0;example.last_frame=500;example.frame_rate=250;
        grid.Initialize(TV_INT(10*cells+1,10*cells+1,10*cells+1),RANGE<TV>(TV(0,0,0),TV((T).1,(T).1,(T).1)));}
    if(test_number==12){
        fluids_parameters.domain_walls[2][2]=false;
        grid.Initialize(TV_INT(15*cells+1,15*cells+1,10*cells+1),RANGE<TV>(TV(0,-.5,0),TV(1.5,1,1)));}
    if(test_number==13){
        fluids_parameters.domain_walls[2][2]=true;
        example.first_frame=0;example.last_frame=500;example.frame_rate=100;
        grid.Initialize(TV_INT(10*cells+1,10*cells+1,10*cells+1),RANGE<TV>(TV(0,0,0),TV(2,2,2)));}
    if(test_number==14){
        fluids_parameters.domain_walls[2][2]=true;
        example.first_frame=0;example.last_frame=20;example.frame_rate=400;
        grid.Initialize(TV_INT(10*cells,20*cells,10*cells),RANGE<TV>(TV((T)-.01,(T)-.01,(T)-.01),TV((T).01,(T).02,(T).01)));}
    if(test_number==15){
        example.first_frame=0;example.last_frame=500;example.frame_rate=60;
        grid.Initialize(TV_INT(20*cells+1,10*cells+1,16*cells+1),RANGE<TV>(TV(0,0,(T).1),TV(1,(T).5,(T).9)));
        world_to_source.Append(MATRIX<T,4>::Identity_Matrix());
        sources.Append(CYLINDER<T>(TV((T).8,(T).45,(T).5),TV((T).8,(T).6,(T).5),(T).075));
        source_velocity.Append(TV(0,(T)-.6,0));
        source_region.Append(3);
        world_to_source.Append(MATRIX<T,4>::Identity_Matrix());
        sources.Append(CYLINDER<T>(TV((T).625,(T).45,(T).5),TV((T).625,(T).6,(T).5),(T).075));
        source_velocity.Append(TV());
        source_region.Append(4);}
    if(test_number==16){
        example.first_frame=0;example.last_frame=1000;example.frame_rate=120;
        grid.Initialize(TV_INT(10*cells+1,11*cells+1,10*cells+1),RANGE<TV>(TV(0,0,0),TV(1,(T)1.1,1)));}
    if(test_number==17){
        example.first_frame=0;example.last_frame=500;example.frame_rate=250;
        grid.Initialize(TV_INT(10*cells+1,10*cells+1,10*cells+1),RANGE<TV>(TV(0,0,0),TV((T).1,(T).1,(T).1)));}
    example.output_directory=LOG::sprintf("Standard_Tests_Multiphase/Test_%d__Resolution_%d_%d_%d",test_number,(grid.counts.x-1),(grid.counts.y-1),(grid.counts.z-1));
    LOG::cout<<"output directory="<<example.output_directory<<std::endl;
}

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    WATER_STANDARD_TESTS_MULTIPHASE<TV,WATER_STANDARD_TESTS_3D<TV> >::Initialize_Bodies();
    if(test_number==15){
        int ground=rigid_body_collection.Add_Rigid_Body(example.stream_type,example.data_directory+"/Rigid_Bodies/ground",(T).01,true,true,false);
        rigid_body_collection.rigid_body_particles.frame(ground).t=TV((T).1,0,0);
        rigid_body_collection.rigid_body_particles.frame(ground).r=ROTATION<TV>::From_Rotation_Vector(TV(0,0,(T)pi/8));
        rigid_body_collection.rigid_body_particles.kinematic(ground)=true;
        fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection);}
}
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2005, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class WATER_STANDARD_TESTS_MULTIPHASE_2D
//#####################################################################
// Test descriptions:
//   11. Two drops colliding
//   12. Sphere splashing into a pool of two liquid phases
//   13. 4 phase splash
//   14. Rising air bubble in water
// Also supports a variety of standard resolutions in powers of 2.
//#####################################################################
#ifndef __WATER_STANDARD_TESTS_MULTIPHASE_2D__
#define __WATER_STANDARD_TESTS_MULTIPHASE_2D__

#include <Core/Random_Numbers/NOISE.h>
#include <Dynamics/Standard_Tests/WATER_STANDARD_TESTS_2D.h>
#include <Dynamics/Standard_Tests/WATER_STANDARD_TESTS_MULTIPHASE.h>
namespace PhysBAM{

template<class TV>
class WATER_STANDARD_TESTS_MULTIPHASE_2D:public WATER_STANDARD_TESTS_MULTIPHASE<TV,WATER_STANDARD_TESTS_2D<TV> >
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,2> TV_INT;
public:
    typedef WATER_STANDARD_TESTS_MULTIPHASE<TV,WATER_STANDARD_TESTS_2D<TV> > BASE;
    using BASE::rigid_body_collection;using BASE::fluids_parameters;using BASE::grid;using BASE::example;using BASE::sphere;using BASE::test_number;using BASE::world_to_source;
    using BASE::source_velocity;using BASE::source_region;using BASE::sources;

    WATER_STANDARD_TESTS_MULTIPHASE_2D(SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>& example_input,FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters_input,FLUID_COLLECTION<TV>& fluid_collection,
        RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
        :WATER_STANDARD_TESTS_MULTIPHASE<TV,WATER_STANDARD_TESTS_2D<TV> >(example_input,fluids_parameters_input,fluid_collection,rigid_body_collection_input)
    {
    }

void Initialize(const int test_number_input,const int resolution,const int restart_frame=-1)
{
    BASE::Initialize(test_number_input,resolution,restart_frame);
    LOG::cout<<"Running Multiphase Standard Test Number "<<test_number<<" at resolution "<<resolution<<std::endl;
    int cells=1*resolution;
    if(test_number==11){
        fluids_parameters.domain_walls[2][2]=true;
        if(!example.user_last_frame) example.last_frame=500;
        if(!example.user_frame_rate) example.frame_rate=250;
        grid.Initialize(TV_INT(10*cells+1,10*cells+1),RANGE<TV>(TV(0,0),TV((T).1,(T).1)));}
    if(test_number==12){
        fluids_parameters.domain_walls[2][2]=false;
        grid.Initialize(TV_INT(20*cells+1,15*cells+1),RANGE<TV>(TV((T)-.5,(T)-.5),TV((T)1.5,1)));}
    if(test_number==13){
        fluids_parameters.domain_walls[2][2]=true;
        if(!example.user_last_frame) example.last_frame=500;
        if(!example.user_frame_rate) example.frame_rate=100;
        grid.Initialize(TV_INT(10*cells+1,10*cells+1),RANGE<TV>(TV(0,0),TV(2,2)));}
    if(test_number==14){
        fluids_parameters.domain_walls[2][2]=true;
        if(!example.user_last_frame) example.last_frame=20;
        if(!example.user_frame_rate) example.frame_rate=400;
        grid.Initialize(TV_INT(20*cells+1,30*cells+1),RANGE<TV>(TV((T)-.01,(T)-.01),TV((T).01,(T).02)));}
    if(test_number==15){
        fluids_parameters.domain_walls[2][2]=true;
        if(!example.user_last_frame) example.last_frame=500;
        if(!example.user_frame_rate) example.frame_rate=60;
        grid.Initialize(TV_INT(20*cells+1,10*cells+1),RANGE<TV>(TV(0,0),TV(1,(T).5)));
        world_to_source.Append(MATRIX<T,3>::Identity_Matrix());
        sources.Append(RANGE<TV>(TV((T).75,(T).45),TV((T).875,(T).6)));
        source_velocity.Append(TV(0,(T)-.6));
        source_region.Append(3);
        world_to_source.Append(MATRIX<T,3>::Identity_Matrix());
        sources.Append(RANGE<TV>(TV((T).6,(T).45),TV((T).725,(T).6)));
        source_velocity.Append(TV(0,0));
        source_region.Append(4);}
    if(test_number==16){
        if(!example.user_last_frame) example.last_frame=500;
        if(!example.user_frame_rate) example.frame_rate=120;
        grid.Initialize(TV_INT(10*cells+1,11*cells+1),RANGE<TV>(TV(0,0),TV(1,(T)1.1)));}
    if(test_number==17){
        if(!example.user_last_frame) example.last_frame=500;
        if(!example.user_frame_rate) example.frame_rate=250;
        grid.Initialize(TV_INT(10*cells+1,10*cells+1),RANGE<TV>(TV(0,0),TV((T).1,(T).1)));}
    if(!example.user_output_directory)
        example.viewer_dir.output_directory=LOG::sprintf("Standard_Tests_Multiphase/Test_%d__Resolution_%d_%d",test_number,(grid.counts.x-1),(grid.counts.y-1));
    LOG::cout<<"output directory="<<example.viewer_dir.output_directory<<std::endl;
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    WATER_STANDARD_TESTS_MULTIPHASE<TV,WATER_STANDARD_TESTS_2D<TV> >::Initialize_Bodies();
    if(test_number==15){
        int ground=rigid_body_collection.Add_Rigid_Body(example.data_directory+"/Rigid_Bodies_2D/ground",(T).1,true,true,false);
        rigid_body_collection.rigid_body_particles.frame(ground).t=TV((T).1,0);
        rigid_body_collection.rigid_body_particles.frame(ground).r=ROTATION<TV>::From_Angle((T)pi/8);
        rigid_body_collection.rigid_body_particles.kinematic(ground)=true;}
    fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection);
}
//#####################################################################
};
}
#endif

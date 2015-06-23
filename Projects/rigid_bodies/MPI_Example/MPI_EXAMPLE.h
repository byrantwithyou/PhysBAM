//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_EXAMPLE
//#####################################################################
//   1. Pyramid of boxes
//   2. Stacked boxes
//   3. Partition Test
//   4. Contact Test 1
//   5. Contact Test 2
//   6. Simple Collision Test
//   7. Collision Test
//   8. Pushout Test
//   9. 320k Analytic Spheres
//#####################################################################
#ifndef __MPI_EXAMPLE__
#define __MPI_EXAMPLE__

#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Tessellation/RANGE_TESSELLATION.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <Rigids/Forces_And_Torques/RIGID_GRAVITY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE.h>
#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include "../../rigid_bodies/RANDOM_PLACEMENT.h"
namespace PhysBAM{

template<class T_input>
class MPI_EXAMPLE:public SOLIDS_EXAMPLE<VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T,3> TV;
public:
    int width, height, num_bodies;

    typedef SOLIDS_EXAMPLE<TV> BASE;
    using BASE::solids_parameters;using BASE::solid_body_collection;using BASE::solids_evolution;using BASE::test_number;using BASE::frame_rate;
    using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::stream_type;

    SOLIDS_STANDARD_TESTS<TV> tests;

    int kinematic_body_id;
    INTERPOLATION_CURVE<T,FRAME<TV> > curve;

    MPI_EXAMPLE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :BASE(stream_type_input,parse_args),width(2),height(4),num_bodies(6),tests(stream_type_input,data_directory,solid_body_collection)
    {
        LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.rigid_body_collision_parameters.use_legacy_push_out=true;
        //solids_parameters.rigid_body_collision_parameters.use_shock_propagation=false;
        solids_parameters.cfl=1;
        parse_args.Add("-width",&width,"value","number of stacks");
        parse_args.Add("-height",&height,"value","height of each stack");
        parse_args.Add("-num_bodies",&num_bodies,"value","number of total bodies");
        parse_args.Parse();
        tests.data_directory=data_directory;
        output_directory=LOG::sprintf("MPI_Example/Test_%d",test_number);
    }

    ~MPI_EXAMPLE()
    {}

    // Unused callbacks
    void Postprocess_Solids_Substep(const T time,const int substep) override {}
    void Apply_Constraints(const T dt,const T time) override {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() override {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) override {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) override {}
    void Preprocess_Solids_Substep(const T time,const int substep) override {}
    void Limit_Solids_Dt(T& dt,const T time) override {}
    void Update_Time_Varying_Material_Properties(const T time) override {}
    void Preprocess_Frame(const int frame) override {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) override {}
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) override {}
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) override {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) override {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) override {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) override {}
    void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt) override {}
    void Postprocess_Substep(const T dt,const T time) override {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) override {}
    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) override {return false;}
    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) override {}
    void Preprocess_Substep(const T dt,const T time) override {}

//#####################################################################
// Function Post_Initialization
//#####################################################################
void Post_Initialization() override
{
    RIGID_BODY_COLLISIONS<TV>& collisions=*solids_evolution->rigid_body_collisions;
    collisions.Set_Push_Out_Level_Iterations(1);
    if(test_number==9||test_number==10){
        int level=1;
        solids_parameters.rigid_body_collision_parameters.collision_iterations=level;
        collisions.Set_Collision_Pair_Iterations(level*2);
        solids_parameters.rigid_body_collision_parameters.contact_iterations=level;
        collisions.Set_Contact_Level_Iterations(level);
        collisions.Set_Contact_Pair_Iterations(level*2);
        solids_parameters.rigid_body_collision_parameters.use_shock_propagation=true;
        collisions.Set_Shock_Propagation_Iterations(level);
        collisions.Set_Shock_Propagation_Level_Iterations(level);
        collisions.Set_Shock_Propagation_Pair_Iterations(level*2);
        solids_parameters.rigid_body_collision_parameters.use_push_out=true;
        solids_parameters.rigid_body_collision_parameters.use_legacy_push_out=true;
        collisions.Set_Push_Out_Iterations(level);
        collisions.Set_Push_Out_Level_Iterations(level);
        collisions.Set_Push_Out_Pair_Iterations(level*2);
        if(test_number==9) collisions.Register_Analytic_Collisions();}
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) override
{
    RIGID_BODY_COLLISIONS<TV>& collisions=*solids_evolution->rigid_body_collisions;
    if(test_number==8){
        solids_parameters.rigid_body_collision_parameters.perform_collisions=false;
        solids_parameters.rigid_body_collision_parameters.collision_iterations=0;
        solids_parameters.rigid_body_collision_parameters.contact_iterations=0;
        collisions.Set_Shock_Propagation_Iterations(0);}
}
void After_Initialization() override {BASE::After_Initialization();}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{
    switch(test_number){
        case 1: Pyramid_Of_Boxes(); break;
        case 2: Stacked_Boxes(); break;
        case 3: Partition_Test(); break;
        case 4: Contact_Test_1(); break;
        case 5: Contact_Test_2(); break;
        case 6: Simple_Collision_Test(); break;
        case 7: Collision_Test(); break;
        case 8: Pushout_Test(); break;
        case 9: Many_Sphere_Test(); break;
        case 10: Break_Levelset(); break;
        default: PHYSBAM_FATAL_ERROR(LOG::sprintf("Unrecognized test number %d",test_number));}
}
//#####################################################################
// Function Pyramid_Of_Boxes
//#####################################################################
void Pyramid_Of_Boxes()
{
    int height = 21;

    T current_x;
    T first_x = 0;
    for (int i = 1; i < height; i++) {
        current_x = first_x;
        for (int j =  1; j <= i ; j++) {
            RIGID_BODY<TV>* rigid_body = &tests.Add_Rigid_Body("subdivided_box",(T)2, (T).1);
            rigid_body->Frame().t=TV((T)current_x, (height+1-i)*5+80,0);
            rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
            current_x += 5;}
        first_x -= 2.5;}

    for(int i=0;i<2;i++){
        RIGID_BODY<TV>* rigid_body0=&tests.Add_Rigid_Body("subdivided_box",(T)10,(T).1);
        rigid_body0->Frame().t=TV(first_x-10,i*20,0);
        rigid_body0->is_static = true;
        rigid_body0->name="left_box";

        RIGID_BODY<TV>* rigid_body1=&tests.Add_Rigid_Body("subdivided_box",(T)10,(T).1);
        rigid_body1->Frame().t=TV(current_x+7.5,i*20,0);
        rigid_body1->is_static = true;
        rigid_body1->name="right_box";}

    tests.Add_Ground(1, -10);
    last_frame = 400;
    solid_body_collection.rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection, true));
}
//#####################################################################
// Function Stacked_Boxes
//#####################################################################
void Stacked_Boxes() {
    int height = 50;

    for (int i = 1; i <= height; i++) {
        for (int j =  1; j < 20 ; j++) {
            RIGID_BODY<TV>* rigid_body = &tests.Add_Rigid_Body("subdivided_box",(T)2, (T).1);
            rigid_body->Frame().t=TV(5*j - 50, 5*i+10,0);
            rigid_body->Set_Coefficient_Of_Restitution((T)0.5);}}

    for (int i = 0; i < height/4; i++) {
        RIGID_BODY<TV>* rigid_body0=&tests.Add_Rigid_Body("subdivided_box",(T)10,(T).1);
        rigid_body0->Frame().t=TV(-58,i*20,0);
        rigid_body0->is_static = true;
        rigid_body0->name="left_box";
    
        RIGID_BODY<TV>* rigid_body1=&tests.Add_Rigid_Body("subdivided_box",(T)10,(T).1);
        rigid_body1->Frame().t=TV(58,i*20,0);
        rigid_body1->is_static = true;
        rigid_body1->name="right_box";}

    tests.Add_Ground(1, -10);
    last_frame = 250;
    solid_body_collection.rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection, true));
}
//#####################################################################
// Function Partition_Test
//#####################################################################
void Partition_Test() {
    RIGID_BODY<TV>* rigid_body = &tests.Add_Rigid_Body("box", (T)1, (T).1);
    rigid_body->Frame().t = TV(0,0,0);
    rigid_body->Twist().linear = TV(1,0,0);

    RIGID_BODY<TV>* rigid_body1 = &tests.Add_Rigid_Body("box", (T)1, (T).1);
    rigid_body1->Frame().t = TV(10,3,0);
    rigid_body1->Twist().linear = TV(-1,0,0);

    last_frame = 200;
}
//#####################################################################
// Function Contact_Test_1
//#####################################################################
void Contact_Test_1() {
    for (int i=0;i<width;i++) {
        for (int j=0;j<height;j++) {
            RIGID_BODY<TV>* rigid_body = &tests.Add_Rigid_Body("subdivided_box",(T)1,(T).1);
            rigid_body->Frame().t = TV(i*4,j*2,0);}}

    tests.Add_Ground(1, -1);
    solid_body_collection.rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection, true));
}
//#####################################################################
// Function Contact_Test_2
//#####################################################################
void Contact_Test_2() {
    if (height>width)
        height=width;

    for (int i=0;i<height;i++) {
        for (int j=0;j<width;j++) {
            RIGID_BODY<TV>* rigid_body = &tests.Add_Rigid_Body("subdivided_box",(T)1,(T).1);
            rigid_body->Frame().t = TV((T)j*2+i*0.5,i*2,0);}
        width--;}

    tests.Add_Ground(1, -1);
    solid_body_collection.rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection, true));
}
//#####################################################################
// Function Simple_Collision_Test
//#####################################################################
void Simple_Collision_Test() {
    RIGID_BODY<TV>* rigid_body = &tests.Add_Rigid_Body("box", (T)1, (T).1);
    rigid_body->Frame().t = TV(0,0,0);
    rigid_body->Twist().linear = TV(2,0,0);

    RIGID_BODY<TV>* rigid_body1 = &tests.Add_Rigid_Body("box", (T)1, (T).1);
    rigid_body1->Frame().t = TV(10,3,0);
    rigid_body1->Twist().linear = TV(-1,0,0);

    RIGID_BODY<TV>* rigid_body2 = &tests.Add_Rigid_Body("box", (T)1, (T).1);
    rigid_body2->Frame().t = TV(10,0,0);
    rigid_body2->Twist().linear = TV(-1,0,0);

    last_frame = 200;
}
//#####################################################################
// Function Collision_Test
//#####################################################################
void Collision_Test() {
    for (int i=0;i<height;i++) {
        RIGID_BODY<TV>* left_body = &tests.Add_Rigid_Body("subdivided_box",(T)1,(T).1);
        RIGID_BODY<TV>* right_body = &tests.Add_Rigid_Body("subdivided_box",(T)1,(T).1);
        left_body->Frame().t = TV(0,i*3,0);
        right_body->Frame().t = TV((T)(width+1)*2.15,i*3,0);
        left_body->Set_Coefficient_Of_Restitution((T)1.0);
        right_body->Set_Coefficient_Of_Restitution((T)1.0);
        left_body->is_static=true;
        right_body->is_static=true;

        for (int j=1;j<=width;j++) {
            RIGID_BODY<TV>* rigid_body = &tests.Add_Rigid_Body("subdivided_box",(T)1,(T).1);
            rigid_body->Frame().t = TV((T)j*2.15,i*3,0);
            rigid_body->Twist().linear = TV(4*(.5-j%2),0,0);
            rigid_body->Set_Coefficient_Of_Restitution((T)1.0);}}

    last_frame = 300;
}
//#####################################################################
// Function Pushout_Test
//#####################################################################
void Pushout_Test() {
    for (int i=0;i<num_bodies;i++) {
        RIGID_BODY<TV>* rigid_body =  &tests.Add_Rigid_Body("subdivided_box",(T)1,(T).1);
        rigid_body->Frame().t = TV((T)1.75*i,(T)1.75*i,0);}
}
//#####################################################################
// Function Analytic_Contact
//#####################################################################
void Many_Sphere_Test()
{
    last_frame=300;
     VECTOR<T,3> num_bodies=TV(40,200,40);
     for(int i=0;i<num_bodies.x;i++) for(int j=0;j<num_bodies.y;j++) for(int k=0;k<num_bodies.z;k++){
        RIGID_BODY_PARTICLES<TV>& particles=solid_body_collection.rigid_body_collection.rigid_body_particles;    
        RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(solid_body_collection.rigid_body_collection);
        T radius=1;TV center;
        T offset=.1,offset2=.1;if(j%2==0){offset=-.1;offset2=-.1;}
        center=VECTOR<T,3>(3*i+offset,3*j+1,3*k+offset2);
        particles.frame(rigid_body.particle_index)=FRAME<TV>(center);
        particles.twist(rigid_body.particle_index)=TWIST<TV>();
        particles.mass(rigid_body.particle_index)=1;
        particles.inertia_tensor(rigid_body.particle_index)=DIAGONAL_MATRIX<T,TV::SPIN::m>();particles.inertia_tensor(rigid_body.particle_index)+=1;
        SPHERE<TV> sphere(TV(0,0,0),radius);
        rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >(sphere));
        solid_body_collection.rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
        solid_body_collection.rigid_body_collection.collision_body_list->Add_Body(new RIGID_COLLISION_GEOMETRY<TV>(rigid_body),rigid_body.particle_index,true);
        solid_body_collection.rigid_body_collection.collision_body_list->Get_Collision_Geometry(rigid_body.particle_index)->add_to_spatial_partition=true;}
    tests.Add_Ground((T).5,0,0);
    RIGID_GRAVITY<TV> *gravity=new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,true);
    solid_body_collection.rigid_body_collection.Add_Force(gravity);
}
//#####################################################################
// Function Cluster_Fracture
//#####################################################################
void Break_Levelset()
{
    typedef VECTOR<int,TV::dimension> TV_INT;
    ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> children;
    ARRAY<int>* referenced_rigid_particles=new ARRAY<int>;
    RIGID_BODY_PARTICLES<TV>& particles=solid_body_collection.rigid_body_collection.rigid_body_particles;    
    GRID<TV>& grid=*new GRID<TV>;
    ARRAY<T,VECTOR<int,3> >& phi=*new ARRAY<T,VECTOR<int,3> >;
    LEVELSET<TV> levelset(grid,phi);
    data_directory="../../../Public_Data/";
    FILE_UTILITIES::Read_From_File(stream_type,LOG::sprintf("%s/Rigid_Bodies/sphere_66k.phi",data_directory.c_str()),levelset);
    TV_INT counts=TV_INT::All_Ones_Vector()*25;
    GRID<TV> body_grid(counts,levelset.grid.Domain());
    TRIANGULATED_SURFACE<T>* surface=0;    
    for(CELL_ITERATOR<TV> iterator(body_grid);iterator.Valid();iterator.Next()){
        TV_INT cell=levelset.grid.Cell(iterator.Location(),0);
        if(!phi.Valid_Index(cell)) continue;
        if(phi(cell)>0) continue;
        RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(solid_body_collection.rigid_body_collection);
        children.Append(rigid_body.particle_index);
        referenced_rigid_particles->Append(rigid_body.particle_index);
        particles.frame(rigid_body.particle_index)=FRAME<TV>(iterator.Location());
        particles.twist(rigid_body.particle_index)=TWIST<TV>();
        particles.mass(rigid_body.particle_index)=1;
        particles.inertia_tensor(rigid_body.particle_index)=DIAGONAL_MATRIX<T,TV::SPIN::m>();particles.inertia_tensor(rigid_body.particle_index)+=1;
        RANGE<TV> box(body_grid.dX*-.5,body_grid.dX*.5);
        //SPHERE<TV> box(TV(),body_grid.dX.Min()*.5);
        rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(box));
        //rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >(box));
        if(!surface) surface=TESSELLATION::Generate_Triangles(box);
        rigid_body.Add_Structure(*surface);  
        solid_body_collection.rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
        solid_body_collection.rigid_body_collection.collision_body_list->Add_Body(new RIGID_COLLISION_GEOMETRY<TV>(rigid_body),rigid_body.particle_index,true);
        solid_body_collection.rigid_body_collection.collision_body_list->Get_Collision_Geometry(rigid_body.particle_index)->add_to_spatial_partition=true;}
    tests.Add_Ground((T).5,-2,0);
    
    solid_body_collection.rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,referenced_rigid_particles));
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) override
{
    solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Positions();
}
//#####################################################################
};
}
#endif

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

#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_FRAME.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <PhysBAM_Geometry/Tessellation/RANGE_TESSELLATION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGID_GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/RIGIDS_EXAMPLE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Standard_Tests/RIGIDS_STANDARD_TESTS.h>
#include "../../rigid_bodies/RANDOM_PLACEMENT.h"
namespace PhysBAM{

template<class T_input>
class MPI_EXAMPLE:public RIGIDS_EXAMPLE<VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T,3> TV;
public:
    int width, height, num_bodies;

    typedef RIGIDS_EXAMPLE<TV> BASE;
    using BASE::rigids_parameters;using BASE::rigid_body_collection;using BASE::rigids_evolution;using BASE::test_number;using BASE::frame_rate;
    using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::stream_type;using BASE::parse_args;using BASE::processes_per_dimension;

    RIGIDS_STANDARD_TESTS<TV> tests;

    int kinematic_body_id;
    INTERPOLATION_CURVE<T,FRAME<TV> > curve;

    MPI_EXAMPLE(const STREAM_TYPE stream_type)
        :BASE(stream_type),tests(*this,rigid_body_collection)
    {
        LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
        rigids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        rigids_parameters.rigid_body_collision_parameters.use_legacy_push_out=true;
        //rigids_parameters.rigid_body_collision_parameters.use_shock_propagation=false;
        rigids_parameters.cfl=1;
    }

    ~MPI_EXAMPLE()
    {}

    // Unused callbacks
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return false;}
    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE {}
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Post_Initialization
//#####################################################################
void Post_Initialization() PHYSBAM_OVERRIDE
{
    RIGID_BODY_COLLISIONS<TV>& collisions=*rigids_evolution->rigid_body_collisions;
    collisions.Set_Push_Out_Level_Iterations(1);
    if(test_number==9||test_number==10){
        int level=1;
        rigids_parameters.rigid_body_collision_parameters.collision_iterations=level;
        collisions.Set_Collision_Pair_Iterations(level*2);
        rigids_parameters.rigid_body_collision_parameters.contact_iterations=level;
        collisions.Set_Contact_Level_Iterations(level);
        collisions.Set_Contact_Pair_Iterations(level*2);
        rigids_parameters.rigid_body_collision_parameters.use_shock_propagation=true;
        collisions.Set_Shock_Propagation_Iterations(level);
        collisions.Set_Shock_Propagation_Level_Iterations(level);
        collisions.Set_Shock_Propagation_Pair_Iterations(level*2);
        rigids_parameters.rigid_body_collision_parameters.use_push_out=true;
        rigids_parameters.rigid_body_collision_parameters.use_legacy_push_out=true;
        collisions.Set_Push_Out_Iterations(level);
        collisions.Set_Push_Out_Level_Iterations(level);
        collisions.Set_Push_Out_Pair_Iterations(level*2);
        if(test_number==9) collisions.Register_Analytic_Collisions();}
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
    RIGID_BODY_COLLISIONS<TV>& collisions=rigid_body_collection.rigid_body_collisions;
    if(test_number==8){
        rigids_parameters.rigid_body_collision_parameters.perform_collisions=false;
        rigids_parameters.rigid_body_collision_parameters.collision_iterations=0;
        rigids_parameters.rigid_body_collision_parameters.contact_iterations=0;
        collisions.Set_Shock_Propagation_Iterations(0);}
}
//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Integer_Argument("-xprocs",0,"procs in the x direction");
    parse_args->Add_Integer_Argument("-yprocs",0,"procs in the y direction");
    parse_args->Add_Integer_Argument("-zprocs",0,"procs in the y direction");
    parse_args->Add_Integer_Argument("-width",2,"number of stacks");
    parse_args->Add_Integer_Argument("-height",4,"height of each stack");
    parse_args->Add_Integer_Argument("-num_bodies",6,"number of total bodies");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    output_directory=STRING_UTILITIES::string_sprintf("MPI_Example/Test_%d",test_number);
    processes_per_dimension.x=parse_args->Get_Integer_Value("-xprocs");
    processes_per_dimension.y=parse_args->Get_Integer_Value("-yprocs");
    processes_per_dimension.z=parse_args->Get_Integer_Value("-zprocs");
    width=parse_args->Get_Integer_Value("-width");
    height=parse_args->Get_Integer_Value("-height");
    num_bodies=parse_args->Get_Integer_Value("-num_bodies");
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
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
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}
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
            rigid_body->X()=TV((T)current_x, (height+1-i)*5+80,0);
            rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
            current_x += 5;}
        first_x -= 2.5;}

    for(int i=0;i<2;i++){
        RIGID_BODY<TV>* rigid_body1=&tests.Add_Rigid_Body("subdivided_box",(T)10,(T).1);
        rigid_body1->X()=TV(first_x-10,i*20,0);
        rigid_body1->is_static = true;
        rigid_body1->Set_Name("left_box");

        RIGID_BODY<TV>* rigid_body2=&tests.Add_Rigid_Body("subdivided_box",(T)10,(T).1);
        rigid_body2->X()=TV(current_x+7.5,i*20,0);
        rigid_body2->is_static = true;
        rigid_body2->Set_Name("right_box");}

    tests.Add_Ground(1, -10);
    last_frame = 400;
    rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(rigid_body_collection, true));
}
//#####################################################################
// Function Stacked_Boxes
//#####################################################################
void Stacked_Boxes() {
    int height = 50;

    for (int i = 1; i <= height; i++) {
        for (int j =  1; j < 20 ; j++) {
            RIGID_BODY<TV>* rigid_body = &tests.Add_Rigid_Body("subdivided_box",(T)2, (T).1);
            rigid_body->X()=TV(5*j - 50, 5*i+10,0);
            rigid_body->Set_Coefficient_Of_Restitution((T)0.5);}}

    for (int i = 0; i < height/4; i++) {
        RIGID_BODY<TV>* rigid_body1=&tests.Add_Rigid_Body("subdivided_box",(T)10,(T).1);
        rigid_body1->X()=TV(-58,i*20,0);
        rigid_body1->is_static = true;
        rigid_body1->Set_Name("left_box");
    
        RIGID_BODY<TV>* rigid_body2=&tests.Add_Rigid_Body("subdivided_box",(T)10,(T).1);
        rigid_body2->X()=TV(58,i*20,0);
        rigid_body2->is_static = true;
        rigid_body2->Set_Name("right_box");}

    tests.Add_Ground(1, -10);
    last_frame = 250;
    rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(rigid_body_collection, true));
}
//#####################################################################
// Function Partition_Test
//#####################################################################
void Partition_Test() {
    RIGID_BODY<TV>* rigid_body = &tests.Add_Rigid_Body("box", (T)1, (T).1);
    rigid_body->X() = TV(0,0,0);
    rigid_body->Twist().linear = TV(1,0,0);

    RIGID_BODY<TV>* rigid_body2 = &tests.Add_Rigid_Body("box", (T)1, (T).1);
    rigid_body2->X() = TV(10,3,0);
    rigid_body2->Twist().linear = TV(-1,0,0);

    last_frame = 200;
}
//#####################################################################
// Function Contact_Test_1
//#####################################################################
void Contact_Test_1() {
    for (int i=0;i<width;i++) {
        for (int j=0;j<height;j++) {
            RIGID_BODY<TV>* rigid_body = &tests.Add_Rigid_Body("subdivided_box",(T)1,(T).1);
            rigid_body->X() = TV(i*4,j*2,0);}}

    tests.Add_Ground(1, -1);
    rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(rigid_body_collection, true));
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
            rigid_body->X() = TV((T)j*2+i*0.5,i*2,0);}
        width--;}

    tests.Add_Ground(1, -1);
    rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(rigid_body_collection, true));
}
//#####################################################################
// Function Simple_Collision_Test
//#####################################################################
void Simple_Collision_Test() {
    RIGID_BODY<TV>* rigid_body = &tests.Add_Rigid_Body("box", (T)1, (T).1);
    rigid_body->X() = TV(0,0,0);
    rigid_body->Twist().linear = TV(2,0,0);

    RIGID_BODY<TV>* rigid_body2 = &tests.Add_Rigid_Body("box", (T)1, (T).1);
    rigid_body2->X() = TV(10,3,0);
    rigid_body2->Twist().linear = TV(-1,0,0);

    RIGID_BODY<TV>* rigid_body3 = &tests.Add_Rigid_Body("box", (T)1, (T).1);
    rigid_body3->X() = TV(10,0,0);
    rigid_body3->Twist().linear = TV(-1,0,0);

    last_frame = 200;
}
//#####################################################################
// Function Collision_Test
//#####################################################################
void Collision_Test() {
    for (int i=0;i<height;i++) {
        RIGID_BODY<TV>* left_body = &tests.Add_Rigid_Body("subdivided_box",(T)1,(T).1);
        RIGID_BODY<TV>* right_body = &tests.Add_Rigid_Body("subdivided_box",(T)1,(T).1);
        left_body->X() = TV(0,i*3,0);
        right_body->X() = TV((T)(width+1)*2.15,i*3,0);
        left_body->Set_Coefficient_Of_Restitution((T)1.0);
        right_body->Set_Coefficient_Of_Restitution((T)1.0);
        left_body->is_static=true;
        right_body->is_static=true;

        for (int j=1;j<=width;j++) {
            RIGID_BODY<TV>* rigid_body = &tests.Add_Rigid_Body("subdivided_box",(T)1,(T).1);
            rigid_body->X() = TV((T)j*2.15,i*3,0);
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
        rigid_body->X() = TV((T)1.75*i,(T)1.75*i,0);}
}
//#####################################################################
// Function Analytic_Contact
//#####################################################################
void Many_Sphere_Test()
{
    last_frame=300;
    typedef typename RIGID_BODY_POLICY<TV>::INERTIA_TENSOR T_INERTIA_TENSOR;    
     VECTOR<T,3> num_bodies=TV(40,200,40);
     for(int i=1;i<=num_bodies.x;i++) for(int j=1;j<=num_bodies.y;j++) for(int k=1;k<=num_bodies.z;k++){
        RIGID_BODY_PARTICLES<TV>& particles=rigid_body_collection.rigid_body_particle;    
        RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection);
        T radius=1;TV center;
        T offset=.1,offset2=.1;if(j%2==0){offset=-.1;offset2=-.1;}
        center=VECTOR<T,3>(3*i+offset,3*j+1,3*k+offset2);
        particles.X(rigid_body.particle_index)=center;
        particles.rotation(rigid_body.particle_index)=ROTATION<TV>();
        particles.twist(rigid_body.particle_index)=TWIST<TV>();
        particles.mass(rigid_body.particle_index)=1;
        particles.inertia_tensor(rigid_body.particle_index)=T_INERTIA_TENSOR();particles.inertia_tensor(rigid_body.particle_index)+=1;
        SPHERE<TV> sphere(TV(0,0,0),radius);
        rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >(sphere));
        rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
        rigid_body_collection.rigid_geometry_collection.collision_body_list->Add_Body(new RIGID_COLLISION_GEOMETRY<TV>(rigid_body),rigid_body.particle_index,true);
        rigid_body_collection.rigid_geometry_collection.collision_body_list->Get_Collision_Geometry(rigid_body.particle_index)->add_to_spatial_partition=true;}
    tests.Add_Ground((T).5,0,0);
    RIGID_GRAVITY<TV> *gravity=new RIGID_GRAVITY<TV>(rigid_body_collection,true);
    rigid_body_collection.Add_Force(gravity);
}
//#####################################################################
// Function Cluster_Fracture
//#####################################################################
void Break_Levelset()
{
    typedef VECTOR<int,TV::dimension> TV_INT;
    typedef typename RIGID_BODY_POLICY<TV>::INERTIA_TENSOR T_INERTIA_TENSOR;
    ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> children;
    ARRAY<int>* referenced_rigid_particles=new ARRAY<int>;
    RIGID_BODY_PARTICLES<TV>& particles=rigid_body_collection.rigid_body_particle;    
    GRID<TV>& grid=*new GRID<TV>;
    ARRAY<T,VECTOR<int,3> >& phi=*new ARRAY<T,VECTOR<int,3> >;
    LEVELSET_3D<GRID<TV> > levelset(grid,phi);
    data_directory="../../../Public_Data/";
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Rigid_Bodies/sphere_66k.phi",data_directory.c_str()),levelset);
    TV_INT counts=TV_INT::All_Ones_Vector()*25;
    GRID<TV> body_grid(counts,levelset.grid.Domain());
    TRIANGULATED_SURFACE<T>* surface=0;    
    for(typename GRID<TV>::CELL_ITERATOR iterator(body_grid);iterator.Valid();iterator.Next()){
        TV_INT cell=levelset.grid.Cell(iterator.Location(),0);
        if(!phi.Valid_Index(cell)) continue;
        if(phi(cell)>0) continue;
        RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection);
        children.Append(rigid_body.particle_index);
        referenced_rigid_particles->Append(rigid_body.particle_index);
        particles.X(rigid_body.particle_index)=iterator.Location();
        particles.rotation(rigid_body.particle_index)=ROTATION<TV>();
        particles.twist(rigid_body.particle_index)=TWIST<TV>();
        particles.mass(rigid_body.particle_index)=1;
        particles.inertia_tensor(rigid_body.particle_index)=T_INERTIA_TENSOR();particles.inertia_tensor(rigid_body.particle_index)+=1;
        BOX<TV> box(body_grid.dX*-.5,body_grid.dX*.5);
        //SPHERE<TV> box(TV(),body_grid.min_dX*.5);
        rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<BOX<TV> >(box));
        //rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >(box));
        if(!surface) surface=TESSELLATION::Generate_Triangles(box);
        rigid_body.Add_Structure(*surface);  
        rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
        rigid_body_collection.rigid_geometry_collection.collision_body_list->Add_Body(new RIGID_COLLISION_GEOMETRY<TV>(rigid_body),rigid_body.particle_index,true);
        rigid_body_collection.rigid_geometry_collection.collision_body_list->Get_Collision_Geometry(rigid_body.particle_index)->add_to_spatial_partition=true;}
    tests.Add_Ground((T).5,-2,0);
    
    rigid_body_collection.Add_Force(new RIGID_GRAVITY<TV>(rigid_body_collection,referenced_rigid_particles));
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    rigid_body_collection.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Positions();
}
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2004, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CLUSTER_EXAMPLE
//##################################################################### 
#ifndef __CLUSTER_EXAMPLE__
#define __CLUSTER_EXAMPLE__

#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_IMPULSE_ACCUMULATOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/RIGID_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Bindings/RIGID_BODY_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/FRACTURE_EVOLUTION_CLUSTER_3D.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include "../ARB_PARAMETERS.h"
#include <Rigid_Bodies/CONSTRAINED_POINT_IN_RIGID_BODY.h>
#include <Rigid_Bodies/RIGID_BODY_CLUSTER_3D.h>
namespace PhysBAM{

template<class T,class RW>
class CLUSTER_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM< GRID<TV>,RW>
{
public:
    typedef VECTOR<T,3> TV;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM< GRID<TV>,RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM< GRID<TV>,RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM< GRID<TV>,RW>::frame_rate;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM< GRID<TV>,RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM< GRID<TV>,RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM< GRID<TV>,RW>::output_directory;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM< GRID<TV>,RW>::solids_parameters;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM< GRID<TV>,RW>::write_last_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM< GRID<TV>,RW>::data_directory;
    
    ARTICULATED_RIGID_BODY<TV>* arb;
    RIGID_BODY<TV> *shelf11,*shelf12,*shelf21,*shelf22;
    int current_frame,start_move,end_move;
    T increment;
    int selection;
    int cluster1,cluster2;
    bool add_ground;
    PARAMETER_LIST parameter_list;
    RIGID_BODY<TV> *speed_up_test1,*speed_up_test2;
    bool springs;

    CLUSTER_EXAMPLE()
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>(0,FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >::NONE)
    {
        //solids_parameters.perform_self_collision=false;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.cfl=(T).1;
        //solids_parameters.perform_collision_body_collisions=false;
        //restart=true;
        restart_frame=30;
        last_frame=2000;
        frame_rate=24;
        output_directory="Clusters/output";
        std::cout << "Frame rate: "<<frame_rate<<std::endl;

        arb=new ARTICULATED_RIGID_BODY<TV>(this->solids_parameters.rigid_body_parameters.list);
        this->solids_parameters.rigid_body_parameters.Set_Articulated_Rigid_Body(arb); 
        current_frame=0;
        increment=(T).05;
        start_move=5;end_move=40;
        shelf11=shelf12=shelf21=shelf22=0;
        arb->Set_Iterative_Tolerance((T)1e-3);
        arb->Set_Extra_Iterations_Per_Contact_Level_Factor(100);
        arb->Set_Extra_Iterations_Per_Shock_Propagation_Level_Factor(100);
        arb->Set_Poststabilization_Iterations(100);
        arb->Set_Use_Shock_Propagation(false);
        arb->Set_Do_Final_Pass(false);
        arb->Use_PD_Actuators();
        write_last_frame=true;
        selection=6;
        springs=false;

//        ARB_PARAMETERS::Read_Common_Parameters("Blocks/example.param",*this,parameter_list);
    }

    ~CLUSTER_EXAMPLE()
    {
        delete arb;
    }

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    add_ground=true;
    int id=0;
    RIGID_BODY<TV> *rigid_body=0;
    int num_joints=0,num_rigid_bodies=0;
    switch(selection){
      case 1:
        Rigid_Constraint_Speed_Up_Test(num_joints,num_rigid_bodies,VECTOR<T,3>(0,0,0),false);
        break;
      case 2:
        Two_Large_Cubes(num_joints,num_rigid_bodies,VECTOR<T,3>(0,7,0));
        break;
      case 3:
        Cubes_With_Constraints(num_joints,num_rigid_bodies,VECTOR<T,3>(0,4,0));
        break;
      case 4:
        Simple_Spring();
        break;
      case 5:
        Weird_Cluster(num_joints,num_rigid_bodies,VECTOR<T,3>(0,10,0),false);
        break;
      case 6:
        One_Large_Cube(num_joints,num_rigid_bodies,VECTOR<T,3>());
        break;
      case 7:
        Cluster_Collisions(num_joints,num_rigid_bodies,VECTOR<T,3>());
        break;
      case 8:
        Simple_Test();
        break;
      default:
        break;
    }

    if(add_ground){
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"ground");
        rigid_body=arb->rigid_bodies_list(id);
        rigid_body->frame.t=VECTOR<T,3>(0,-5,0);
        rigid_body->twist.linear=VECTOR<T,3>(0,0,0);
        rigid_body->Set_Coefficient_Of_Restitution(1);
        rigid_body->Set_Coefficient_Of_Friction(.5);
        rigid_body->Set_Name("ground");
        rigid_body->is_static=true;
        rigid_body->add_to_spatial_partition=false;}

    RIGID_BODY_LIST<VECTOR<T,3> >& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    for(int i=1;i<=rigid_body_list.Number_Of_Active_Elements();i++) if(!rigid_body_list.rigid_bodies(i)->is_static){
        rigid_body_list.rigid_bodies(i)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
    }
    //solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM< GRID<TV>,RW>::Initialize_Bodies();
    std::cout<<"done initializing example\n";
}


//#####################################################################
// Rigid Constriant With 2 Blocks
//#####################################################################
void Rigid_Constraint_Speed_Up_Test(int& num_joints,int& num_rigid_bodies,VECTOR<T,3> shift,const bool velocity_on=true)
{
    // CLUSTER 1
    int id3,id4,id5,id6,id7;
    id3=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"subdivided_box");
    speed_up_test1=arb->rigid_bodies_list.rigid_bodies(id3);
    speed_up_test1->frame.t=VECTOR<T,3>(0,2,0)+shift;
    speed_up_test1->Set_Coefficient_Of_Restitution(0.5);
    speed_up_test1->Set_Coefficient_Of_Friction(0.5);
    speed_up_test1->Set_Name("parent");

    id4=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"subdivided_box");
    speed_up_test2=arb->rigid_bodies_list.rigid_bodies(id4);
    speed_up_test2->frame.t=VECTOR<T,3>(2,2,0)+shift;
    speed_up_test2->twist.linear=VECTOR<T,3>(0,0,0);
    speed_up_test2->Set_Coefficient_Of_Restitution(0.5);
    speed_up_test2->Set_Coefficient_Of_Friction(.5);
    speed_up_test2->Set_Name("child");

    // CLUSTER 2
    id5=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"subdivided_box");
    speed_up_test1=arb->rigid_bodies_list.rigid_bodies(id5);
    speed_up_test1->frame.t=VECTOR<T,3>(0,4,0)+shift;
    speed_up_test1->Set_Coefficient_Of_Restitution(0.5);
    speed_up_test1->Set_Coefficient_Of_Friction(0.5);
    speed_up_test1->Set_Name("parent2");

    id6=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"subdivided_box");
    speed_up_test2=arb->rigid_bodies_list.rigid_bodies(id6);
    speed_up_test2->frame.t=VECTOR<T,3>(0,6,0)+shift;
    speed_up_test2->Set_Coefficient_Of_Restitution(0.5);
    speed_up_test2->Set_Coefficient_Of_Friction(.5);
    speed_up_test2->Set_Name("child2");

    id7=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"subdivided_box");
    speed_up_test2=arb->rigid_bodies_list.rigid_bodies(id7);
    speed_up_test2->frame.t=VECTOR<T,3>(2,8,0)+shift;
    speed_up_test2->Set_Coefficient_Of_Restitution(0.5);
    speed_up_test2->Set_Coefficient_Of_Friction(.5);
    speed_up_test2->Set_Name("child2.5");

    ARRAY<RIGID_BODY<TV>*>* bodies=new ARRAY<RIGID_BODY<TV>*>();
    bodies->Append(solids_parameters.rigid_body_parameters.list(id3));
    bodies->Append(solids_parameters.rigid_body_parameters.list(id4));
    cluster1=solids_parameters.rigid_body_parameters.list.Add_Cluster_Body(bodies,solids_parameters.collision_body_list);
    arb->rigid_bodies_list(cluster1)->Set_Name("combo_square");

    ARRAY<RIGID_BODY<TV>*>* bodies2=new ARRAY<RIGID_BODY<TV>*>();
    bodies2->Append(solids_parameters.rigid_body_parameters.list(id5));
    bodies2->Append(solids_parameters.rigid_body_parameters.list(id6));
    bodies2->Append(solids_parameters.rigid_body_parameters.list(id7));
    cluster2=solids_parameters.rigid_body_parameters.list.Add_Cluster_Body(bodies2,solids_parameters.collision_body_list);
    arb->rigid_bodies_list(cluster2)->Set_Name("combo_square2");

}
//#####################################################################
// Rigid Constriant With 2 Blocks
//#####################################################################
void Cubes_With_Constraints(int& num_joints,int& num_rigid_bodies,VECTOR<T,3> shift,const bool velocity_on=true)
{
    DEFORMABLE_OBJECT<TV>& deformable_object=solid_body_collection.deformable_object;
    PARTICLES<TV>& particles=deformable_object.particles;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=deformable_object.rigid_body_particles;
    ARRAY<RIGID_BODY<TV>*>& rigid_bodies=solids_parameters.rigid_body_parameters.list.rigid_bodies;
    SEGMENTED_CURVE<TV>& segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);

    particles.Store_Velocity();
    // CLUSTER 1
    RIGID_BODY<TV>* rigid_body;
    int ids[9];

    int block=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"subdivided_box",2);
    rigid_body=arb->rigid_bodies_list.rigid_bodies(block);
    rigid_body->frame.t=VECTOR<T,3>(16,0,0)+shift;
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(.5);
    rigid_body->Set_Name("block");

    int block_root=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"subdivided_box");
    rigid_body=arb->rigid_bodies_list.rigid_bodies(block_root);
    rigid_body->frame.t=VECTOR<T,3>(13,3,1)+shift;
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(.5);
    rigid_body->Set_Name("block_root");
    rigid_body->is_static=true;


    for(int i=0;i<9;i++){
        ids[i]=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"subdivided_box");
        rigid_body=arb->rigid_bodies_list.rigid_bodies(ids[i]);
        rigid_body->Set_Coefficient_Of_Restitution(0.5);
        rigid_body->Set_Coefficient_Of_Friction(0.5);
        rigid_body->Set_Name("child");
        rigid_body->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

    }
    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list); // we add all the bodies now because then they're not double added
    solids_parameters.rigid_body_parameters.list.Add_Rigid_Body_Particles(rigid_body_particles);

    arb->rigid_bodies_list.rigid_bodies(ids[0])->frame.t=VECTOR<T,3>(-1,-1,-1)+shift;
    arb->rigid_bodies_list.rigid_bodies(ids[1])->frame.t=VECTOR<T,3>(1,-1,-1)+shift;
    arb->rigid_bodies_list.rigid_bodies(ids[2])->frame.t=VECTOR<T,3>(-1,1,-1)+shift;
    arb->rigid_bodies_list.rigid_bodies(ids[3])->frame.t=VECTOR<T,3>(1,1,-1)+shift;
    arb->rigid_bodies_list.rigid_bodies(ids[4])->frame.t=VECTOR<T,3>(-1,-1,1)+shift;
    arb->rigid_bodies_list.rigid_bodies(ids[5])->frame.t=VECTOR<T,3>(1,-1,1)+shift;
    arb->rigid_bodies_list.rigid_bodies(ids[6])->frame.t=VECTOR<T,3>(-1,1,1)+shift;
    arb->rigid_bodies_list.rigid_bodies(ids[7])->frame.t=VECTOR<T,3>(1,1,1)+shift;

    arb->rigid_bodies_list.rigid_bodies(ids[8])->frame.t=VECTOR<T,3>(-3,3,1)+shift;
    arb->rigid_bodies_list.rigid_bodies(ids[8])->is_static=true;

    if(springs){
        // add particles and set up segmented curve
        for(int i=0;i<4;i++){
            particles.array_collection->Add_Element();
            particles.mass(i)=1;}

        // probably don't need these -- can just sync wiht bindings, but doesn't hurt
        particles.X(1)=rigid_bodies(block)->frame*VECTOR<T,3>(-1,1,1);
        particles.X(2)=rigid_bodies(block_root)->frame*VECTOR<T,3>();
        particles.X(3)=rigid_bodies(ids[6])->frame*VECTOR<T,3>();
        particles.X(4)=rigid_bodies(ids[8])->frame*VECTOR<T,3>();
        segmented_curve.mesh.elements.Append(1,2);
        segmented_curve.mesh.elements.Append(3,4);
        solid_body_collection.deformable_object.Add_Structure(&segmented_curve);
        
        solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new RIGID_BODY_BINDING<TV>(particles,1,rigid_body_particles,rigid_body_particles.id.array.Find(block),VECTOR<T,3>(-1,1,1)));
        solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new RIGID_BODY_BINDING<TV>(particles,2,rigid_body_particles,rigid_body_particles.id.array.Find(block_root),VECTOR<T,3>(0,0,0)));
        solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new RIGID_BODY_BINDING<TV>(particles,3,rigid_body_particles,rigid_body_particles.id.array.Find(ids[6]),VECTOR<T,3>(0,0,0)));
        solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new RIGID_BODY_BINDING<TV>(particles,4,rigid_body_particles,rigid_body_particles.id.array.Find(ids[8]),VECTOR<T,3>(0,0,0)));
        segmented_curve.Update_Number_Nodes();
        deformable_object.particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
        
        //solid_body_collection.Add_Force(new GRAVITY<TV>(particles));
        solid_body_collection.Add_Force(Create_Edge_Springs(segmented_curve,(T)1e2,(T)1));
        solid_body_collection.Update_Fragments();
        std::cout<<"fragments: "<<deformable_object.rigid_body_particles_of_fragment<<std::endl;
        
        ARRAY<RIGID_BODY<TV>*>* bodies=new ARRAY<RIGID_BODY<TV>*>();
        for(int i=0;i<8;i++){
            bodies->Append(solids_parameters.rigid_body_parameters.list(ids[i]));}
        cluster1=solids_parameters.rigid_body_parameters.list.Add_Cluster_Body(bodies,solids_parameters.collision_body_list);
        arb->rigid_bodies_list(cluster1)->Set_Name("combo_square");
        
        ((RIGID_BODY_CLUSTER_3D<T>*)solids_parameters.rigid_body_parameters.list(cluster1))->Fix_Bindings_For_Cluster(solid_body_collection.deformable_body_collection.binding_list,&rigid_body_particles,&solid_body_collection.deformable_object);
        solid_body_collection.Update_Fragments();
        std::cout<<"fragments: "<<deformable_object.rigid_body_particles_of_fragment<<std::endl;
    }
    else{
        // add in a joint constraint
        JOINT<TV>* joint;
        if(!this->restart){
            joint=new POINT_JOINT<TV>();arb->joint_mesh.Add_Joint(joint);
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(VECTOR<T,3>(1,-1,1)));    
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(VECTOR<T,3>(-2,2,2)));
            arb->joint_mesh.Add_Articulation(block_root,block,1);} 
        if(!this->restart){
            joint=new POINT_JOINT<TV>();arb->joint_mesh.Add_Joint(joint);
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(VECTOR<T,3>(1,-1,1)));    
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(VECTOR<T,3>(-1,1,1)));
            arb->joint_mesh.Add_Articulation(ids[8],ids[6],2);} 
        if(!this->restart){
            joint=new POINT_JOINT<TV>();arb->joint_mesh.Add_Joint(joint);
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(VECTOR<T,3>(1,1,1)));    
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(VECTOR<T,3>(-1,1,1)));
            arb->joint_mesh.Add_Articulation(ids[6],ids[7],3);}
        
        ARRAY<RIGID_BODY<TV>*>* bodies=new ARRAY<RIGID_BODY<TV>*>();
        for(int i=0;i<8;i++){
            bodies->Append(solids_parameters.rigid_body_parameters.list(ids[i]));
            std::cout<<"just added body "<<solids_parameters.rigid_body_parameters.list(ids[i])->name<<" to cluster"<<std::endl;}
        cluster1=solids_parameters.rigid_body_parameters.list.Add_Cluster_Body(bodies,solids_parameters.collision_body_list);
        arb->rigid_bodies_list(cluster1)->Set_Name("combo_square");
        ((RIGID_BODY_CLUSTER_3D<T>*)solids_parameters.rigid_body_parameters.list(cluster1))->Fix_Joint_Constraints_For_Cluster(arb);
    }
}
//#####################################################################
// Rigid Constriant With 2 Blocks
//#####################################################################
void One_Large_Cube(int& num_joints,int& num_rigid_bodies,VECTOR<T,3> shift,const bool velocity_on=true)
{
    solids_parameters.fracture=true;
    FRACTURE_EVOLUTION_CLUSTER_3D<TV>* fracture_evolution=new FRACTURE_EVOLUTION_CLUSTER_3D<TV>(solids_parameters);
    solids_parameters.Set_Fracture_Evolution(fracture_evolution);

    // CLUSTER 1
    RIGID_BODY<TV>* rigid_body;
    int ids[8];
    for(int i=0;i<8;i++){
        ids[i]=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"subdivided_box");
        rigid_body=arb->rigid_bodies_list.rigid_bodies(ids[i]);
        rigid_body->Set_Coefficient_Of_Restitution(0.5);
        rigid_body->Set_Coefficient_Of_Friction(0.5);
        rigid_body->Set_Name(STRING_UTILITIES::string_sprintf("child::%d",ids[i]));
        // rigid_body->frame=FRAME<TV>(QUATERNION<T>::From_Euler_Angles(79.9655*pi/180,12.0032*pi/180,50.4023*pi/180))*rigid_body->frame;
        rigid_body->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
    }
    arb->rigid_bodies_list.rigid_bodies(ids[0])->frame.t=VECTOR<T,3>(-1,-1,-1);
    arb->rigid_bodies_list.rigid_bodies(ids[1])->frame.t=VECTOR<T,3>(1,-1,-1);
    arb->rigid_bodies_list.rigid_bodies(ids[2])->frame.t=VECTOR<T,3>(-1,1,-1);
    arb->rigid_bodies_list.rigid_bodies(ids[3])->frame.t=VECTOR<T,3>(1,1,-1);
    arb->rigid_bodies_list.rigid_bodies(ids[4])->frame.t=VECTOR<T,3>(-1,-1,1);
    arb->rigid_bodies_list.rigid_bodies(ids[5])->frame.t=VECTOR<T,3>(1,-1,1);
    arb->rigid_bodies_list.rigid_bodies(ids[6])->frame.t=VECTOR<T,3>(-1,1,1);
    arb->rigid_bodies_list.rigid_bodies(ids[7])->frame.t=VECTOR<T,3>(1,1,1);

    FRAME<TV> shift_frame(shift);
    for(int i=0;i<8;i++){
        rigid_body=arb->rigid_bodies_list.rigid_bodies(ids[i]);
        rigid_body->frame=shift_frame*FRAME<TV>(QUATERNION<T>::From_Euler_Angles(29.9655*pi/180,32.0032*pi/180,50.4023*pi/180))*rigid_body->frame;}


    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    ARRAY<RIGID_BODY<TV>*>* bodies=new ARRAY<RIGID_BODY<TV>*>();
    for(int i=0;i<8;i++) bodies->Append(solids_parameters.rigid_body_parameters.list(ids[i]));
    cluster1=solids_parameters.rigid_body_parameters.list.Add_Cluster_Body(bodies,solids_parameters.collision_body_list);
    arb->rigid_bodies_list(cluster1)->Set_Name("combo_square");
    arb->rigid_bodies_list(cluster1)->impulse_accumulator=new RIGID_BODY_IMPULSE_ACCUMULATOR<TV,3>(*rigid_body);
    ((RIGID_BODY_CLUSTER_3D<T>*)arb->rigid_bodies_list(cluster1))->perform_cluster_breaks=true;
    ((RIGID_BODY_CLUSTER_3D<T>*)arb->rigid_bodies_list(cluster1))->Initialize_Initial_Strain();
    ARRAY<ARRAY<T> > strain;strain.Resize(((RIGID_BODY_CLUSTER_3D<T>*)arb->rigid_bodies_list(cluster1))->constituent_bodies.m);
//    for(int i=0;i<strain.m;i++){strain(i).Resize(((RIGID_BODY_CLUSTER_3D<T>*)arb->rigid_bodies_list(cluster1))->constituent_bodies.m);strain(i).Fill((T)1.01);}
    for(int i=0;i<strain.m;i++){strain(i).Resize(((RIGID_BODY_CLUSTER_3D<T>*)arb->rigid_bodies_list(cluster1))->constituent_bodies.m);strain(i).Fill((T)0);}
    for(int i=0;i<4;i++) for(int j=0;j<4;j++) strain(i)(j)=1.5;
    ((RIGID_BODY_CLUSTER_3D<T>*)arb->rigid_bodies_list(cluster1))->Initialize_Allowable_Strain(&strain);
}
//#####################################################################
// Rigid Constriant With 2 Blocks
//#####################################################################
void Two_Large_Cubes(int& num_joints,int& num_rigid_bodies,VECTOR<T,3> shift,const bool velocity_on=true)
{
    // CLUSTER 1
    RIGID_BODY<TV>* rigid_body;
    int ids[8];
    for(int i=0;i<8;i++){
        ids[i]=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"subdivided_box");
        rigid_body=arb->rigid_bodies_list.rigid_bodies(ids[i]);
        rigid_body->Set_Coefficient_Of_Restitution(0.5);
        rigid_body->Set_Coefficient_Of_Friction(0.5);
        rigid_body->Set_Name(STRING_UTILITIES::string_sprintf("child::%d",ids[i]));
        // rigid_body->frame=FRAME<TV>(QUATERNION<T>::From_Euler_Angles(79.9655*pi/180,12.0032*pi/180,50.4023*pi/180))*rigid_body->frame;
        rigid_body->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
    }
    arb->rigid_bodies_list.rigid_bodies(ids[0])->frame.t=VECTOR<T,3>(-1,-1,-1);
    arb->rigid_bodies_list.rigid_bodies(ids[1])->frame.t=VECTOR<T,3>(1,-1,-1);
    arb->rigid_bodies_list.rigid_bodies(ids[2])->frame.t=VECTOR<T,3>(-1,1,-1);
    arb->rigid_bodies_list.rigid_bodies(ids[3])->frame.t=VECTOR<T,3>(1,1,-1);
    arb->rigid_bodies_list.rigid_bodies(ids[4])->frame.t=VECTOR<T,3>(-1,-1,1);
    arb->rigid_bodies_list.rigid_bodies(ids[5])->frame.t=VECTOR<T,3>(1,-1,1);
    arb->rigid_bodies_list.rigid_bodies(ids[6])->frame.t=VECTOR<T,3>(-1,1,1);
    arb->rigid_bodies_list.rigid_bodies(ids[7])->frame.t=VECTOR<T,3>(1,1,1);

    FRAME<TV> shift_frame(shift);
    for(int i=0;i<8;i++){
        rigid_body=arb->rigid_bodies_list.rigid_bodies(ids[i]);
        rigid_body->frame=shift_frame*FRAME<TV>(QUATERNION<T>::From_Euler_Angles(79.9655*pi/180,12.0032*pi/180,50.4023*pi/180))*rigid_body->frame;}

//    arb->rigid_bodies_list.rigid_bodies(ids[4])->Set_Name("heavy");

    int block=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"subdivided_box",2);
    rigid_body=arb->rigid_bodies_list.rigid_bodies(block);
    rigid_body->frame.t=VECTOR<T,3>(6,0,0)+shift;
    rigid_body->frame.r=QUATERNION<T>::From_Euler_Angles(79.9655*pi/180,12.0032*pi/180,50.4023*pi/180);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(.5);
    rigid_body->Set_Name("block");

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    ARRAY<RIGID_BODY<TV>*>* bodies=new ARRAY<RIGID_BODY<TV>*>();
    for(int i=0;i<8;i++) bodies->Append(solids_parameters.rigid_body_parameters.list(ids[i]));
    cluster1=solids_parameters.rigid_body_parameters.list.Add_Cluster_Body(bodies,solids_parameters.collision_body_list);
    arb->rigid_bodies_list(cluster1)->Set_Name("combo_square");

//    arb->rigid_bodies_list(cluster1)->frame.r=QUATERNION<T>(.5,VECTOR<T,3>(0,0,1));
//    arb->rigid_bodies_list(block)->frame.r=QUATERNION<T>(.5,VECTOR<T,3>(0,0,1));
    //arb->rigid_bodies_list(cluster1)->frame.r=QUATERNION<T>::From_Euler_Angles(79.9655*pi/180,12.0032*pi/180,50.4023*pi/180);
    //arb->rigid_bodies_list(block)->frame.r=QUATERNION<T>::From_Euler_Angles(79.9655*pi/180,12.0032*pi/180,50.4023*pi/180);

}
//#####################################################################
// Rigid Constriant With 2 Blocks
//#####################################################################
void Cluster_Collisions(int& num_joints,int& num_rigid_bodies,VECTOR<T,3> shift,const bool velocity_on=true)
{
    // CLUSTER 1
    RIGID_BODY<TV>* rigid_body;
    int ids[16];
    for(int i=0;i<16;i++){
        ids[i]=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"box");
        rigid_body=arb->rigid_bodies_list.rigid_bodies(ids[i]);
        rigid_body->Set_Coefficient_Of_Restitution(0.5);
        rigid_body->Set_Coefficient_Of_Friction(0.5);
        rigid_body->Set_Name(STRING_UTILITIES::string_sprintf("block::%d",ids[i]));
        // rigid_body->frame=FRAME<TV>(QUATERNION<T>::From_Euler_Angles(79.9655*pi/180,12.0032*pi/180,50.4023*pi/180))*rigid_body->frame;
              rigid_body->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
    }
    for(int i=0;i<4;i++) for(int j=0;j<4;j++) arb->rigid_bodies_list.rigid_bodies(ids[i*4+j])->frame.t=VECTOR<T,3>(i*3,-1,j*3);
//    arb->rigid_bodies_list.rigid_bodies(ids[0])->frame.t=VECTOR<T,3>(0,-1,0);
//    arb->rigid_bodies_list.rigid_bodies(ids[1])->frame.t=VECTOR<T,3>(3,-1,0);

    int ids_ball[16];
    for(int i=0;i<16;i++){
        ids_ball[i]=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"sphere");
        rigid_body=arb->rigid_bodies_list.rigid_bodies(ids_ball[i]);
        rigid_body->Set_Coefficient_Of_Restitution(0.5);
        rigid_body->Set_Coefficient_Of_Friction(0.5);
        rigid_body->Set_Name(STRING_UTILITIES::string_sprintf("sphere::%d",ids_ball[i]));
    }
    for(int i=0;i<4;i++) for(int j=0;j<4;j++) arb->rigid_bodies_list.rigid_bodies(ids_ball[i*4+j])->frame.t=VECTOR<T,3>(i*3,2,j*3);
    //arb->rigid_bodies_list.rigid_bodies(ids_ball[0])->frame.t=VECTOR<T,3>(0,2,0);
    //arb->rigid_bodies_list.rigid_bodies(ids_ball[1])->frame.t=VECTOR<T,3>(3,2,0);

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    ARRAY<RIGID_BODY<TV>*>* bodies=new ARRAY<RIGID_BODY<TV>*>();
    for(int i=0;i<16;i++) bodies->Append(solids_parameters.rigid_body_parameters.list(ids[i]));
    cluster1=solids_parameters.rigid_body_parameters.list.Add_Cluster_Body(bodies,solids_parameters.collision_body_list);
    arb->rigid_bodies_list(cluster1)->Set_Name("combo_squares");
    ((RIGID_BODY_CLUSTER_3D<T>*)arb->rigid_bodies_list(cluster1))->perform_cluster_breaks=false;

}
//#####################################################################
// Simple Test
//#####################################################################
void Simple_Test()
{
    // CLUSTER 1
    RIGID_BODY<TV>* rigid_body;
    int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"box");
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(0.5);
    rigid_body->Set_Name(STRING_UTILITIES::string_sprintf("block::%d",id));
    // rigid_body->frame=FRAME<TV>(QUATERNION<T>::From_Euler_Angles(79.9655*pi/180,12.0032*pi/180,50.4023*pi/180))*rigid_body->frame;
    rigid_body->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
        
    int ball=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"sphere");
    rigid_body=arb->rigid_bodies_list.rigid_bodies(ball);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(0.5);
    rigid_body->frame.t.y=3;
    rigid_body->Set_Name(STRING_UTILITIES::string_sprintf("sphere::%d",ball));
    
    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
    
    if(true){
        ARRAY<RIGID_BODY<TV>*>* bodies=new ARRAY<RIGID_BODY<TV>*>();
        bodies->Append(solids_parameters.rigid_body_parameters.list(id));
        cluster1=solids_parameters.rigid_body_parameters.list.Add_Cluster_Body(bodies,solids_parameters.collision_body_list);
        arb->rigid_bodies_list(cluster1)->Set_Name("combo_squares");
        ((RIGID_BODY_CLUSTER_3D<T>*)arb->rigid_bodies_list(cluster1))->perform_cluster_breaks=false;}

}
//#####################################################################
// Rigid Constriant With 2 Blocks
//#####################################################################
void Weird_Cluster(int& num_joints,int& num_rigid_bodies,VECTOR<T,3> shift,const bool velocity_on=true)
{
    // CLUSTER 1
    RIGID_BODY<TV>* rigid_body;
    int id1,id2;

    id1=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"subdivided_box");
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id1);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(1);
    rigid_body->Set_Name("child1");
    rigid_body->frame.t=VECTOR<T,3>(6,0,0)+shift;
    rigid_body->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

    id2=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"subdivided_box");
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id2);
    rigid_body->frame.t=VECTOR<T,3>(0,0,0)+shift;
    rigid_body->frame.r=QUATERNION<T>(pi/2,VECTOR<T,3>(1,0,0));
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(1);
    rigid_body->Set_Name("child2");
    rigid_body->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);


    ARRAY<RIGID_BODY<TV>*>* bodies=new ARRAY<RIGID_BODY<TV>*>();
    bodies->Append(solids_parameters.rigid_body_parameters.list(id1));
    bodies->Append(solids_parameters.rigid_body_parameters.list(id2));
    cluster1=solids_parameters.rigid_body_parameters.list.Add_Cluster_Body(bodies,solids_parameters.collision_body_list);
    arb->rigid_bodies_list(cluster1)->Set_Name("combo_square");


}
//#####################################################################
// Simple Spring
//#####################################################################
void Simple_Spring()
{
    // CLUSTER 1
    int id3,id4;
    id3=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"subdivided_box");
    speed_up_test1=arb->rigid_bodies_list.rigid_bodies(id3);
    speed_up_test1->frame.t=VECTOR<T,3>(0,5,0);
    speed_up_test1->Set_Coefficient_Of_Restitution(0.5);
    speed_up_test1->Set_Coefficient_Of_Friction(0.5);
    speed_up_test1->Set_Name("parent");
    speed_up_test1->is_static=true;

    id4=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"subdivided_box");
    speed_up_test2=arb->rigid_bodies_list.rigid_bodies(id4);
    speed_up_test2->frame.t=VECTOR<T,3>(0,2,0);
    speed_up_test2->Set_Coefficient_Of_Restitution(0.5);
    speed_up_test2->Set_Coefficient_Of_Friction(.5);
    speed_up_test2->Set_Name("child");

/*
    ARRAY<VECTOR<int,2> > segments(0);

    soft_constraints.Add_Constrained_Point(new CONSTRAINED_POINT_IN_RIGID_BODY<T,VECTOR<T,3> >(*speed_up_test1,VECTOR<T,3>(0,0,0)));
    soft_constraints.Add_Constrained_Point(new CONSTRAINED_POINT_IN_RIGID_BODY<T,VECTOR<T,3> >(*speed_up_test2,VECTOR<T,3>(0,0,0)));
    segments.Append(1,2);

    SEGMENTED_CURVE<T,VECTOR<T,3> > & segmented_curve = soft_constraints.deformable_object.template Find_Structure<SEGMENTED_CURVE<T,VECTOR<T,3> >&>(); 
    segmented_curve.mesh.Initialize_Mesh(2,segments);segmented_curve.mesh.Initialize_Connected_Segments();
    soft_constraints.Synchronize_Particles_With_Constrained_Points();
    soft_constraints.deformable_object.particles.mass.Compute_Effective_Mass(soft_constraints.solid_body_collection.deformable_body_collection.binding_list);
    LINEAR_SPRINGS<T,VECTOR<T,3> >* ls=Create_Edge_Springs(segmented_curve.mesh,soft_constraints.deformable_object.particles,(T)2e4);
    soft_constraints.solid_body_collection.Add_Force(ls);
    ls->Clamp_Restlength(.1);
    soft_constraints.solid_body_collection.Update_Fragments();
*/
}
//#####################################################################
// Apply_Constraints
//#####################################################################
void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE
{
    std::cout << "Applying constraint " << std::endl;
//    soft_constraints.Apply_Constraint_Correcting_Impulses(dt,time);
}
//#####################################################################
// Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{

    //if(frame==35)  solids_parameters.rigid_body_parameters.list.Remove_Cluster_Body(cluster1,solids_parameters.collision_body_list);
    
  //if(frame==45) solids_parameters.rigid_body_parameters.list.Remove_Cluster_Body(cluster2);
}
//#####################################################################
// Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
}
//#####################################################################
};
}
#endif

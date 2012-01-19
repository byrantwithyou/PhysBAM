//#####################################################################
// Copyright 2006, Mike Rodgers, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SIMPLE_MUSCLE_EXAMPLE
//##################################################################### 
#ifndef __SIMPLE_MUSCLE_EXAMPLE__
#define __SIMPLE_MUSCLE_EXAMPLE__

#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/PARTICLE_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TETRAHEDRON_COLLISION_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_S3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/ANGLE_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE_FORCE_CURVE.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Dynamics/Motion/FRAME_TRACK_3D.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <fstream>
#include "../ARB_PARAMETERS.h"

namespace PhysBAM{

template<class T,class RW>
class SIMPLE_MUSCLE_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef VECTOR<T,3> TV;

    ARRAY<ARRAY<int> > muscle_tets;
    ARRAY<ARRAY<VECTOR<T,3> > > muscle_fibers;
    ARRAY<ARRAY<T> > muscle_densities;
    ARRAY<T> muscle_activations;

    enum GEOMETRY_TYPE {TETRAHEDRALIZED_VOLUME_TYPE,TRIANGULATED_SURFACE_TYPE,EMBEDDED_TETRAHEDRALIZED_VOLUME_TYPE};
    int number_of_deformable_bodies;
    ARRAY<std::string> deformable_body_geometry_filenames;
    ARRAY<GEOMETRY_TYPE> deformable_body_geometry_types;
    ARRAY<RIGID_BODY_STATE<TV>* > deformable_body_initial_states;
    ARRAY<ARRAY<TV> > deformable_body_rest_positions;
    bool use_ground_plane;

    VECTOR<int,3> mattress_dimensions;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    // typedef CONSTRAINED_POINT_IN_RIGID_BODY<T,VECTOR<T,3>,RIGID_BODY<TV> > T_CONSTRAINED_POINT_IN_RIGID_BODY;    
    typedef CONSTRAINED_POINT_IN_RIGID_BODY<T,VECTOR<T,3> > T_CONSTRAINED_POINT_IN_RIGID_BODY;

    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;
    using BASE::restart;using BASE::restart_frame;using BASE::output_directory;
    using BASE::solids_parameters;using BASE::write_last_frame;using BASE::data_directory;

    ARTICULATED_RIGID_BODY<TV>* arb;
    MUSCLE<TV>* muscle;
    bool add_ground;
    PARAMETER_LIST parameter_list;
    T peak_force;


    // Data pertaining to enslaving flesh to bone
    int num_planks;
    ARRAY<int> plank_ids;
    ARRAY<ARRAY<int> > enslaved_nodes; // indicies of enslaved nodes. 
    ARRAY<ARRAY<VECTOR<T,3> > > positions_relative_to_plank_frames; // corresponding positions


    SIMPLE_MUSCLE_EXAMPLE()
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>(0,FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >::NONE),number_of_deformable_bodies(0),use_ground_plane(false),mattress_dimensions(5,35,6)/*(3,7,3)*/
    {
        // Mattress is 5x43x6 spatially
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.cfl=(T).1;
        solids_parameters.collisions_repulsion_thickness=(T)1e-2;
        solids_parameters.collisions_repulsion_clamp_fraction=(T).9;
        solids_parameters.collision_repulsion_spring_multiplier=100;
        solids_parameters.collisions_repulsion_spring_constant_over_mass_times_length=30;    
        solids_parameters.collisions_nonrigid_collision_attempts=3;
        solids_parameters.collide_with_interior=false;
        solids_parameters.check_initial_mesh_for_self_intersection=false;
        solids_parameters.collisions_disable_repulsions_based_on_proximity_factor=1.5;
        solids_parameters.quasistatic=true;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        solids_parameters.implicit_solve_parameters.cg_iterations=900;
        solids_parameters.newton_iterations=10;
        solids_parameters.newton_tolerance=(T)1e-2;
        solids_parameters.use_partially_converged_result=true;
        solids_parameters.collision_tolerance=(T)1e-4; 
        solids_parameters.perform_self_collision=false;
        solids_parameters.perform_collision_body_collisions=false;

        Add_Deformable_Body("",TETRAHEDRALIZED_VOLUME_TYPE,0);

        last_frame=3000;
        frame_rate=24;
        output_directory="Simple_Muscle/output";
        std::cout << "Frame rate: "<<frame_rate<<std::endl;

        arb=new ARTICULATED_RIGID_BODY<TV>(this->solids_parameters.rigid_body_parameters.list);
        this->solids_parameters.rigid_body_parameters.Set_Articulated_Rigid_Body(arb);

        arb->Set_Do_Final_Pass(false);
        arb->Use_Epsilon_Scale(false);
        arb->Set_Use_Shock_Propagation(false);

        arb->Use_Muscle_Actuators();

        muscle=0;
        
        ARB_PARAMETERS::Read_Common_Parameters("Simple_Muscle/example.param",*this,parameter_list);
        peak_force=parameter_list.Get_Parameter("peak_force",(T)10);
    }

    ~SIMPLE_MUSCLE_EXAMPLE()
    {
        for(int i=0;i<deformable_body_initial_states.m;i++)
            delete deformable_body_initial_states(i);
    }

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    // Initialize deformable bodies
    Get_Initial_Data();

    DEFORMABLE_OBJECT<TV>& deformable_object=solid_body_collection.deformable_object;
    PARTICLES<TV>& particles=deformable_object.particles;

    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_object.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
    solid_body_collection.deformable_object.Add_Force(Create_Quasistatic_Diagonalized_Finite_Volume(tetrahedralized_volume,new DIAGONALIZED_NEO_HOOKEAN<T,3>((T)1e5,(T).45)));

    solid_body_collection.deformable_object.Update_Fragments();


    // Initialize rigid bodies and muscles
    arb->muscle_list->muscle_force_curve.Initialize(data_directory); // initialize here rather than constructor since data directory might be set after constructor

    add_ground=true;
    //Arm();
    Simple_Muscle_Across_Joint();
    //Stand_Test();

    if(add_ground){
        int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/ground",(T).1);
        RIGID_BODY<TV>* rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
        rigid_body->frame.t=VECTOR<T,3>(0,0,0);
        rigid_body->Set_Coefficient_Of_Restitution(0.5);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name("ground");
        rigid_body->is_static=true;
        rigid_body->add_to_spatial_partition=false;}

    std::cout <<"\nnumber of muscles is: "<<arb->muscle_list->muscles.m<<std::endl;

    RIGID_BODY_LIST<T,VECTOR<T,3> >& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    tests.Add_Gravity();

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Initialize_Bodies();

    Get_Constrained_Particle_Data();
}
//#####################################################################
// Function Simple_Muscle_Across_Joint
//#####################################################################
MUSCLE<TV>* Add_Basic_Muscle(const std::string& name,RIGID_BODY<TV>& origin_body,const VECTOR<T,3>& origin,RIGID_BODY<TV>& insertion_body,const VECTOR<T,3>& insertion)
{
    MUSCLE<TV>* muscle=new MUSCLE<TV>(arb->muscle_list->muscle_force_curve);
    muscle->Set_Name(name);
    muscle->Set_Attachment_Point_1(new T_CONSTRAINED_POINT_IN_RIGID_BODY(origin_body,origin));
    muscle->Set_Attachment_Point_2(new T_CONSTRAINED_POINT_IN_RIGID_BODY(insertion_body,insertion));
    T total_length=muscle->Total_Length();
    muscle->Set_Optimal_Length((T).8*total_length);
    muscle->Set_Tendon_Slack_Length((T).2*total_length);
    muscle->Set_Peak_Force(peak_force);
    muscle->Set_Max_Shortening_Velocity(1);
    arb->muscle_list->Add_Muscle(muscle);
    return muscle;
}

void Simple_Muscle_Across_Joint()
{
    RIGID_BODY_LIST<T,VECTOR<T,3> >& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    assert(rigid_body_list.Number_Of_Elements()==0);

    add_ground=false;
    solids_parameters.gravity=0;

    int num_bodies=num_planks=parameter_list.Get_Parameter("num_bodies",(int)2);
    plank_ids.Resize(num_planks);
    bool use_bend_joint=parameter_list.Get_Parameter("use_bend_joint",false);
    VECTOR<T,3> plank_rescale=parameter_list.Get_Parameter("plank_rescale",VECTOR<T,3>(1,1,1));
    T target_angle=parameter_list.Get_Parameter("target_angle",-(T)pi/2);
    T initial_angle=parameter_list.Get_Parameter("initial_angle",-(T)0);
    T k_p=parameter_list.Get_Parameter("k_p",(T)100);

    FRAME_TRACK_3D<T>* frame_track=0;
    if(parameter_list.Get_Parameter("use_frame_track",false)){
        int samples=1000;T period=10;
        frame_track=new FRAME_TRACK_3D<T>(samples,0,period);frame_track->periodic=true;
        for(int i=0;i<samples;i++){
            frame_track->trajectory(i)=FRAME_3D<T>(QUATERNION<T>(initial_angle+(target_angle-initial_angle)*0.5*(1-cos(2*pi*(i-1)/(samples-1))),VECTOR<T,3>(1,0,0)));}
        for(int i=0;i<arb->joint_mesh.joints.m;i++) arb->joint_mesh.joints(i)->joint_function->track=frame_track;
    }

    
    for(int i=0;i<num_bodies;i++){
        int id=plank_ids(i)=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/plank",(T).2);
        RIGID_BODY<TV>* rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
        rigid_body->triangulated_surface->Rescale(plank_rescale.x,plank_rescale.y,plank_rescale.z);
        rigid_body->frame.r=QUATERNION<T>(pi/2,VECTOR<T,3>(0,1,0));
        rigid_body->Set_Coefficient_Of_Restitution(0.5);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name(STRING_UTILITIES::string_sprintf("body_%d",i));
        rigid_body->Set_Mass(1);

        if(i>1){
            JOINT<TV>* joint=0;
            if(use_bend_joint) joint=new ANGLE_JOINT<TV>();
            else joint=new POINT_JOINT<TV>();
            joint->Set_Joint_To_Parent_Frame(FRAME_3D<T>(VECTOR<T,3>(0,0,1.05)));
            joint->Set_Joint_To_Child_Frame(FRAME_3D<T>(VECTOR<T,3>(0,.0,-1.05)));
            int joint_id=arb->joint_mesh.Add_Joint(joint);
            arb->Add_Articulation(i-1,i,joint_id);
            JOINT_FUNCTION<TV>* jfunc=new JOINT_FUNCTION<TV>(joint,rigid_body_list(i-1),rigid_body_list(i));
            joint->Set_Joint_Function(jfunc);
            if(frame_track) jfunc->track=frame_track;
            else jfunc->Set_Target_Angle(QUATERNION<T>(target_angle,VECTOR<T,3>(1,0,0)));
            //jfunc->Set_Target_Angle(QUATERNION<T>::From_Euler_Angles((T)-pi/2,0.5,0));
            jfunc->muscle_control=true;
            joint->joint_function->Set_k_p(k_p);
            joint->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(initial_angle,VECTOR<T,3>(1,0,0))));
        }
    }

    rigid_body_particles.Rigid_Body(1).frame.r=QUATERNION<T>(-pi/2,VECTOR<T,3>(0,0,1))*rigid_body_particles.Rigid_Body(1).frame.r;


    for(int i=2;i<=num_bodies;i++){
        std::string suffix=STRING_UTILITIES::string_sprintf("_j%d",i-1);
        RIGID_BODY<TV> *parent=rigid_body_list(i-1),*child=rigid_body_list(i);

        Add_Basic_Muscle("flexor"+suffix,*parent,VECTOR<T,3>(0,0.05,-0.3),*child,VECTOR<T,3>(0,0.05,-0.2));
        if(parameter_list.Get_Parameter("add_extensor",true)){
            MUSCLE<TV> *extensor=Add_Basic_Muscle("extensor"+suffix,*parent,VECTOR<T,3>(0,-0.05,0.5),*child,VECTOR<T,3>(0,-0.05,-0.5));
            extensor->Add_Via_Point(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*parent,VECTOR<T,3>(0,-0.2,1)));
            extensor->Add_Via_Point(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*child,VECTOR<T,3>(0,-0.2,-1)));}
        if(parameter_list.Get_Parameter("add_redundant_muscles",false)){ // EXTRA MUSCLES
            Add_Basic_Muscle("redundant_flexor"+suffix,*parent,VECTOR<T,3>(0,0.05,0.5),*child,VECTOR<T,3>(0,0.05,0.3));
            MUSCLE<TV> *redundant_extensor=Add_Basic_Muscle("redundant_extensor"+suffix,*parent,VECTOR<T,3>(0,-0.05,-0.5),*child,VECTOR<T,3>(0,-0.05,0.5));
            redundant_extensor->Add_Via_Point(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*parent,VECTOR<T,3>(0,-0.5,1)));
            redundant_extensor->Add_Via_Point(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*child,VECTOR<T,3>(0,-0.4,-1)));}
        if(parameter_list.Get_Parameter("add_diagonal_muscles",false)){
            // MUSCLE<TV> *diag1=Add_Basic_Muscle("diag1"+suffix,*parent,VECTOR<T,3>(0.2,0.05,-0.5),*child,VECTOR<T,3>(-0.2,0.05,-0.5));
            // MUSCLE<TV> *diag2=Add_Basic_Muscle("diag2"+suffix,*parent,VECTOR<T,3>(-0.2,0.05,-0.5),*child,VECTOR<T,3>(0.2,0.05,-0.5));
            MUSCLE<TV> *side1=Add_Basic_Muscle("side1"+suffix,*parent,VECTOR<T,3>(0.2,0,-0.5),*child,VECTOR<T,3>(0.2,0,0.7));
            // side1->Add_Via_Point(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*parent,VECTOR<T,3>(0.4,0,0.95)));
            side1->Add_Via_Point(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*child,VECTOR<T,3>(0.4,0,-0.95)));
            MUSCLE<TV> *side2=Add_Basic_Muscle("side2"+suffix,*parent,VECTOR<T,3>(-0.2,0,-0.5),*child,VECTOR<T,3>(-0.2,0,0.7));
            // side2->Add_Via_Point(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*parent,VECTOR<T,3>(-0.4,0,0.95)));
            side2->Add_Via_Point(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*child,VECTOR<T,3>(-0.4,0,-0.95)));
        }
    }

    if(parameter_list.Get_Parameter("add_multijoint_muscle",false)){
//        MUSCLE<TV> *muscle=Add_Basic_Muscle("multijoint",*rigid_body_list(1),VECTOR<T,3>(0,0.05,-0.5),*rigid_body_list(rigid_body_list.rigid_bodies.m),VECTOR<T,3>(0,0.05,-0.5));
    }

    arb->Update_With_Breadth_First_Directed_Graph(1);

    if(parameter_list.Get_Parameter("add_obstacle",false)){
        int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/box",0.3);
        RIGID_BODY<TV>* rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
        rigid_body->frame.t=VECTOR<T,3>(1,3,0);
        rigid_body->Set_Coefficient_Of_Restitution(0.5);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name("obstacle");
        rigid_body->Set_Mass(1);
        rigid_body->is_static=true;}
}
//#####################################################################
// Function Stand_Test
//#####################################################################
void Stand_Test()
{
    RIGID_BODY_LIST<T,VECTOR<T,3> >& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    assert(rigid_body_list.Number_Of_Elements()==0);

    add_ground=parameter_list.Get_Parameter("add_ground",false);
    int num_bodies=parameter_list.Get_Parameter("num_bodies",(int)2);
    T k_p=parameter_list.Get_Parameter("k_p",(T)100);

    for(int i=0;i<num_bodies;i++){
        int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/plank",(T).2);
        RIGID_BODY<TV>* rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
        rigid_body->frame.r=QUATERNION<T>(pi/2,VECTOR<T,3>(0,1,0));
        rigid_body->Set_Coefficient_Of_Restitution(0.5);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name(STRING_UTILITIES::string_sprintf("body_%d",i));
        rigid_body->Set_Mass(1);

        if(i>1){
            JOINT<TV>* joint=new ANGLE_JOINT<TV>();
            joint->Set_Joint_To_Parent_Frame(FRAME_3D<T>(VECTOR<T,3>(0,0,1)));
            joint->Set_Joint_To_Child_Frame(FRAME_3D<T>(VECTOR<T,3>(0,0,-1)));
            int joint_id=arb->joint_mesh.Add_Joint(joint);
            arb->Add_Articulation(i-1,i,joint_id);
            JOINT_FUNCTION<TV>* jfunc=new JOINT_FUNCTION<TV>(joint,rigid_body_list(i-1),rigid_body_list(i));
            joint->Set_Joint_Function(jfunc);
            jfunc->Set_Target_Angle(QUATERNION<T>(-(T)3*pi/4,VECTOR<T,3>(1,0,0)));
            //jfunc->Set_Target_Angle(QUATERNION<T>::From_Euler_Angles((T)-pi/2,0.5,0));
            jfunc->muscle_control=true;
            joint->joint_function->Set_k_p(k_p);
            joint->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(-(T)pi/2,VECTOR<T,3>(1,0,0))));
        }
    }

    bool use_track=parameter_list.Get_Parameter("use_track",false);

    if(num_bodies==5){
        if(parameter_list.Get_Parameter("push_up",false)){
            rigid_body_particles.Rigid_Body(1).frame.r=QUATERNION<T>(3*pi/4,VECTOR<T,3>(0,0,1))*rigid_body_particles.Rigid_Body(1).frame.r;
            rigid_body_particles.Rigid_Body(1).frame.t=VECTOR<T,3>(0,.8,0);

            rigid_body_particles.Rigid_Body(3).Set_Mass(parameter_list.Get_Parameter("mass_scale",(T)100)*rigid_body_particles.Rigid_Body(3).mass);

            arb->joint_mesh.joints(1)->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(pi/2,VECTOR<T,3>(1,0,0))));
            arb->joint_mesh.joints(2)->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(-3*pi/4,VECTOR<T,3>(1,0,0))));
            arb->joint_mesh.joints(3)->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(-pi/4,VECTOR<T,3>(1,0,0))));
            arb->joint_mesh.joints(4)->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(0,VECTOR<T,3>(1,0,0))));

            //arb->joint_mesh.joints(1)->joint_function->Set_Target_Angle(QUATERNION<T>(pi/2,VECTOR<T,3>(1,0,0)));
            //arb->joint_mesh.joints(2)->joint_function->Set_Target_Angle(QUATERNION<T>(-3*pi/4,VECTOR<T,3>(1,0,0)));
            //arb->joint_mesh.joints(3)->joint_function->Set_Target_Angle(QUATERNION<T>(-pi/4,VECTOR<T,3>(1,0,0)));
            //arb->joint_mesh.joints(4)->joint_function->Set_Target_Angle(QUATERNION<T>(0,VECTOR<T,3>(1,0,0)));

            arb->joint_mesh.joints(1)->joint_function->Set_Target_Angle(QUATERNION<T>(0,VECTOR<T,3>(1,0,0)));
            arb->joint_mesh.joints(2)->joint_function->Set_Target_Angle(QUATERNION<T>(-pi/2,VECTOR<T,3>(1,0,0)));
            arb->joint_mesh.joints(3)->joint_function->Set_Target_Angle(QUATERNION<T>(-pi/4,VECTOR<T,3>(1,0,0)));
            arb->joint_mesh.joints(4)->joint_function->Set_Target_Angle(QUATERNION<T>(0,VECTOR<T,3>(1,0,0)));
        }
        else{
            rigid_body_particles.Rigid_Body(1).frame.r=QUATERNION<T>(3*pi/4,VECTOR<T,3>(0,0,1))*rigid_body_particles.Rigid_Body(1).frame.r;
            rigid_body_particles.Rigid_Body(1).frame.t=VECTOR<T,3>(0,.8,0);

            rigid_body_particles.Rigid_Body(3).Set_Mass(parameter_list.Get_Parameter("mass_scale",(T)100)*rigid_body_particles.Rigid_Body(3).mass);

            arb->joint_mesh.joints(1)->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(0,VECTOR<T,3>(1,0,0))));
            arb->joint_mesh.joints(2)->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(-pi/4,VECTOR<T,3>(1,0,0))));
            arb->joint_mesh.joints(3)->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(-pi/4,VECTOR<T,3>(1,0,0))));
            arb->joint_mesh.joints(4)->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(0,VECTOR<T,3>(1,0,0))));

            arb->joint_mesh.joints(1)->joint_function->Set_Target_Angle(QUATERNION<T>(0,VECTOR<T,3>(1,0,0)));
            arb->joint_mesh.joints(2)->joint_function->Set_Target_Angle(QUATERNION<T>(-pi/4,VECTOR<T,3>(1,0,0)));
            arb->joint_mesh.joints(3)->joint_function->Set_Target_Angle(QUATERNION<T>(-pi/4,VECTOR<T,3>(1,0,0)));
            arb->joint_mesh.joints(4)->joint_function->Set_Target_Angle(QUATERNION<T>(0,VECTOR<T,3>(1,0,0)));

            //arb->joint_mesh.joints(1)->joint_function->Set_Target_Angle(QUATERNION<T>(0,VECTOR<T,3>(1,0,0)));
            //arb->joint_mesh.joints(2)->joint_function->Set_Target_Angle(QUATERNION<T>(-pi/2,VECTOR<T,3>(1,0,0)));
            //arb->joint_mesh.joints(3)->joint_function->Set_Target_Angle(QUATERNION<T>(-pi/4,VECTOR<T,3>(1,0,0)));
            //arb->joint_mesh.joints(4)->joint_function->Set_Target_Angle(QUATERNION<T>(0,VECTOR<T,3>(1,0,0)));
        }
    }
    else if(num_bodies==3){
        rigid_body_particles.Rigid_Body(1).frame.r=QUATERNION<T>(pi,VECTOR<T,3>(0,0,1))*rigid_body_particles.Rigid_Body(1).frame.r;
        rigid_body_particles.Rigid_Body(1).frame.t=VECTOR<T,3>(0,.05,0);

        rigid_body_particles.Rigid_Body(2).Set_Mass(parameter_list.Get_Parameter("mass_scale",(T)100)*rigid_body_particles.Rigid_Body(3).mass);

        arb->joint_mesh.joints(1)->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(0,VECTOR<T,3>(1,0,0))));
        arb->joint_mesh.joints(2)->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(0,VECTOR<T,3>(1,0,0))));

        arb->joint_mesh.joints(1)->joint_function->Set_Target_Angle(QUATERNION<T>(-pi/4,VECTOR<T,3>(1,0,0)));
        arb->joint_mesh.joints(2)->joint_function->Set_Target_Angle(QUATERNION<T>(-pi/4,VECTOR<T,3>(1,0,0)));
    }

    arb->Update_With_Breadth_First_Directed_Graph(1);

    if(parameter_list.Get_Parameter("add_obstacle",false)){
        int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/box",0.3);
        RIGID_BODY<TV>* rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
        rigid_body->frame.t=VECTOR<T,3>(1,3,0);
        rigid_body->Set_Coefficient_Of_Restitution(0.5);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name("obstacle");
        rigid_body->Set_Mass(1);
        rigid_body->is_static=true;
    }
}
//#####################################################################
// Function Arm
//#####################################################################
void Arm()
{
    //PRISMATIC_JOINT_2D<T>* joint=new PRISMATIC_JOINT_2D<T>();arb->joint_mesh.Add_Joint(joint);
    //joint->Set_Prismatic_Component_Bounds(VECTOR_2D<T>(0,-.5),VECTOR_2D<T>(0,.5));
    JOINT<TV>* joint=new POINT_JOINT<TV>();arb->joint_mesh.Add_Joint(joint);
    std::cout<<" building ARM\n";
    RIGID_BODY<TV> *rigid_body1 = 0;
    RIGID_BODY<TV> *rigid_body2 = 0;

    int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/box",1);
    rigid_body1=arb->rigid_bodies_list.rigid_bodies(id);
    rigid_body1->frame.t=VECTOR<T,3>(0,4,0);
    rigid_body1->Set_Coefficient_Of_Restitution(0.5);
    rigid_body1->Set_Coefficient_Of_Friction(0);
    rigid_body1->Set_Name("square1");
    rigid_body1->Set_Mass(1);
    joint->Set_Joint_To_Parent_Frame(FRAME_3D<T>(VECTOR<T,3>(0,2,0)));
    
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/box",1);
    rigid_body2=arb->rigid_bodies_list.rigid_bodies(id);
    rigid_body2->frame.t=VECTOR<T,3>(0,(T)6.7,0);
    rigid_body2->Set_Coefficient_Of_Friction(0);
    rigid_body2->Set_Coefficient_Of_Restitution((T)0.5);
    rigid_body2->Set_Name("square2");
    rigid_body2->is_static=true;
    joint->Set_Joint_To_Child_Frame(FRAME_3D<T>(VECTOR<T,3>(0,-1,0)));

    muscle=new MUSCLE<TV>(arb->muscle_list->muscle_force_curve);
    muscle->Set_Attachment_Point_1(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*rigid_body1,VECTOR<T,3>(0,1,0))); //don't think these are exactly right . . .
    muscle->Set_Attachment_Point_2(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*rigid_body2,VECTOR<T,3>(0,-1,0))); //don't think these are exactly right . . .
    muscle->Set_Optimal_Length((T).5);
    muscle->Set_Peak_Force(peak_force);
    muscle->Set_Max_Shortening_Velocity(1);
    arb->muscle_list->Add_Muscle(muscle);

    Add_Basic_Muscle("test",*rigid_body1,VECTOR<T,3>(0,1,0),*rigid_body2,VECTOR<T,3>(0,-1,0));

    JOINT_FUNCTION<TV>* jfunc=new JOINT_FUNCTION<TV>(joint,rigid_body1,rigid_body2);
    joint->Set_Joint_Function(jfunc);
//    jfunc->Set_Desired_Angle(-.3);
    jfunc->Set_Target_Angle(QUATERNION<T>(-.3,VECTOR<T,3>(1,0,0)));
    joint->joint_function->Set_k_p(10);
    jfunc->muscle_control=true;

    arb->Add_Articulation(1,2,1);
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    std::cout<<"PREPROCESSING\n"<<std::endl;
    //T activation=fabs(muscle->goal_length/muscle->optimal_length-muscle->Length()/muscle->optimal_length);
    //if(activation>1) activation=1;
    //if(muscle) arb->muscle_list->Set_Muscle_Activation(muscle->id,muscle->Calculate_Activation(muscle2->Force(1)));
    
}
void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Write_Output_Files(frame);
#if 0
    const RIGID_BODY_LIST<T,VECTOR<T,3> >& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    VECTOR<T,3> total_linear_momentum,total_angular_momentum;
    for(int i=0;i<rigid_body_list.Number_Of_Elements();i++){
        LOG::cout << i << ": " << rigid_body_particles.Rigid_Body(i).Angular_Momentum() << std::endl;
        total_linear_momentum+=rigid_body_particles.Rigid_Body(i).mass*rigid_body_particles.Rigid_Body(i).velocity;
        total_angular_momentum+=VECTOR<T,3>::Cross_Product(rigid_body_particles.Rigid_Body(i).frame.t,rigid_body_particles.Rigid_Body(i).mass*rigid_body_particles.Rigid_Body(i).velocity)+rigid_body_particles.Rigid_Body(i).Angular_Momentum();}
    LOG::cout << "MOMENTA === linear " << total_linear_momentum << ", angular " << total_angular_momentum << std::endl;
#endif
}
//#####################################################################


// Begin what was originally class SOLIDS_STANDARD_TESTS_3D
public:
    TV attachment_velocity;

//#####################################################################
// Function Add_Deformable_Body
//#####################################################################
void Add_Deformable_Body(const std::string& filename,const GEOMETRY_TYPE type,RIGID_BODY_STATE<TV>* initial_state)
{
    number_of_deformable_bodies++;
    deformable_body_geometry_filenames.Append(filename);deformable_body_geometry_types.Append(type);deformable_body_initial_states.Append(initial_state);
}
//#####################################################################
// Function Set_Initial_Particle_Configuration
//#####################################################################
void Set_Initial_Particle_Configuration(PARTICLES<TV>& particles,const int index)
{
    if(deformable_body_initial_states(index)){
        LOG::cout<<"Deformable body "<<index<<" - Total Particles : "<<particles.array_collection->Size()<<std::endl;
        BOX_3D<T> bounding_box(particles.X(1));for(int i=2;i<=particles.array_collection->Size();i++) bounding_box.Enlarge_To_Include_Point(particles.X(i));TV center=bounding_box.Center();
        RIGID_BODY_STATE<TV>& state=*deformable_body_initial_states(index);
        for(int p=0;p<particles.array_collection->Size();p++){
            particles.X(p)=state.frame*(particles.X(p)-center);
            particles.V(p)=state.velocity+TV::Cross_Product(state.angular_velocity,particles.X(p)-state.frame.t);}}
}
//#####################################################################
// Function Create_Tetrahedralized_Volume
//#####################################################################
STRUCTURE<TV>* Create_Tetrahedralized_Volume(int index)
{
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*TETRAHEDRALIZED_VOLUME<T>::Create();
    PARTICLES<TV>& particles=tetrahedralized_volume.particles;
    GRID<TV> mattress_grid(mattress_dimensions.x,mattress_dimensions.y,mattress_dimensions.z,(T)-.25,(T).25,(T)-3.2,(T)1.10,(T)-.30,(T).30);
    tetrahedralized_volume.Initialize_Cube_Mesh_And_Particles(mattress_grid);
    deformable_body_rest_positions(index)=particles.X.array;
    LOG::cout<<"Deformable body "<<index<<" - Total Tetrahedra : "<<tetrahedralized_volume.mesh.elements.m<<std::endl;
    particles.Update_Velocity();particles.Store_Mass();
    tetrahedralized_volume.Set_Density(1000);tetrahedralized_volume.Set_Mass_Of_Particles(true);
    Set_Initial_Particle_Configuration(particles,index);
    return &tetrahedralized_volume;
}
//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
void Initialize_Rigid_Bodies()
{
    if(use_ground_plane){
        int ground_index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/ground");
        solids_parameters.rigid_body_parameters.list.rigid_bodies(ground_index)->Set_Coefficient_Of_Friction((T).3);
        solids_parameters.rigid_body_parameters.list.rigid_bodies(ground_index)->is_static=true;
        solids_parameters.rigid_body_parameters.list.rigid_bodies(ground_index)->add_to_spatial_partition=false;}
}
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    // deformable bodies
    DEFORMABLE_OBJECT<TV>& deformable_object=solid_body_collection;
    deformable_body_rest_positions.Resize(number_of_deformable_bodies);
    // initialize geometry
    STRUCTURE<TV>* structure=Create_Tetrahedralized_Volume(1);
    // add to deformable_object
    deformable_object.Add_Structure(structure->Append_Particles_And_Create_Copy(deformable_object.particles));
    deformable_object.collisions.collision_structures.Append(deformable_object.structures.Last());
    delete structure;

    // correct number nodes
    for(int i=0;i<deformable_object.structures.m;i++) deformable_object.structures(i)->Update_Number_Nodes();

    // Rigid Bodies
    Initialize_Rigid_Bodies();

    // Collisions
    solids_parameters.perform_self_collision=false;Initialize_Tetrahedron_Collisions(1,deformable_object.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>());                                                                                                                                                                
    // solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
    solids_parameters.perform_self_collision=false;
}
//#####################################################################
// Function Update_Collision_Body_Positions_And_Velocities
//#####################################################################
void Update_Collision_Body_Positions_And_Velocities(const T time) PHYSBAM_OVERRIDE{
    BASE::Update_Collision_Body_Positions_And_Velocities(time);
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<TV> X,const T time)
{
    assert(id_number==1);
    for(int p=0;p<num_planks;p++){
        RIGID_BODY<TV>& plank_rigid_body=*arb->rigid_bodies_list.rigid_bodies(plank_ids(p));
        FRAME_3D<T> inverted_frame=plank_rigid_body.frame.Inverse();
        for(int i=0;i<enslaved_nodes(p).m;i++){
            X(enslaved_nodes(p)(i))=inverted_frame.Local_Coordinate(positions_relative_to_plank_frames(p)(i));}}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time)
{
    assert(fragment_id==FRAGMENT_ID(1));
    for(int p=0;p<num_planks;p++)for(int i=0;i<enslaved_nodes(p).m;i++){
        X(enslaved_nodes(p)(i))=TV();} 
}
//#####################################################################
// Function Add_External_Impulses
//#####################################################################
void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt)
{} // TODO: Overload these for asynchronous
//#####################################################################
// Function Add_External_Impulse
//#####################################################################
void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt)
{} // TODO: Overload these for asynchronous
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
void Update_Time_Varying_Material_Properties(const T time)
{
    if (muscle_activations.m < 1)
        return;

}
//#####################################################################
// Function Get_Constrained_Particle_Data
//#####################################################################
void Get_Constrained_Particle_Data()
{
    PARTICLES<TV>& particles=(solid_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>()).particles;
 
    enslaved_nodes.Resize(num_planks);
    positions_relative_to_plank_frames.Resize(num_planks);

    for(int i=1;i<=mattress_dimensions.x*mattress_dimensions.y*mattress_dimensions.z;i++){
        for(int p=0;p<num_planks;p++){
            RIGID_BODY<TV>& plank_rigid_body=*arb->rigid_bodies_list.rigid_bodies(plank_ids(p));
            if (!plank_rigid_body.Implicit_Geometry_Lazy_Outside(particles.X(i))) {
                enslaved_nodes(p).Append(i);
                positions_relative_to_plank_frames(p).Append(plank_rigid_body.frame.Local_Coordinate(particles.X(i)));}}}
}
//#####################################################################
};
}
#endif

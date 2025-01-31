//#####################################################################
// Copyright 2007-2008, Micahel Lentine, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//   1. Point joint between 2 blocks
//   2. Rigid joint between 2 blocks
//   3. Hinge (bend) joint between 2 blocks
//   4. Twist joint between 2 blocks
//   5. String of bodies connected by rigid joints, colliding with ground
//   6. Closed loop of nonconvex objects
//   7. Planks PD controlled to curl to form an "O" shape
//   8. An 8-block cluster and mimicking solid block of same size
//   9. An 8-block cluster breaking
//  10. Prismatic joint between 2 blocks
//  11. Constrained hinge
//  12. Heavy bottom link
//  13. Universal joint
//  14. Pre-stabilization test system
//  15. Test 2 synthesized from 3 point joints
//  16. Test 4 synthesized from 2 point joints
//  17. Joint with muscle
//  18. Muscle bend
//  19. Nested clusters
//  20. Post-stabilization test system
//  21. Sphere mesh
//  22. Overconstrained Joints
//  23. Clusters & springs
//  24. Clusters & springs with cluster breaking
//  25. Pendulum with point joint
//  26. Pendulum with angle joint
//  27. Normal joint
//  28. Angle joint with a kinematic object
//  29. Target angles
//  30. Push out and arb
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <Core/Log/DEBUG_UTILITIES.h>
#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <Tools/Ordinary_Differential_Equations/ODE_SOLVER.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER.h>
#include <Rigids/Joints/ANGLE_JOINT.h>
#include <Rigids/Joints/CONSTRAINT_FUNCTION.h>
#include <Rigids/Joints/JOINT_FUNCTION.h>
#include <Rigids/Joints/NORMAL_JOINT.h>
#include <Rigids/Joints/POINT_JOINT.h>
#include <Rigids/Joints/PRISMATIC_TWIST_JOINT.h>
#include <Rigids/Joints/RIGID_JOINT.h>
#include <Rigids/Muscles/ATTACHMENT_POINT.h>
#include <Rigids/Muscles/MUSCLE.h>
#include <Rigids/Muscles/MUSCLE_LIST.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Deformables/Bindings/BINDING_LIST.h>
#include <Deformables/Bindings/RIGID_BODY_BINDING.h>
#include <Deformables/Bindings/SOFT_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/FINITE_VOLUME.h>
#include <Deformables/Forces/LINEAR_SPRINGS.h>
#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS:public SOLIDS_EXAMPLE<VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_EXAMPLE<TV> BASE;
    typedef typename TV::SPIN T_SPIN;

    using BASE::viewer_dir;using BASE::solids_parameters;using BASE::data_directory;using BASE::last_frame;
    using BASE::stream_type;using BASE::frame_rate;using BASE::solid_body_collection;using BASE::test_number;
    using BASE::user_last_frame;
    
    SOLIDS_STANDARD_TESTS<TV> tests;
    int cluster_id;
    int kinematic_id;
    T peak_force;
    INTERPOLATION_CURVE<T,FRAME<TV> > curve;
    int parameter;
    bool use_prestab_iterations;

    struct NONLINEAR_PENDULUM:public NONLINEAR_FUNCTION<VECTOR<T,2>(T,VECTOR<T,2>)>
    {
        T g_over_l;
        VECTOR<T,2> y_not;

        VECTOR<T,2> operator()(const T t,const VECTOR<T,2> v) const override
        {return VECTOR<T,2>(v.y,-g_over_l*sin(v.x));}
    };
    NONLINEAR_PENDULUM nonlinear_pendulum;

    STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :BASE(stream_type_input,parse_args),tests(stream_type_input,data_directory,solid_body_collection),peak_force(10),parameter(3),use_prestab_iterations(false)
    {
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;

        solids_parameters.cfl=(T).9;

        solids_parameters.rigid_body_collision_parameters.use_push_out=true;
        solids_parameters.use_rigid_deformable_contact=true;
        ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
        parse_args.Add("-parameter",&parameter,"value","general use parameter");
        parse_args.Add("-print_energy",&solid_body_collection.print_energy,"print energy statistics");
        parse_args.Add_Not("-disable_prestab",&arb.use_prestab,"disable prestabilization");
        parse_args.Add_Not("-disable_poststab",&arb.use_poststab,"disable poststabilization");
        parse_args.Add("-use_krylov_poststab",&arb.use_krylov_poststab,"use krylov poststabilization");
        parse_args.Add("-use_krylov_prestab",&arb.use_krylov_prestab,"use krylov prestabilization");
        parse_args.Add("-test_arb_system",&solids_parameters.implicit_solve_parameters.test_system,"test arb system properties");
        parse_args.Add("-print_arb_matrix",&solids_parameters.implicit_solve_parameters.print_matrix,"print arb system");
        parse_args.Add("-prestab_iterations",&arb.max_iterations,&use_prestab_iterations,"iterations","prestabilization iterations");
        parse_args.Parse();

        if(!this->fixed_dt && !this->max_dt && !this->min_dt) this->fixed_dt=1;
        tests.data_directory=data_directory;
        if(!this->user_output_directory)
            viewer_dir.output_directory=LOG::sprintf("Standard_Tests/Test_%d",test_number);
        if(arb.use_krylov_poststab) arb.use_poststab_in_cg=false;
        if(use_prestab_iterations) solids_parameters.rigid_body_collision_parameters.contact_iterations=arb.max_iterations;
    }

    ~STANDARD_TESTS()
    {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    JOINT<TV>* joint=0;
    RIGID_BODY<TV>* rigid_body0=0,*rigid_body1=0;
    ROTATION<TV> rotation=ROTATION<TV>::From_Euler_Angles((T)79.9655*(T)pi/180,(T)12.0032*(T)pi/180,(T)50.4023*(T)pi/180);

    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    arb.Set_Iterative_Tolerance((T)1e-4);
    // prestabilization settings
    arb.Set_Contact_Level_Iterations(5);
    arb.Set_Shock_Propagation_Level_Iterations(5);
    arb.Set_Use_Shock_Propagation(false);
    arb.Set_Do_Final_Pass(false);
    // poststabilization settings
    arb.Set_Poststabilization_Iterations(5);
    arb.poststabilization_projection_iterations=2;
    if(!user_last_frame) last_frame=96;

    if(test_number<=4 || test_number==10 || test_number==11 || (test_number>=13 && test_number<=16) || test_number==18 || test_number==29){
        rigid_body0=&tests.Add_Rigid_Body("subdivided_box",1,(T).5);
        rigid_body1=&tests.Add_Rigid_Body("subdivided_box",1,(T).5);
        rigid_body0->Frame().t=TV(0,2,0);
        rigid_body0->name="parent";
        rigid_body1->name="child";}

    switch(test_number){
        case 1: // point joint
            rigid_body1->Frame().t=TV(0,4,2);
            joint=new POINT_JOINT<TV>();arb.joint_mesh.Add_Articulation(rigid_body0->particle_index,rigid_body1->particle_index,joint);
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,1)));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,-1,-1)));
            break;
        case 14:if(!user_last_frame) last_frame=1;
        case 2: // rigid with prismatic translation
            rigid_body1->Frame().t=TV((T)2.5,2,0);
            joint=new RIGID_JOINT<TV>();((RIGID_JOINT<TV>*)joint)->Set_Prismatic_Component_Translation(TV((T).5,0,0));
            arb.joint_mesh.Add_Articulation(rigid_body0->particle_index,rigid_body1->particle_index,joint);
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,1)));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,1,1)));
            break;
        case 3: // hinge
            rigid_body1->Frame().t=TV(2,4,0);
            joint=new ANGLE_JOINT<TV>();arb.joint_mesh.Add_Articulation(rigid_body0->particle_index,rigid_body1->particle_index,joint);
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,1),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1,1),ROTATION<TV>((T)pi/2,TV(0,0,1))*ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            break;
        case 29:if(!user_last_frame) last_frame=240;
        case 4: // twist
            rigid_body1->Frame().t=TV((T)2.1,2,0);
            joint=new ANGLE_JOINT<TV>();arb.joint_mesh.Add_Articulation(rigid_body0->particle_index,rigid_body1->particle_index,joint);
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,0,0)));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-(T)1.1,0,0)));
            if(test_number==4){
                rigid_body1->Angular_Momentum()=TV(15,0,0);
                rigid_body1->Update_Angular_Velocity();}
            else{
                rigid_body1->is_static=true;
                arb.Use_PD_Actuators();
                JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
                joint_function->Set_k_p(1000);
                joint_function->Set_Target_Angular_Velocity(TV(15,0,0));}
            break;
        case 20:
        case 5: // multiple rigid contraints
            if(!user_last_frame) last_frame=240;
            for(int i=0;i<8;i++){
                RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("subdivided_box",1,(T).5);
                rigid_body.Frame().t=TV((T)2*i,(T)2*i,(T)2*i);
                rigid_body.Set_Coefficient_Of_Restitution((T).9);
                if(i>0){
                    joint=new RIGID_JOINT<TV>();arb.joint_mesh.Add_Articulation(rigid_body.particle_index-1,rigid_body.particle_index,joint);
                    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,1)));
                    joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1,-1)));}}
            break;
        case 6: // closed loop
            if(!user_last_frame) last_frame=240;
            tests.Make_Lathe_Chain(FRAME<TV>(TV(0,10,0)));
            break;
        case 7: // pd control
            PD_Curl(TV(),ROTATION<TV>(),50);
            arb.Use_PD_Actuators();
            arb.global_post_stabilization=true;
            break;
        case 8:{ // cluster with rigid block of same size
            if(!user_last_frame) last_frame=132;
            Large_Cluster_Cube(FRAME<TV>(TV(0,2,0),rotation),(T)1,(T)0);
            RIGID_BODY<TV>*rigid_body=&tests.Add_Rigid_Body("subdivided_box",2,(T)0);
            rigid_body->Frame()=FRAME<TV>(TV(10,2,0),rotation);
            break;}
        case 9:{ // cluster break
            break;}
        case 10:{ // prismatic twist joint
            rigid_body0->Frame().t=TV(2,4,0);
            rigid_body1->Frame().t=TV(0,2,0);
            PRISMATIC_TWIST_JOINT<TV>* joint=new PRISMATIC_TWIST_JOINT<TV>();
            arb.joint_mesh.Add_Articulation(rigid_body0->particle_index,rigid_body1->particle_index,joint);
            joint->Set_Prismatic_Constraints(VECTOR<bool,3>(true,true,true),TV(0,-2,0),TV(0,2,0));
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,0,0)));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,0,0)));
            break;}
        case 11:{ // constrained twist
            rigid_body1->Frame().t=TV((T)2.1,2,0);
            rigid_body1->Angular_Momentum()=TV(3,0,0);
            rigid_body1->Update_Angular_Velocity();
            ANGLE_JOINT<TV>* joint=new ANGLE_JOINT<TV>();arb.joint_mesh.Add_Articulation(rigid_body0->particle_index,rigid_body1->particle_index,joint);
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,0,0)));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-(T)1.1,0,0)));
            joint->Set_Angle_Constraints(true,-(T)pi/4,(T)pi/4);
            break;}
        case 12: Heavy_Bottom_Link_Test(); break;
        case 13: // universal joint (no twist or translation components)
            rigid_body1->Frame().t=TV(0,2,3);
            rigid_body0->Angular_Momentum()=TV(5,5,5);
            rigid_body0->Update_Angular_Velocity();
            rigid_body1->Update_Angular_Velocity();
            joint=new POINT_JOINT<TV>();arb.joint_mesh.Add_Articulation(rigid_body0->particle_index,rigid_body1->particle_index,joint);
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(0,0,(T)1.5),ROTATION<TV>((T)pi/2,TV(0,1,0))));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,-(T)1.5),ROTATION<TV>((T)pi/2,TV(0,1,0))));
            static_cast<POINT_JOINT<TV>*>(joint)->Use_Twist_Constraint(0,0);
            LOG::cout<<"need both: "<<(joint->Has_Prismatic_Constraint() && joint->Has_Angular_Constraint())<<std::endl;
            break;
        case 15: // rigid as 3 point joints
            rigid_body1->Frame().t=TV((T)2.5,2,0);
            joint=new POINT_JOINT<TV>();arb.joint_mesh.Add_Articulation(rigid_body0->particle_index,rigid_body1->particle_index,joint);
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T)1.25,1,1)));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-(T)1.25,1,1)));
            joint=new POINT_JOINT<TV>();arb.joint_mesh.Add_Articulation(rigid_body0->particle_index,rigid_body1->particle_index,joint);
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T)1.25,-1,1)));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-(T)1.25,-1,1)));
            joint=new POINT_JOINT<TV>();arb.joint_mesh.Add_Articulation(rigid_body0->particle_index,rigid_body1->particle_index,joint);
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T)1.25,1,-1)));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-(T)1.25,1,-1)));
            break;
        case 16:{ // twist as 2 point joints
            T separation=(T)2.1,distance=0;
            rigid_body1->Frame().t=TV(separation,2,0);
            rigid_body1->Angular_Momentum()=TV(15,0,0);
            rigid_body1->Update_Angular_Velocity();
            joint=new POINT_JOINT<TV>();arb.joint_mesh.Add_Articulation(rigid_body0->particle_index,rigid_body1->particle_index,joint);
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-distance,0,0)));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-(distance+separation),0,0)));
            joint=new POINT_JOINT<TV>();arb.joint_mesh.Add_Articulation(rigid_body0->particle_index,rigid_body1->particle_index,joint);
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(distance+separation,0,0)));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(distance,0,0)));
            break;}
        case 17: // muscle control
            Muscle_Curl(TV(),ROTATION<TV>((T)pi/6,TV(0,0,1)));
            break;
        case 18:{ // muscle hinge
            arb.muscle_list->muscle_force_curve.Initialize(data_directory);
            rigid_body1->Frame().t=TV(2,4,0);
            joint=new ANGLE_JOINT<TV>();arb.joint_mesh.Add_Articulation(rigid_body0->particle_index,rigid_body1->particle_index,joint);
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,1),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1,1),ROTATION<TV>((T)pi/2,TV(0,0,1))*ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            arb.Create_Joint_Function(joint->id_number)->muscle_control=true;
            peak_force=1000;
            MUSCLE<TV>* muscle1=Add_Basic_Muscle("muscle",*rigid_body0,TV(1,-1,0),*rigid_body1,TV(1,-1,0));
            muscle1->Set_Optimal_Length(sqrt((T)8.0));
            MUSCLE<TV>* muscle2=Add_Basic_Muscle("muscle",*rigid_body0,TV(-1,1,0),*rigid_body1,TV(-1,1,0));
            muscle2->Set_Optimal_Length(sqrt((T)8.0));
            Initialize_Muscle_Segments();
            arb.Use_Muscle_Actuators();
            break;}
        case 19:{ // nested clusters
#if 0
            solids_parameters.fracture=true;
            FRACTURE_EVOLUTION_CLUSTER_3D<TV>* fracture_evolution=new FRACTURE_EVOLUTION_CLUSTER_3D<TV>(solids_parameters);
            solids_parameters.Set_Fracture_Evolution(fracture_evolution);
            Nested_Clusters_Test();
#endif
            break;}
        case 21:
            Sphere_Mesh();break;
        case 22:
            Overconstrained_Joint();break;
        case 23:{ // spring with clusters
            Spring_Cluster_Test(FRAME<TV>(TV(0,3,0)),false);
            break;}
        case 24:{// spring with clusters & fracture
#if 0
            solids_parameters.fracture=true;
            FRACTURE_EVOLUTION_CLUSTER_3D<TV>* fracture_evolution=new FRACTURE_EVOLUTION_CLUSTER_3D<TV>(solids_parameters);
            solids_parameters.Set_Fracture_Evolution(fracture_evolution);
            Spring_Cluster_Test(FRAME<TV>(TV(0,3,0)),true);
#endif
            break;}
        case 25:
        case 26:{
            if(!user_last_frame) last_frame=480;
            RIGID_BODY<TV>& box=tests.Add_Rigid_Body("subdivided_box",1,(T)0);
            box.Frame().t.y=20;box.is_static=true;
            RIGID_BODY<TV>& sphere=tests.Add_Rigid_Body("sphere",(T)1,(T)0);
            sphere.Frame().t=TV(0,10,10);
            TV offset=box.Frame().t-sphere.Frame().t;
            JOINT<TV>* joint=0;
            if(test_number==25) joint=new POINT_JOINT<TV>();else joint=new ANGLE_JOINT<TV>();
            arb.joint_mesh.Add_Articulation(box.particle_index,sphere.particle_index,joint);
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>());
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(offset));
            T length=offset.Normalize(),inertia_adjustment=1+sphere.Inertia_Tensor()(1,1)/(sphere.Mass()*sqr(length));
            nonlinear_pendulum.g_over_l=(T)9.8/(length*inertia_adjustment);
            nonlinear_pendulum.y_not=VECTOR<T,2>(-asin(offset.z),0);
            break;}
        case 27: Normal_Joint_Test();break;
        case 28: Kinematic_Angle_Joint_Test();break;
        case 30:{
            if(!user_last_frame) last_frame=240;
            RIGID_BODY<TV>* box=&tests.Add_Rigid_Body("subdivided_box",5,(T)0);
            box->Frame().t=TV(0,3,0);
            box=&tests.Add_Rigid_Body("subdivided_box",5,(T)0);
            box->Frame().t=TV(0,24,0);
            ARRAY<RIGID_BODY<TV>*> connections;
            connections.Append(&tests.Add_Rigid_Body("subdivided_box",1,(T)0));
            connections.Last()->Frame().t=TV(-3,13.5,0);
            connections.Append(&tests.Add_Rigid_Body("subdivided_box",1,(T)0));
            connections.Last()->Frame().t=TV(3,13.5,0);
            connections.Append(&tests.Add_Rigid_Body("subdivided_box",1,(T)0));
            connections.Last()->Frame().t=TV(0,13.5,-3);
            connections.Append(&tests.Add_Rigid_Body("subdivided_box",1,(T)0));
            connections.Last()->Frame().t=TV(0,13.5,3);
            for(int i=0;i<2;i++){
                T y=i==0?10.5:16.5;
                for(int j=0;j<8;j++){
                    T x=j<3?-3:j>4?3:0;T z=(j==2||j==4||j==7)?-3:(j==1||j==6)?0:3;
                    RIGID_BODY<TV>* box=&tests.Add_Rigid_Body("subdivided_box",1,(T)0);
                    box->Frame().t=TV(x,y,z);
                    for(int k=0;k<connections.m;k++){
                        JOINT<TV>* joint=new POINT_JOINT<TV>();
                        arb.joint_mesh.Add_Articulation(box->particle_index,connections(k)->particle_index,joint);
                        FRAME<TV> joint_frame;joint_frame.t=(box->Frame().t+connections(k)->Frame().t)/2.;
                        joint->Set_Joint_To_Parent_Frame(connections(k)->Frame()*joint_frame.Inverse());
                        joint->Set_Joint_To_Child_Frame(box->Frame()*joint_frame.Inverse());}}}
            break;}
        default:PHYSBAM_FATAL_ERROR(LOG::sprintf("Unrecognized test number %d",test_number));}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    Get_Initial_Data();

    if(test_number!=27 && test_number!=28) tests.Add_Ground((T).5,-2,1);

    // add forces
    if(test_number!=28) tests.Add_Gravity();

    // correct number nodes
    for(int i=0;i<deformable_body_collection.structures.m;i++) deformable_body_collection.structures(i)->Update_Number_Nodes();

    // binding_list.Distribute_Mass_To_Parents(particles.mass.array);
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();
}
//#####################################################################
// Function Sphere_Mesh
//#####################################################################
void Sphere_Mesh()
{
    T structure_scale=20;
    TRIANGULATED_SURFACE<T>* pattern=TESSELLATION::Generate_Triangles(SPHERE<TV>(),parameter);
    SEGMENT_MESH& segment_mesh=pattern->mesh.Get_Segment_Mesh();

    T edge_length=10;
    for(int i=0;i<segment_mesh.elements.m;i++)
        edge_length=min(edge_length,(pattern->particles.X(segment_mesh.elements(i).x)-pattern->particles.X(segment_mesh.elements(i).y)).Magnitude()*structure_scale);
    T sphere_size=edge_length*(T).18;
    T cylinder_size=edge_length*(T).15;
    T offset_scale=1.5;
    
    for(int i=0;i<pattern->particles.Size();i++){
        RIGID_BODY<TV>& small_sphere=tests.Add_Rigid_Body("sphere",sphere_size,(T)0);
        small_sphere.Frame().t=structure_scale*pattern->particles.X(i) + offset_scale*structure_scale*TV(0,1,0);}

    for(int i=0;i<segment_mesh.elements.m;i++){
        RIGID_BODY<TV>& cylinder=tests.Add_Rigid_Body("cyllink",cylinder_size,(T)0);
        VECTOR<TV,2> endpoints(pattern->particles.X.Subset(segment_mesh.elements(i)));
        TV average=((T).5)*(endpoints.x + endpoints.y);
        cylinder.Frame().t=structure_scale*average + offset_scale*structure_scale*TV(0,1,0);
        cylinder.Frame().r=(ROTATION<VECTOR<T,3> >)MATRIX<T,3>::Rotation_Matrix(TV(0,1,0),(endpoints(1)-endpoints(0)).Normalized());
        RIGID_BODY<TV>& body_x=solid_body_collection.rigid_body_collection.Rigid_Body(segment_mesh.elements(i).x);
        RIGID_BODY<TV>& body_y=solid_body_collection.rigid_body_collection.Rigid_Body(segment_mesh.elements(i).y);
        tests.Connect_With_Point_Joint(body_x,cylinder,body_x.Frame().t);
        tests.Connect_With_Point_Joint(body_y,cylinder,body_y.Frame().t);}

    delete pattern;
}
//#####################################################################
// Function Sphere_Mesh
//#####################################################################
void Overconstrained_Joint()
{
    RIGID_BODY<TV>& left=tests.Add_Rigid_Body("box",1,(T).5);
    RIGID_BODY<TV>& right=tests.Add_Rigid_Body("box",1,(T).5);
    RIGID_BODY<TV>& sphere=tests.Add_Rigid_Body("sphere",1,(T).5);
    tests.Connect_With_Point_Joint(left,sphere,TV());
    tests.Connect_With_Point_Joint(right,sphere,TV());
    left.Frame().t=TV(-3,1,0);
    right.Frame().t=TV(3,1,0);
    sphere.Frame().t=TV(0,1,0);
    left.is_static=true;
    right.is_static=true;
}
//#####################################################################
// Function PD Curl
//#####################################################################
void PD_Curl(TV shift,ROTATION<TV> orient,const T k_p)
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    JOINT<TV>* joint=0;RIGID_BODY<TV> *parent_body=0,*child_body=0;
    T cheight=(T)0;

    // Create first body
    parent_body=&tests.Add_Rigid_Body("miniplank25wide2",1,(T).5);
    parent_body->Frame().t=shift;
    parent_body->Frame().r=orient;
    parent_body->name="parent";
    parent_body->Set_Mass(5);
    parent_body->is_static=true;

    // Add children and joints
    T desired_x=((T)pi*2)/((T)(parameter+1));
    for(int i=0;i<parameter;i++){
        cheight+=(T)1.25;
        child_body=&tests.Add_Rigid_Body("miniplank25wide2",(T).8,(T).5);
        child_body->Frame().t=orient.Rotate(TV(cheight,0,0))+shift;
        child_body->Frame().r=orient;
        child_body->Set_Coefficient_Of_Restitution((T)0.5);
        child_body->name=LOG::sprintf("child_%d",i);

        joint=new POINT_JOINT<TV>();arb.joint_mesh.Add_Articulation(child_body->particle_index-1,child_body->particle_index,joint);
        ROTATION<TV> desired_rotation=ROTATION<TV>(desired_x,TV(1,0,0));
        JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
        joint_function->Set_k_p(k_p);joint_function->Set_Target_Angle(desired_rotation);
        joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T).625,0,0),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
        joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-(T).625,0,0),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));

        parent_body=child_body;}
}
//#####################################################################
// Function Muscle_Curl
//#####################################################################
void Muscle_Curl(TV shift,ROTATION<TV> orient)
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    arb.muscle_list->muscle_force_curve.Initialize(data_directory);

    // Create first body
    RIGID_BODY<TV>* parent=&tests.Add_Rigid_Body("miniplank25wide2",1,(T).5);
    parent->Frame().t=shift;
    parent->Frame().r=orient;
    parent->name="parent";
    parent->Set_Mass(5);
    parent->is_static=true;

    int njoints=100;

    TV radius((T).625,0,0);
    ROTATION<TV> twist((T)pi/6,TV(0,1,0));

    // Add children and joints
    for(int i=0;i<njoints;i++){
        shift+=orient.Rotate(radius);
        orient=twist*orient;
        shift+=orient.Rotate(radius);

        RIGID_BODY<TV>* child=&tests.Add_Rigid_Body("miniplank25wide2",1,(T).5);
        child->Frame().r=orient;
        child->Frame().t=shift;
        child->Set_Coefficient_Of_Restitution(0.5);
        child->name=LOG::sprintf("child_%d",i);

        POINT_JOINT<TV>& joint=*new POINT_JOINT<TV>();arb.joint_mesh.Add_Articulation(child->particle_index-1,child->particle_index,&joint);
        joint.Set_Joint_To_Parent_Frame(FRAME<TV>(radius,ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
        joint.Set_Joint_To_Child_Frame(child->Frame().Inverse()*parent->Frame()*joint.F_pj()*FRAME<TV>());
        arb.Create_Joint_Function(joint.id_number)->muscle_control=true;
        Add_Basic_Muscle("muscle",*parent,TV(),*child,TV());

        parent=child;}

    Initialize_Muscle_Segments();
    arb.Use_Muscle_Actuators();
}
//#####################################################################
// Function Add_Basic_Muscle
//#####################################################################
MUSCLE<TV>* Add_Basic_Muscle(const std::string& name,RIGID_BODY<TV>& origin_body,const TV& origin,RIGID_BODY<TV>& insertion_body,const TV& insertion)
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    MUSCLE<TV>* muscle=new MUSCLE<TV>(arb.muscle_list->muscle_force_curve);
    muscle->name=name;
    muscle->Set_Attachment_Point_1(new ATTACHMENT_POINT<TV>(origin_body,origin));
    muscle->Set_Attachment_Point_2(new ATTACHMENT_POINT<TV>(insertion_body,insertion));
    T total_length=muscle->Total_Length();
    muscle->Set_Optimal_Length((T).8*total_length);
    muscle->Set_Tendon_Slack_Length((T).2*total_length);
    muscle->Set_Peak_Force(peak_force);
    muscle->Set_Max_Shortening_Velocity(1);
    arb.muscle_list->Add_Muscle(muscle);
    return muscle;
}
//#####################################################################
// Function Initialize Muscle_Segments
//#####################################################################
void Initialize_Muscle_Segments()
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    for(int i=0;i<arb.muscle_list->muscles.m;i++) arb.muscle_list->muscles(i)->Initialize();
}
//#####################################################################
// Function Large_Cluster_Cube
//#####################################################################
int Large_Cluster_Cube(FRAME<TV>shift_frame,T scale,const T friction)
{
    // CLUSTER 1
    ARRAY<RIGID_BODY<TV>*>& bodies=*new ARRAY<RIGID_BODY<TV>*>(8);
    int count=0;
    for(int i=0;i<8;i++){
//        bodies(i)->name=LOG::sprintf("child::%d",bodies(i)));}
        bodies(i)=&tests.Add_Rigid_Body("subdivided_box",1,friction);
        bodies(i)->name=LOG::sprintf("child::%d",bodies(i));}
    for(int i=-1;i<=1;i+=2) for(int j=-1;j<=1;j+=2) for(int k=-1;k<=1;k+=2) bodies(++count)->Frame()=shift_frame*FRAME<TV>(TV((T)k,(T)j,(T)i));
    tests.Add_Gravity();
//    cluster_id=rigid_body_particles.Add_Cluster_Body(&bodies);
    solid_body_collection.rigid_body_collection.Rigid_Body(cluster_id).name="combo_square";
    return cluster_id;
}
//#####################################################################
// Function Nested_Clusters_Test
//#####################################################################
int Nested_Clusters_Test()
{
    PHYSBAM_FATAL_ERROR("RIGID_BODY_CLUSTER_3D has been replaced by RIGID_BODY_CLUSTER_BINDING");
    return int();
}
//#####################################################################
// Function Spring_Cluster_Test
//#####################################################################
int Spring_Cluster_Test(FRAME<TV> shift_frame,bool fracture)
{
    PHYSBAM_FATAL_ERROR("RIGID_BODY_CLUSTER_3D has been replaced by RIGID_BODY_CLUSTER_BINDING");
    return int();
}
//#####################################################################
// Draped Block Rope
//#####################################################################
void Heavy_Bottom_Link_Test()
{
    if(!user_last_frame) last_frame=240;
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    JOINT<TV> *joint;
    RIGID_BODY<TV>* rigid_body[3];

    arb.Set_Contact_Level_Iterations(50);
    arb.Set_Shock_Propagation_Level_Iterations(50);

    for(int i=0;i<3;i++){
        rigid_body[i]=&tests.Add_Rigid_Body("subdivided_box",1,(T).5);
        rigid_body[i]->Frame().r=ROTATION<TV>(-(T)pi/4,TV(0,0,1));
        rigid_body[i]->Frame().t=rigid_body[i]->Frame().r.Rotate(TV(0,(T)-3*i,0))+TV(0,5,0);}

    for(int i=1;i<3;i++){
        joint=new POINT_JOINT<TV>();arb.joint_mesh.Add_Articulation(rigid_body[i-1]->particle_index,rigid_body[i]->particle_index,joint);
        joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(0,-(T)1.5,0)));
        joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,(T)1.5,0)));}

    rigid_body[0]->is_static=true;
    rigid_body[2]->Mass()*=10000;
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) override
{
    if(test_number==14) Test_System_Prestabilization(JOINT_ID(0));
    else if(test_number==20) Test_System_Poststabilization();
//    else if((test_number==8 && frame==35) || (test_number==23 && frame==10)){
//        rigid_body_particles.Remove_Cluster_Body(cluster_id);
//    if(test_number==24 && frame==10){
//        ((RIGID_BODY_CLUSTER_3D<T>*)rigid_body_list(cluster_id))->Initialize_Allowable_Strain((T)0);}
}
//#####################################################################
// Function Test_System_Prestabilization_Print
//#####################################################################
template<class T_VECTOR> void
Test_System_Prestabilization_Print(const T_VECTOR& j)
{
    for(int i=0;i<j.Size();i++) LOG::cout<<(i>1?" ":"")<<j(i);
}
//#####################################################################
// Function Test_System_Prestabilization_Print
//#####################################################################
template<class T_VECTOR,class T_MATRIX> void
Test_System_Prestabilization_Print(const T_VECTOR& f1,const T_VECTOR& f2,const T_MATRIX& m1,const T_MATRIX& m2)
{
    LOG::cout<<"Constraint Error 1:  ";Test_System_Prestabilization_Print(f1);LOG::cout<<std::endl;
    LOG::cout<<"Constraint Error 2:  ";Test_System_Prestabilization_Print(f2);LOG::cout<<std::endl;
    LOG::cout<<"Difference: ";Test_System_Prestabilization_Print(f1-f2);LOG::cout<<std::endl;
    LOG::cout<<"Jacobian 1:\n"<<m1<<std::endl;
    LOG::cout<<"Jacobian 2:\n"<<m2<<std::endl;
    LOG::cout<<"Difference:\n"<<(m1-m2)<<std::endl;
}
//####################################################################################
// Function Post_Stabilization_Constraint_Matrix
//####################################################################################
template<class T> MATRIX_MXN<T>
Post_Stabilization_Constraint_Matrix(const JOINT_ID joint_id)
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    JOINT<TV>& joint=*arb.joint_mesh(joint_id);
    RIGID_BODY<TV>* rigid_bodies[]={arb.Parent(joint_id),arb.Child(joint_id)};

    // get linear and angular constrained axes
    MATRIX_MXN<T> prismatic_constraints,angular_constraints;
    joint.Prismatic_Constraint_Matrix(rigid_bodies[0]->Frame(),prismatic_constraints);
    joint.Angular_Constraint_Matrix(rigid_bodies[0]->Frame(),angular_constraints);
    int p=prismatic_constraints.Columns(),a=angular_constraints.Columns();

    int d=TV::m,s=T_SPIN::m;
    MATRIX_MXN<T> R_D[2];
    TV location=joint.Location(*rigid_bodies[0],*rigid_bodies[1]);
    for(int i=0;i<2;i++){
        R_D[i].Resize(d+s,p+a);
        R_D[i].Set_Submatrix(0,0,prismatic_constraints);
        R_D[i].Set_Submatrix(d,0,prismatic_constraints.Cross_Product_Matrix_Times(location-rigid_bodies[i]->Frame().t));
        R_D[i].Set_Submatrix(d,p,angular_constraints);}

    // constraint matrix
    MATRIX_MXN<T> constraint(2*(d+s),p+a);
    constraint.Set_Submatrix(0,0,R_D[0]);constraint.Set_Submatrix(d+s,0,-R_D[1]);

    return constraint;
}
//#####################################################################
// Function Test_System_Poststabilization
//#####################################################################
void Test_System_Poststabilization()
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    arb.Initialize_Poststabilization_Projection();

    for(int i=0;i<arb.joint_mesh.Num_Joints();i++){JOINT_ID joint_id=arb.joint_mesh.Joints(i)->id_number;
        // check that the projection matrices are projections
        LOG::cout<<"--------------- JOINT = "<<joint_id<<"---------------"<<std::endl;
        MATRIX_MXN<T> P=-arb.lambda_to_delta_v(joint_id)*arb.v_to_lambda(joint_id)+1;

        LOG::cout<<"P=\n"<<P<<std::endl;

        LOG::cout<<"P.Frobenius_Norm="<<sqrt(P.Frobenius_Norm_Squared())<<std::endl;
        LOG::cout<<"P.Infinity_Norm="<<P.Infinity_Norm()<<std::endl;
        MATRIX_MXN<T> difference=P*P-P;
        LOG::cout<<"difference.Frobenius_Norm="<<sqrt(difference.Frobenius_Norm_Squared())<<std::endl;
        LOG::cout<<"difference.Infinity_Norm="<<difference.Infinity_Norm()<<std::endl;
        LOG::cout<<"P*P-P=\n"<<difference<<std::endl;

        // check that a projected vector satisfies the constraint, i.e., C^T P v = 0 forall v
        MATRIX_MXN<T> C=Post_Stabilization_Constraint_Matrix<T>(joint_id);
        MATRIX_MXN<T> CT_P=C.Transpose_Times(P);
        LOG::cout<<"C.Infinity_Norm="<<C.Infinity_Norm()<<std::endl;
        LOG::cout<<"CT_P.Infinity_Norm="<<CT_P.Infinity_Norm()<<std::endl;
        LOG::cout<<"CT_P=\n"<<CT_P<<std::endl;}
}
//#####################################################################
// Function Test_System_Prestabilization
//#####################################################################
void Test_System_Prestabilization(const JOINT_ID joint_id)
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    RIGID_BODY<TV>* parent=const_cast<RIGID_BODY<TV>*>(arb.Parent(joint_id)),*child=const_cast<RIGID_BODY<TV>*>(arb.Child(joint_id));
    T dt=(T).1;
    T epsilon_scale=1;
    T epsilon=(T)1e-3;
    T impulse_magnitude=1;
    TV location;

    parent->is_static=true;
    RIGID_BODY<TV>::Apply_Impulse(*parent,*child,location,TV(-(T)3.4,(T)5.8,-(T).9),TV((T).73,-(T).43,(T)4.9));
    parent->is_static=false;
    child->is_static=true;
    RIGID_BODY<TV>::Apply_Impulse(*parent,*child,location,TV((T).25,(T)9.2,-(T)1.0),TV((T)3.7,(T)0.2,-(T)8.2));
    child->is_static=false;

    LOG::cout<<"Linear test"<<std::endl;
    for(int i=0;i<3;i++){
        LOG::cout<<"--------------------------------------------------------------------------------"<<std::endl;
        TV j;j(i)=impulse_magnitude;
        LINEAR_CONSTRAINT_FUNCTION<TV> f_error(arb,joint_id,dt,epsilon_scale,location);
        typename LINEAR_CONSTRAINT_FUNCTION<TV>::T_CONSTRAINT_ERROR f_error_result=f_error.F(j);
        MATRIX<T,3> jacobian=f_error.Jacobian(j);
        RIGID_BODY<TV>::Apply_Impulse(*parent,*child,location,j/dt);
        LINEAR_CONSTRAINT_FUNCTION<TV> f_error_2(arb,joint_id,dt,epsilon_scale,location);
        typename LINEAR_CONSTRAINT_FUNCTION<TV>::T_CONSTRAINT_ERROR f_error_2_result=f_error_2.F(TV());
        MATRIX<T,3> jacobian_2=f_error_2.Jacobian(TV());
        RIGID_BODY<TV>::Apply_Impulse(*parent,*child,location,-j/dt);
        Test_System_Prestabilization_Print(f_error_result,f_error_2_result,jacobian,jacobian_2);}

    LOG::cout<<"Angular test"<<std::endl;
    for(int i=0;i<3;i++){
        LOG::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;
        TV j;j(i)=impulse_magnitude;
        ANGULAR_CONSTRAINT_FUNCTION<TV> f_error(arb,joint_id,dt,epsilon_scale);
        typename ANGULAR_CONSTRAINT_FUNCTION<TV>::T_CONSTRAINT_ERROR f_error_result=f_error.F(j);
        MATRIX<T,3> jacobian=f_error.Jacobian(j);
        RIGID_BODY<TV>::Apply_Impulse(*parent,*child,TV(),TV(),j/dt);
        ANGULAR_CONSTRAINT_FUNCTION<TV> f_error_2(arb,joint_id,dt,epsilon_scale);
        typename ANGULAR_CONSTRAINT_FUNCTION<TV>::T_CONSTRAINT_ERROR f_error_2_result=f_error_2.F(TV());
        MATRIX<T,3> jacobian_2=f_error_2.Jacobian(TV());
        RIGID_BODY<TV>::Apply_Impulse(*parent,*child,TV(),TV(),-j/dt);
        if(TV::Dot_Product(f_error_result,f_error_2_result)<0){
            LOG::cout<<"Flipping sign"<<std::endl;
            f_error_2_result=-f_error_2_result;jacobian_2=-jacobian_2;} // negating both f and j is ok
        Test_System_Prestabilization_Print(f_error_result,f_error_2_result,jacobian,jacobian_2);}

    LOG::cout<<"Combined test"<<std::endl;
    for(int i=0;i<6;i++){
        LOG::cout<<"-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"<<std::endl;
        VECTOR<T,6> j_combined;j_combined(i)=impulse_magnitude;
        TV j,j_tau;
        j_combined.Extract(j,j_tau);
        LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<TV> f_error(arb,joint_id,dt,epsilon_scale,location);
        typename LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<TV>::T_CONSTRAINT_ERROR f_error_result=f_error.F(j_combined);
        MATRIX<T,6> jacobian=f_error.Jacobian(j_combined);
        RIGID_BODY<TV>::Apply_Impulse(*parent,*child,location,j/dt,j_tau/dt);
        LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<TV> f_error_2(arb,joint_id,dt,epsilon_scale,location);
        typename LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<TV>::T_CONSTRAINT_ERROR f_error_2_result=f_error_2.F(VECTOR<T,6>());
        MATRIX<T,6> jacobian_2=f_error_2.Jacobian(VECTOR<T,6>());
        RIGID_BODY<TV>::Apply_Impulse(*parent,*child,location,-j/dt,-j_tau/dt);
        TV f_angular(f_error_result.Array_View(3,3)),f_angular_2(f_error_2_result.Array_View(3,3));
        LOG::cout<<"Testing component "<<i<<std::endl;
        if(TV::Dot_Product(f_angular,f_angular_2)<0){
            LOG::cout<<"Flipping sign"<<std::endl;
            f_error_2_result.Array_View(3,3)=-f_angular_2;
            for(int i=3;i<6;i++) for(int j=0;j<6;j++) jacobian_2(i,j)=-jacobian_2(i,j);}
        Test_System_Prestabilization_Print(f_error_result,f_error_2_result,jacobian,jacobian_2);}

    LOG::cout<<"################################################################################"<<std::endl;
    LOG::cout<<"Check Jacobians"<<std::endl;

    TV jn(-(T)3.4,(T)5.8,-(T).9),j_tau((T).73,-(T).43,(T)4.9);

    {LINEAR_CONSTRAINT_FUNCTION<TV> f_error(arb,joint_id,dt,epsilon_scale,location);
    MATRIX<T,3> jacobian(f_error.Jacobian(jn)),jacobian_2;
    for(int i=0;i<3;i++){TV j;j(i)=epsilon;jacobian_2.Column(i)=(f_error.F(jn+j)-f_error.F(jn-j))/(2*epsilon);}
    LOG::cout<<"Computed Linear Jacobian:\n"<<jacobian<<"Approximated Linear Jacobian:\n"<<jacobian_2;
    LOG::cout<<"Difference:\n"<<(jacobian-jacobian_2)<<std::endl;}

    {ANGULAR_CONSTRAINT_FUNCTION<TV> f_error(arb,joint_id,dt,epsilon_scale);
    MATRIX<T,3> jacobian(f_error.Jacobian(j_tau)),jacobian_2;
    for(int i=0;i<3;i++){TV j;j(i)=epsilon;jacobian_2.Column(i)=(f_error.F(j_tau+j)-f_error.F(j_tau-j))/(2*epsilon);}
    LOG::cout<<"Computed Angular Jacobian:\n"<<jacobian<<"Approximated Angular Jacobian:\n"<<jacobian_2;
    LOG::cout<<"Difference:\n"<<(jacobian-jacobian_2)<<std::endl;}

    {LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<TV> f_error(arb,joint_id,dt,epsilon_scale,location);
    VECTOR<T,6> j_combined;
    j_combined.Combine(jn,j_tau);
    MATRIX<T,6> jacobian(f_error.Jacobian(j_combined)),jacobian_2;
    for(int i=0;i<6;i++){VECTOR<T,6> j;j(i)=epsilon;jacobian_2.Set_Column(i,(f_error.F(j_combined+j)-f_error.F(j_combined-j))/(2*epsilon));}
    LOG::cout<<"Computed Combined Jacobian:\n"<<jacobian<<"Approximated Combined Jacobian:\n"<<jacobian_2;
    LOG::cout<<"Difference:\n"<<(jacobian-jacobian_2)<<std::endl;}

    // post-stabilization
    POINT_JOINT<TV>* point_joint=new POINT_JOINT<TV>();arb.joint_mesh.Add_Articulation(parent->particle_index,child->particle_index,point_joint);
    point_joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,1)));
    point_joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,-1,-1)));
    TWIST<TV> delta_relative_twist;arb.Delta_Relative_Twist(point_joint->id_number,false,location,delta_relative_twist);
    LOG::cout<<"Point joint: Delta relative twist before post-stabilization="<<delta_relative_twist<<std::endl;
    arb.Apply_Poststabilization_To_Joint(point_joint->id_number,false);
    arb.Delta_Relative_Twist(point_joint->id_number,false,location,delta_relative_twist);
    LOG::cout<<"             Delta relative twist after post-stabilization="<<delta_relative_twist<<std::endl;

    parent->is_static=true;
    RIGID_BODY<TV>::Apply_Impulse(*parent,*child,location,TV(-(T)3.4,(T)5.8,-(T).9),TV((T).73,-(T).43,(T)4.9));
    parent->is_static=false;
    child->is_static=true;
    RIGID_BODY<TV>::Apply_Impulse(*parent,*child,location,TV((T).25,(T)9.2,-(T)1.0),TV((T)3.7,(T)0.2,-(T)8.2));
    child->is_static=false;

    arb.Delta_Relative_Twist(joint_id,false,location,delta_relative_twist);
    LOG::cout<<"Rigid joint: Delta relative twist before post-stabilization="<<delta_relative_twist<<std::endl;
    arb.Apply_Poststabilization_To_Joint(joint_id,false);
    arb.Delta_Relative_Twist(joint_id,false,location,delta_relative_twist);
    LOG::cout<<"             Delta relative twist after post-stabilization="<<delta_relative_twist<<std::endl;

    exit(1);
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) override
{
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particles;
    if(test_number==25 || test_number==26){
        ODE_SOLVER<T,VECTOR<T,2> > rk;
        nonlinear_pendulum.y_not=rk.Runge_Kutta_4(nonlinear_pendulum,nonlinear_pendulum.y_not,0,(T)1.0/frame_rate,50);
        TV offset=rigid_body_particles.frame(0).t-rigid_body_particles.frame(1).t;
        offset.Normalize();T theta=-asin(offset.z);
        LOG::cout<<"ANALYTIC TEST:   theoretical: "<<nonlinear_pendulum.y_not.x<<"      computed: "<<theta<<std::endl;}
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) override
{
    if(test_number==28 && id==int(0))
        frame=curve.Value(time);
}
//#####################################################################
// Function Normal_Joint_Test
//#####################################################################
void Normal_Joint_Test()
{
    if(!user_last_frame) last_frame=240;
    solids_parameters.rigid_body_collision_parameters.contact_iterations=0;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;

    // box sliding down an inclined plane
    T mu=(T)0,ground_angle_rad=(T)20*((T)pi/180);
    FRAME<TV> frame(TV(-sin(ground_angle_rad),cos(ground_angle_rad),0),ROTATION<TV>(ground_angle_rad,TV(0,0,1)));

    RIGID_BODY<TV>& inclined_plane=tests.Add_Ground(mu,0,0);
    inclined_plane.name="inclined_plane";
    inclined_plane.Set_Coefficient_Of_Rolling_Friction(1);
    inclined_plane.Frame().r=ROTATION<TV>(ground_angle_rad,TV(0,0,1));

    // try turning pre stab off to see if it behaves physically
    RIGID_BODY<TV>& box=tests.Add_Rigid_Body("subdivided_box",1,mu);
    box.Frame()=frame;
    box.Set_Coefficient_Of_Restitution(0);
    box.name="box";
    box.Twist().angular=frame*TV((T)0,(T)2,(T)0); // with this initial velocity the sliding box is stopped by the static box

    // create the joints
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    JOINT<TV>* joint=new NORMAL_JOINT<TV>();arb.joint_mesh.Add_Articulation(inclined_plane.particle_index,box.particle_index,joint);
    FRAME<TV> J(frame*TV((T)0,-(T)1,(T)0),ROTATION<TV>::From_Rotated_Vector(TV::Axis_Vector(0),frame.r.Rotated_Axis(1)));
    joint->Set_Joint_To_Parent_Frame(inclined_plane.Frame().Inverse()*J);
    joint->Set_Joint_To_Child_Frame(box.Frame().Inverse()*J);

    // put a static box in the constrained box's way
    RIGID_BODY<TV>& box2=tests.Add_Rigid_Body("subdivided_box",1,mu);
    box2.Frame()=frame;box2.Frame().t-=frame.r.Rotated_X_Axis()*(T)10;box2.Frame().t.z-=(T)1.5;
    box2.name="static box";
    box2.is_static=true;

    // compare system without the joint
    T z_offset=(T)10;
    RIGID_BODY<TV>& box3=tests.Add_Rigid_Body("subdivided_box",1,mu);
    box3.Frame()=frame;box3.Frame().t.z+=z_offset;
    box3.Set_Coefficient_Of_Restitution(0);
    box3.name="box3";
    box3.Twist().angular=frame*TV((T)0,(T)2,(T)0); // with this initial velocity the sliding box is stopped by the static box

    // put a static box in the constrained box's way
    RIGID_BODY<TV>& box4=tests.Add_Rigid_Body("subdivided_box",1,mu);
    box4.Frame()=frame;box4.Frame().t-=frame.r.Rotated_X_Axis()*(T)10;box4.Frame().t.z+=-(T)1.5+z_offset;
    box4.name="static box";
    box4.is_static=true;
}
//#####################################################################
// Function Kinematic_Angle_Joint_Test
//#####################################################################
void Kinematic_Angle_Joint_Test()
{
    if(!user_last_frame) last_frame=240;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;

    // master
    RIGID_BODY<TV>& master=*new RIGID_BODY<TV>(rigid_body_collection,true);
    TV master_edges((T)14.0662,(T)7.12134,(T)14.1081);
    RANGE<TV> master_box(-master_edges,master_edges);
    master.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(master_box));
    master.Is_Kinematic()=true;

    // slave
    RIGID_BODY<TV>& slave=*new RIGID_BODY<TV>(rigid_body_collection,true);
    TV slave_edges((T)2.49461,(T)9.85401,(T)3.84579);
    RANGE<TV> slave_box(-slave_edges,slave_edges);
    slave.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(slave_box));
    slave.Mass()=13;
    slave.Inertia_Tensor()=DIAGONAL_MATRIX<T,3>(932,175,861);
    slave.Frame()=FRAME<TV>(TV((T)-335.915,(T)618.673,(T)9.6621),ROTATION<TV>((T)0.9439,TV((T)-0.25355,(T)-0.450801,(T)0.855857)));

    // add bodies
    rigid_body_collection.Add_Rigid_Body_And_Geometry(&master);
    rigid_body_collection.Add_Rigid_Body_And_Geometry(&slave);

    // add joint
    ANGLE_JOINT<TV>* joint=new ANGLE_JOINT<TV>();
    joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV((T)-1.46045,(T)-19.9852,(T)1.12249),ROTATION<TV>((T)0.94073,TV((T)0.683208,(T)0.223715,(T)0.69511))));
    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T)13.774,(T)14.3539,(T)12.1388),ROTATION<TV>((T)0.0458756,TV((T)-0.955585,(T)-0.287603,(T)0.0643592))));
    arb.joint_mesh.Add_Articulation(master.particle_index,slave.particle_index,joint);

    // set up animation cue
    curve.Add_Control_Point((T)0.0416667,FRAME<TV>(TV((T)-314.243,(T)589.525,(T)2.61261),ROTATION<TV>((T)1.51451,TV((T)0.117017,(T)0.122549,(T)0.98554))));
    curve.Add_Control_Point((T)0.166667,FRAME<TV>(TV((T)-297.741,(T)591.505,(T)1.71612),ROTATION<TV>((T)1.51577,TV((T)0.113364,(T)0.135305,(T)0.984297))));
    curve.Add_Control_Point((T)5.70833,FRAME<TV>(TV((T)139.226,(T)606.004,(T)-0.488798),ROTATION<TV>((T)0.598542,TV((T)0.212646,(T)-0.0477872,(T)-0.97596))));
    curve.Add_Control_Point((T)5.875,FRAME<TV>(TV((T)137.305,(T)607.236,(T)-0.67634),ROTATION<TV>((T)0.529389,TV((T)0.229795,(T)-0.0424575,(T)-0.972313))));
    curve.Add_Control_Point((T)6.91667,FRAME<TV>(TV((T)137.998,(T)606.974,(T)-1.08095),ROTATION<TV>((T)0.54646,TV((T)0.198967,(T)-0.0341637,(T)-0.97941))));
    curve.Add_Control_Point((T)7.33333,FRAME<TV>(TV((T)131.233,(T)609.832,(T)-2.01024),ROTATION<TV>((T)0.396043,TV((T)0.206426,(T)-0.0105733,(T)-0.978405))));
    curve.Add_Control_Point((T)7.91667,FRAME<TV>(TV((T)134.205,(T)609.294,(T)-3.03306),ROTATION<TV>((T)0.454507,TV((T)0.108923,(T)0.00161945,(T)-0.994049))));
    curve.Add_Control_Point((T)8.125,FRAME<TV>(TV((T)136.436,(T)608.663,(T)-3.34368),ROTATION<TV>((T)0.497564,TV((T)0.0796029,(T)0.00453764,(T)-0.996816))));
    curve.Add_Control_Point((T)8.29167,FRAME<TV>(TV((T)148.513,(T)608.131,(T)-6.66214),ROTATION<TV>((T)0.529263,TV((T)0.0750016,(T)0.00311249,(T)-0.997179))));
    curve.Add_Control_Point((T)8.5,FRAME<TV>(TV((T)147.714,(T)608.118,(T)-6.66214),ROTATION<TV>((T)0.529263,TV((T)0.0750016,(T)0.00311249,(T)-0.997179))));
    curve.Add_Control_Point((T)8.91667,FRAME<TV>(TV((T)147.714,(T)600.124,(T)-6.66213),ROTATION<TV>((T)0.529263,TV((T)0.0750016,(T)0.00311249,(T)-0.997179))));
}
//#####################################################################
};
}
#endif

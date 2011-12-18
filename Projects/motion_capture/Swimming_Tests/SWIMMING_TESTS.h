//#####################################################################
// Copyright 2008, Michael Lentine, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SWIMMING_TESTS
//#####################################################################
//   1. Floppy Person
//   2. Static Person
//   8. Floppy fish - Tamar and Craig (for reference)
//  14. Sidewinding - Tamar and Craig (for reference)
//  16. Magget - Tamar and Craig (for reference)
//#####################################################################
#ifndef __SWIMMING_TESTS__
#define __SWIMMING_TESTS__

#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/INCOMPRESSIBLE_FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_TET_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/ANGLE_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/CONSTRAINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/PRISMATIC_TWIST_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/RIGID_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Bindings/RIGID_BODY_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/ETHER_DRAG.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Motion/BODY_MOTION_SEQUENCE.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class SWIMMING_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;
    typedef VECTOR<T_input,3> TV;typedef VECTOR<int,3> TV_INT;
    typedef typename RIGID_BODY_POLICY<TV>::INERTIA_TENSOR T_INERTIA_TENSOR;
public:
    SOLIDS_STANDARD_TESTS<TV> tests;

    int parameter;
    int kinematic_id;
    INTERPOLATION_CURVE<T,FRAME<TV> > curve;

    T ground_angle_rad,mu;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;

    T aspect_ratio;
    int number_side_panels;
    ROTATION<TV> torus_rotation;
    TV torus_min,torus_max;
    T period,start_time;
    RIGID_BODY<TV>* ground;
    int number_of_joints;
    T wavelength,height_amplitude,bend_amplitude,joint_separation;
    ROTATION<TV> initial_orientation;
    ARRAY<PAIR<TETRAHEDRALIZED_VOLUME<T>*,T> > structure_clamp_time;
    ARRAY<PAIR<RIGID_BODY<TV>*,T> > rigid_body_clamp_time;
    T maximum_fall_speed;
    ARRAY<VECTOR<int,3> > surface_elements;
    ARRAY<int> surface_particles;
    int subsamples;
    RANDOM_NUMBERS<T> random_numbers;
    ARRAY<ARRAY<int> > triangle_free_particles;
    T refinement_distance;
    bool dynamic_subsampling,temporarily_disable_dynamic_subsampling;
    ARRAY<T> particle_distances;
    int old_number_particles;
    T ring_mass;
    int number_of_maggots;
    ARRAY<T> periods,phases;
    bool fish_mattress;  // low res mattress for testing
    bool plastic,use_finite_volume;
    bool strain_limit,use_implicit;
    int motion_frame_rate;
    BODY_MOTION_SEQUENCE<T> body_motion;
    ARRAY<int,int> id_to_index;
    ARRAY<int> rigid_body_ids;
    FRAME<TV> rigid_base_transform;
    FINITE_VOLUME<TV,3>* finite_volume;
    LINEAR_SPRINGS<TV>* edge_springs;
    LINEAR_TET_SPRINGS<T>* tet_springs;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::restart;
    using BASE::frame_rate;using BASE::stream_type;using BASE::solid_body_collection;using BASE::test_number;using BASE::parse_args;

    SWIMMING_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),parameter(1),rigid_body_collection(solid_body_collection.rigid_body_collection),number_of_joints(2),
        subsamples(8),refinement_distance((T).2),dynamic_subsampling(false),temporarily_disable_dynamic_subsampling(false),old_number_particles(0),ring_mass(10000),fish_mattress(false),
        plastic(false),use_finite_volume(true),strain_limit(true),use_implicit(false),finite_volume(0)
    {
    }

    ~SWIMMING_TESTS()
    {}

    // Unused callbacks
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > F,const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_String_Argument("-mtn","Specify the physbam motion file to use");
    parse_args->Add_Integer_Argument("-parameter",0,"parameter used by multiple tests to change the parameters of the test");
    parse_args->Add_Integer_Argument("-subsamples",8,"number of particle subsamples per triangle");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    fluids_parameters.simulate=false;
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.cfl=1;
    output_directory=STRING_UTILITIES::string_sprintf("Swimming_Tests/Test_%d",test_number);

    if(parse_args->Is_Value_Set("-parameter")){parameter=parse_args->Get_Integer_Value("-parameter");output_directory+=STRING_UTILITIES::string_sprintf("_param%i",parameter);}
    if(parse_args->Is_Value_Set("-subsamples")){subsamples=parse_args->Get_Integer_Value("-subsamples");output_directory+=STRING_UTILITIES::string_sprintf("_s%i",subsamples);}
    if(parse_args->Is_Value_Set("-mtn")) FILE_UTILITIES::Read_From_File(stream_type,parse_args->Get_String_Value("-mtn"),body_motion);
    else if(test_number==1) PHYSBAM_FATAL_ERROR("You must supply a mtn file to simulate a human.");

    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.triangle_collision_parameters.output_interaction_pairs=false;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
    solids_parameters.use_rigid_deformable_contact=true;
    solids_parameters.rigid_body_collision_parameters.collision_bounding_box_thickness=(T)1e-3;
    solids_parameters.triangle_collision_parameters.collisions_output_number_checked=false;
    solid_body_collection.deformable_body_collection.soft_bindings.use_gauss_seidel_for_impulse_based_collisions=true;
    solids_parameters.verbose_dt=true;
    solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;
    solids_parameters.implicit_solve_parameters.project_nullspace_frequency=1;
    frame_rate=30;motion_frame_rate=120;

    switch(test_number){
        case 1: last_frame=1200;break;
        case 2: last_frame=360;break;
        case 8: last_frame=240;break;
        case 14: last_frame=1200;break;
        case 16: last_frame=360;break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
ROTATION<TV> Upright_Orientation(const TV& x,const TV& y)
{
    TV u=(y-x).Normalized(),j(0,1,0);
    MATRIX<T,3> rotation(u,j.Projected_Orthogonal_To_Unit_Direction(u),TV::Cross_Product(u,j));
    return ROTATION<TV>(rotation);
}
TV Sidewinding_Position(const T time,const T segment)
{
    const T k=(T)two_pi/wavelength,w=(T)two_pi/period,theta=segment*k+time*w;
    T y=height_amplitude*cos(theta);
    T x=bend_amplitude*sin(theta);
    return TV(x,y,segment);
}
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    int frame_num=(int)(time*motion_frame_rate)+1;
    T alpha=time*motion_frame_rate-frame_num+1;
    RIGID_BODY<TV>* parent=&arb.rigid_body_collection.Rigid_Body(id);
    int parent_index=id_to_index(parent->particle_index);
    FRAME<TV> rigid_base_transform_parent1=rigid_base_transform,rigid_base_transform_child1=rigid_base_transform;
    FRAME<TV> rigid_base_transform_parent2=rigid_base_transform,rigid_base_transform_child2=rigid_base_transform;
    rigid_base_transform_parent1.t*=body_motion.trajectories(parent_index)(frame_num).length;
    rigid_base_transform_parent2.t*=body_motion.trajectories(parent_index)(frame_num+1).length;
    frame=FRAME<TV>::Interpolation(body_motion.trajectories(parent_index)(frame_num).targeted_transform*rigid_base_transform_parent1,
        body_motion.trajectories(parent_index)(frame_num+1).targeted_transform*rigid_base_transform_parent2,alpha);
}
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{return false;}
//#####################################################################
// Function Self_Collisions_Begin_Callback
//#####################################################################
void Self_Collisions_Begin_Callback(const T time,const int substep) PHYSBAM_OVERRIDE
{
    if(dynamic_subsampling && solids_parameters.triangle_collision_parameters.perform_self_collision) Update_Subsamples();
}
//#####################################################################
// Function Preprocess_Substep
//#####################################################################
void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
{
    TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>(1);
    if(plastic){
        /*if(use_finite_volume){
            finite_volume->strain_measure.Initialize_Dm_Inverse(tetrahedralized_volume->particles.X);
            delete finite_volume->dPi_dFe;finite_volume->dPi_dFe=0;
            delete finite_volume->dP_dFe;finite_volume->dP_dFe=0;
            delete finite_volume->Be_scales_save;finite_volume->Be_scales_save=0;
            delete finite_volume->V;finite_volume->V=0;
            delete finite_volume->node_stiffness;finite_volume->node_stiffness=0;
            delete finite_volume->edge_stiffness;finite_volume->edge_stiffness=0;
            delete finite_volume->semi_implicit_data;finite_volume->semi_implicit_data=0;
            finite_volume->force_segments=0;finite_volume->twice_max_strain_per_time_step=0;
            finite_volume->Update_Be_Scales();}*/
        if(use_finite_volume){
            T stiffness=(T)2e4,damping=(T).01;
            NEO_HOOKEAN<T,3>* constitutive_model=new NEO_HOOKEAN<T,3>(stiffness,(T).45,damping,(T).25);
            STRAIN_MEASURE<TV,3>* strain_measure=new STRAIN_MEASURE<TV,3>(*tetrahedralized_volume);
            finite_volume=new FINITE_VOLUME<TV,3>(false,*strain_measure,*constitutive_model,0);}
        else{edge_springs->Set_Restlength_From_Particles();tet_springs->Set_Restlength_From_Particles();}}
}
//#####################################################################
// Function Preprocess_Solids_Substep
//#####################################################################
void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE
{
    if(dynamic_subsampling && !solids_parameters.triangle_collision_parameters.perform_self_collision) Update_Subsamples();
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    //RIGID_BODY_COLLISIONS<TV>& collisions=*solids_evolution->rigid_body_collisions;
    int frame=(int)(time*motion_frame_rate)+1;
    T alpha=time*motion_frame_rate-frame+1;
    if(test_number==1){
       //motion capture input
       for(int i=1;i<=arb.joint_mesh.joints.m;i++){JOINT<TV>& joint=*arb.joint_mesh.joints(i);
            RIGID_BODY<TV>* parent=arb.Parent(joint.id_number),*child=arb.Child(joint.id_number);
            int parent_index=id_to_index(parent->particle_index),child_index=id_to_index(child->particle_index);
            FRAME<TV> rigid_base_transform_parent1=rigid_base_transform,rigid_base_transform_child1=rigid_base_transform;
            FRAME<TV> rigid_base_transform_parent2=rigid_base_transform,rigid_base_transform_child2=rigid_base_transform;
            rigid_base_transform_parent1.t*=body_motion.trajectories(parent_index)(frame).length;
            rigid_base_transform_child1.t*=body_motion.trajectories(child_index)(frame).length;
            rigid_base_transform_parent2.t*=body_motion.trajectories(parent_index)(frame+1).length;
            rigid_base_transform_child2.t*=body_motion.trajectories(child_index)(frame+1).length;
            FRAME<TV> parent_frame=FRAME<TV>::Interpolation(body_motion.trajectories(parent_index)(frame).targeted_transform*rigid_base_transform_parent1,
                body_motion.trajectories(parent_index)(frame+1).targeted_transform*rigid_base_transform_parent2,alpha);
            FRAME<TV> child_frame=FRAME<TV>::Interpolation(body_motion.trajectories(child_index)(frame).targeted_transform*rigid_base_transform_child1,
                body_motion.trajectories(child_index)(frame+1).targeted_transform*rigid_base_transform_child2,alpha);
            if(joint.joint_function) joint.joint_function->Set_Target_Angle((joint.frame_jp*parent_frame.Inverse()*child_frame*joint.frame_cj).r);}}
    if(test_number==8){
        T desired_x=(T)two_pi/16;
        ROTATION<TV> desired_rotation=ROTATION<TV>(desired_x*sin(4*time),TV(0,1,0));
        for(int i=1;i<=arb.joint_mesh.joints.m;i++){JOINT<TV>& joint=*arb.joint_mesh.joints(i);
            if(joint.joint_function) joint.joint_function->Set_Target_Angle(desired_rotation);}}
    if(test_number==14){
        T sim_time=max((T)0,time-start_time);
        for(int i=1;i<=arb.joint_mesh.joints.m;i++){JOINT<TV>& joint=*arb.joint_mesh.joints(i);
            TV previous=Sidewinding_Position(sim_time,(i-1)*joint_separation);
            TV current=Sidewinding_Position(sim_time,i*joint_separation);
            TV next=Sidewinding_Position(sim_time,(i+1)*joint_separation);
            ROTATION<TV> previous_orientation=Upright_Orientation(previous,current); // parent
            ROTATION<TV> next_orientation=Upright_Orientation(current,next); // child
            ROTATION<TV> joint_frame=joint.frame_jp.r*initial_orientation.Inverse()*previous_orientation.Inverse()*next_orientation*initial_orientation*joint.frame_cj.r;
            if(time<start_time) joint_frame=joint_frame.Scale_Angle(time/start_time);
            if(joint.joint_function) joint.joint_function->Set_Target_Angle(joint_frame);}}
    if(test_number==16){
        T sim_time=max((T)0,time-start_time);
        for(int i=1;i<=arb.joint_mesh.joints.m;i++){JOINT<TV>& joint=*arb.joint_mesh.joints(i);
            ROTATION<TV> previous_orientation;
            ROTATION<TV> next_orientation(sin((T)two_pi*sim_time/period),TV(0,1,0));
            ROTATION<TV> joint_frame=joint.frame_jp.r*initial_orientation.Inverse()*previous_orientation.Inverse()*next_orientation*initial_orientation*joint.frame_cj.r;
            if(time<start_time) joint_frame=joint_frame.Scale_Angle(time/start_time);
            if(joint.joint_function) joint.joint_function->Set_Target_Angle(joint_frame);}}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    switch(test_number){
        case 1: Floppy_Human();break;
        case 2: Static_Human();break;
        case 8: Floppy_Fish();break;
        case 14: Sidewinding();break;
        case 16: Maggot();break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    // correct number nodes
    for(int i=1;i<=deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);soft_bindings.Set_Mass_From_Effective_Mass();

    T restlength_clamp=(T)1e-4;
    T cfl_strain_rate=(T)0.1;
    T stiffness=(T)2e6;
    T damping=(T).01;
    if(test_number==8) stiffness=(T)4e5;
    else if(test_number==14) stiffness=(T)2e5;
    else if(test_number==16) stiffness=(T)2e5;
    if(test_number==8) damping=(T).03;
    if(!use_finite_volume){stiffness=2e4;damping=(T)1;}
    if(use_implicit) solids_parameters.cfl*=10;
    
    for(int i=1;TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>(i);i++){
        if(use_finite_volume){
            finite_volume=Create_Finite_Volume(*tetrahedralized_volume,new NEO_HOOKEAN<T,3>(stiffness,(T).45,damping,(T).25),true,(T).1);
            solid_body_collection.Add_Force(finite_volume);}
        else{
            for(int i=1;i<=particles.array_collection->Size();i++){particles.mass(i)=(T)1;}
            T linear_stiffness=stiffness,linear_damping=damping;
            edge_springs=Create_Edge_Springs(*tetrahedralized_volume,linear_stiffness,linear_damping,strain_limit,cfl_strain_rate,true,(T)0,true,use_implicit);
            edge_springs->Clamp_Restlength(restlength_clamp); 
            solid_body_collection.Add_Force(edge_springs);
            T tet_stiffness=stiffness,tet_damping=damping;
            tet_springs=Create_Tet_Springs(*tetrahedralized_volume,tet_stiffness,tet_damping,false,(T).1,strain_limit,cfl_strain_rate,true,(T)0,true,use_implicit);
            tet_springs->Clamp_Restlength(restlength_clamp); 
            solid_body_collection.Add_Force(tet_springs);}}

    // disable strain rate CFL for all forces
    for(int i=1;i<=solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->limit_time_step_by_strain_rate=false;

    for(int i=1;i<=deformable_body_collection.deformable_geometry.structures.m;i++) if(!dynamic_cast<SEGMENTED_CURVE<TV>*>(deformable_body_collection.deformable_geometry.structures(i))){
        deformable_body_collection.collisions.collision_structures.Append(deformable_body_collection.deformable_geometry.structures(i));
        if(solids_parameters.triangle_collision_parameters.perform_self_collision && (!dynamic_cast<FREE_PARTICLES<TV>*>(deformable_body_collection.deformable_geometry.structures(i))))
            solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(deformable_body_collection.deformable_geometry.structures(i));}
    if(test_number==8){
        // collide structures with the ground only
        deformable_body_collection.collisions.Use_Structure_Collide_Collision_Body(true);
        COLLISION_GEOMETRY_ID ground_id=solid_body_collection.rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(ground->particle_index);
        for(int s=1;s<=deformable_body_collection.deformable_geometry.structures.m;s++) deformable_body_collection.collisions.structure_collide_collision_body(s).Set(ground_id);}
    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();

    // temporary
    if(test_number!=1 && test_number!=2) tests.Add_Gravity();
}
//#####################################################################
// Function Set_Particle_Is_Simulated
//#####################################################################
void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    if(dynamic_subsampling){
        static bool first_time=true;
        if(restart && first_time){ // rebuild the particle deletion list
            PHYSBAM_ASSERT(old_number_particles); // ensure Initialize_Dynamic_Subsampling was called 
            LOG::cout<<"DEBUG rebuilding particle deletion list, particles.array_collection->Size()="<<particles.array_collection->Size()<<", old_number_particles="<<old_number_particles<<std::endl; 
            for(int p=old_number_particles+1;p<=particles.array_collection->Size();p++) particles.array_collection->Add_To_Deletion_List(p);
            binding_list.Clean_Memory();soft_bindings.Clean_Memory(); // clear bindings: assumes all bindings are for dynamic samples
            LOG::cout<<"DEBUG number of particles in the deletion list="<<particles.array_collection->deletion_list.m<<std::endl;
            first_time=false;}
        particle_is_simulated.Subset(particles.array_collection->deletion_list).Fill(false);
        for(int p=1;p<=rigid_body_collection.rigid_body_particle.array_collection->Size();p++){
            int particle_to_exclude=particles.array_collection->Size()+p;
            bool exclude;
            if(particle_to_exclude<=deformable_body_collection.particles.array_collection->Size()) exclude=false;
            else if(!rigid_body_collection.Is_Active(particle_to_exclude-deformable_body_collection.particles.array_collection->Size())) exclude=true;
            else exclude=!rigid_body_collection.Rigid_Body(particle_to_exclude-deformable_body_collection.particles.array_collection->Size()).Is_Simulated();}}
}
//#####################################################################
// Function Initialize_Joint_Between
//#####################################################################
void Initialize_Joint_Between(JOINT<TV>* joint,const RIGID_BODY<TV>& parent,const RIGID_BODY<TV>& child,TV up)
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    arb.joint_mesh.Add_Articulation(parent.particle_index,child.particle_index,joint);
    const FRAME<TV>& pf=parent.Frame();
    const FRAME<TV>& cf=child.Frame();
    TV x=(cf.t-pf.t).Normalized();
    TV y=up.Projected_Orthogonal_To_Unit_Direction(x).Normalized();
    TV z=TV::Cross_Product(x,y);
    FRAME<TV> J((T).5*(cf.t+pf.t),ROTATION<TV>(MATRIX<T,3>(x,y,z)));
    joint->Set_Parent_To_Joint_Frame(J.Inverse_Times(pf));
    joint->Set_Child_To_Joint_Frame(J.Inverse_Times(cf));
}
//#####################################################################
// Function Create_Joints_From_Hierarchy
//#####################################################################
void Create_Joints_From_Hierarchy(int parent)
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    for(int i=1;i<=body_motion.bone_hierarchy(parent).m;i++){
        int child=body_motion.name_to_track_index.Get(body_motion.bone_hierarchy(parent)(i));
        JOINT<TV>* joint=new POINT_JOINT<TV>;
        TV up=(arb.rigid_body_collection.Rigid_Body(rigid_body_ids(parent)).Frame()*TV(0,0,1)+arb.rigid_body_collection.Rigid_Body(rigid_body_ids(child)).Frame()*TV(0,0,1))/2.;
        Initialize_Joint_Between(joint,arb.rigid_body_collection.Rigid_Body(rigid_body_ids(parent)),arb.rigid_body_collection.Rigid_Body(rigid_body_ids(child)),up);
        JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
        joint_function->Set_k_p(1000);joint_function->Set_Target_Angle(ROTATION<TV>());
        Create_Joints_From_Hierarchy(child);}
}
//#####################################################################
// Function Floppy_Human
//#####################################################################
void Floppy_Human()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;

    solids_parameters.enforce_poststabilization_in_cg=false;
    solids_parameters.cfl=(T)1;
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solid_body_collection.print_residuals=true;
    solids_parameters.rigid_body_evolution_parameters.rigid_geometry_evolution_parameters.use_kinematic_keyframes=false;
    solids_parameters.triangle_collision_parameters.perform_self_collision=true;
    plastic=false;use_finite_volume=true;
    strain_limit=true;use_implicit=false;
    bool use_embedding=false; //can't do stacked embeddings....

    arb.Set_Iterative_Tolerance((T)1e-4);
    // prestabilization settings
    arb.Set_Contact_Level_Iterations(5);
    arb.Set_Shock_Propagation_Level_Iterations(5);
    arb.Set_Use_Shock_Propagation(false);
    arb.Set_Do_Final_Pass(false);
    // poststabilization settings
    arb.Set_Poststabilization_Iterations(5);
    arb.poststabilization_projection_iterations=2;

    T density=1000;

    //adding bones
    id_to_index.Resize(int(body_motion.trajectories.m));
    rigid_body_ids.Resize(body_motion.trajectories.m);
    for(int i=1;i<=body_motion.trajectories.m;i++){
        T scale=(T).9,length=body_motion.trajectories(i)(1).length,height=scale*length,radius=(T).18;
        RIGID_BODY<TV>& rigid_body=tests.Add_Analytic_Cylinder(height,radius);
        int id=rigid_body.particle_index;assert(id);
        if(id_to_index.Size()<id) id_to_index.Resize(id);
        id_to_index(id)=i;rigid_body_ids(i)=id;
        T volume=(T)pi*radius*radius*height;
        rigid_body.Set_Mass(density*volume);
        rigid_body.Update_Bounding_Box();
        if(i==1){
            T_INERTIA_TENSOR& inertia_tensor=arb.rigid_body_collection.rigid_body_particle.inertia_tensor(rigid_body_ids(i));
            T max=inertia_tensor(1),translation=rigid_body.axis_aligned_bounding_box.min_corner.x;int axis=1;
            if(max<inertia_tensor(2)){max=inertia_tensor(2);axis=2;}
            if(max<inertia_tensor(3)){max=inertia_tensor(3);axis=3;}
            axis=3; //manually set for now
            if(axis==2){translation=rigid_body.axis_aligned_bounding_box.min_corner.y;rigid_base_transform=FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,0,-1)));}
            if(axis==3){translation=rigid_body.axis_aligned_bounding_box.min_corner.z;rigid_base_transform=FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0)));}
            rigid_base_transform=FRAME<TV>(TV(-1*translation,0,0),ROTATION<TV>())*rigid_base_transform;
            rigid_base_transform.t/=body_motion.trajectories(1)(1).length;}
        FRAME<TV> rigid_base_transform_i=rigid_base_transform;
        rigid_base_transform_i.t*=length;
        rigid_body.Set_Frame(body_motion.trajectories(i)(1).targeted_transform*rigid_base_transform_i);}
    
    Create_Joints_From_Hierarchy(body_motion.name_to_track_index.Get("Root"));

    // add the person
    RIGID_BODY_STATE<TV> human_state(FRAME<TV>(TV(0,0,0),ROTATION<TV>()));
    if(!use_embedding) tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/body_150_scaled_10.tet",human_state,false,false,1000,(T)1);
    else{
        TETRAHEDRALIZED_VOLUME<T>& volume_original=*TETRAHEDRALIZED_VOLUME<T>::Create();
        TETRAHEDRALIZED_VOLUME<T>& volume=*(TETRAHEDRALIZED_VOLUME<T>*)volume_original.Append_Particles_And_Create_Copy(particles);
        FILE_UTILITIES::Read_From_File<T>(data_directory+"/Tetrahedralized_Volumes/body_150_scaled_10.tet",volume);
        int offset=particles.array_collection->Size();

        TRIANGULATED_SURFACE<T>& surface_original=*TRIANGULATED_SURFACE<T>::Create();
        FILE_UTILITIES::Read_From_File<T>(data_directory+"/Triangulated_Surfaces/body_scaled_1.tri",surface_original);
        TRIANGULATED_SURFACE<T>& surface=*(TRIANGULATED_SURFACE<T>*)surface_original.Append_Particles_And_Create_Copy(particles);
        surface.Update_Number_Nodes();volume.Update_Number_Nodes();
        
        ARRAY<PAIR<int,TV> > bindings; 
        if(FILE_UTILITIES::File_Exists(data_directory+"/bindings/body_bindings")) FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/bindings/body_bindings",bindings);
        else{
            ARRAY<int> tets;ARRAY<PAIR<int,TV> > bindings;const T tolerance=(T)1e-4;
            volume.Initialize_Hierarchy();
            for(int p=1;p<=surface_original.particles.array_collection->Size();p++){
                tets.Remove_All();volume.hierarchy->Intersection_List(surface_original.particles.X(p),tets,tolerance);bool got_bind=false;
                for(int tt=1;tt<=tets.m;tt++){int t=tets(tt);
                    TV bary=TETRAHEDRON<T>::First_Three_Barycentric_Coordinates(surface_original.particles.X(p),volume.particles.X.Subset(volume.mesh.elements(t)));
                    if(bary.x>-tolerance && bary.y>-tolerance && bary.z>-tolerance && bary.x+bary.y+bary.z<(T)1+tolerance){bindings.Append(PAIR<int,TV>(t,bary));got_bind=true;break;}}
                if(!got_bind){LOG::cout<<"no binding on particle "<<p<<std::endl;bindings.Append(PAIR<int,TV>(0,TV(0,0,0)));}}
            FILE_UTILITIES::Write_To_File(stream_type,data_directory+"/bindings/body_bindings",bindings);}
       for(int i=1;i<=bindings.m;i++){int p=offset+i;
            if(bindings(i).x==0) continue;
            VECTOR<int,4> nodes=volume.mesh.elements(bindings(i).x);
            binding_list.Add_Binding(new LINEAR_BINDING<TV,4>(particles,p,nodes,bindings(i).y));}
        particles.Store_Velocity();
        LOG::cout<<offset<<std::endl;
        SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(volume,density,true);
        deformable_body_collection.deformable_geometry.Add_Structure(&volume);
        deformable_body_collection.deformable_geometry.Add_Structure(&surface);}
    
    // binding the deformable particles to the rigid bodies
    for(int p=1;p<=rigid_body_collection.rigid_body_particle.array_collection->Size();p++) tests.Bind_Particles_In_Rigid_Body(rigid_body_collection.Rigid_Body(p));

    arb.Use_PD_Actuators();
    arb.global_post_stabilization=true;
    arb.poststabilization_projection_iterations=5;
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=5;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
}
//#####################################################################
// Function Static_Human
//#####################################################################
void Static_Human()
{
    solids_parameters.enforce_poststabilization_in_cg=false;
    solids_parameters.cfl=(T)2;
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
    solid_body_collection.print_residuals=true;
    plastic=false;use_finite_volume=true;
    strain_limit=false;use_implicit=true;

    //T friction=(T).2;
    T density=1000;

    // add the person
    RIGID_BODY_STATE<TV> human_state(FRAME<TV>(TV(0,0,0),ROTATION<TV>()));
    tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/body.tet",human_state,false,false,1000,(T)1);

    RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("sphere",(T).2,(T).1);    
    rigid_body.X()=TV((T)1,(T).2,0);
    rigid_body.Twist().linear=TV(-2.5,0,0);
    rigid_body.Update_Bounding_Box();
    T volume=(T)(4./3.*pi*.1*.1*.1);
    rigid_body.Set_Mass(density*volume);
    
    //ground=&tests.Add_Ground(friction,0,0);

    solids_parameters.implicit_solve_parameters.cg_projection_iterations=5;
}
//#####################################################################
// Function Floppy_Fish
//#####################################################################
void Floppy_Fish()
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;

    solids_parameters.enforce_poststabilization_in_cg=false;
    solids_parameters.cfl=(T)2;

    arb.Set_Iterative_Tolerance((T)1e-4);
    // prestabilization settings
    arb.Set_Contact_Level_Iterations(5);
    arb.Set_Shock_Propagation_Level_Iterations(5);
    arb.Set_Use_Shock_Propagation(false);
    arb.Set_Do_Final_Pass(false);
    // poststabilization settings
    arb.Set_Poststabilization_Iterations(5);
    arb.poststabilization_projection_iterations=2;

    T friction=(T).2;
    TV center(0,3,0);
//    T scale=(T).55;
//    T spacing=(T).25;
//    T plank_length=(T)1;
//    T extent=scale*(number_of_joints*(plank_length+spacing)+plank_length);

    ARRAY<RIGID_BODY<TV>*> bones;

    bones.Append(&tests.Add_Rigid_Body("miniplank25wide2",(T)1,friction));
    bones.Last()->Set_Frame(FRAME<TV>(TV(-(T)2,(T)3,0)));

    bones.Append(&tests.Add_Rigid_Body("miniplank25wide2",(T)1,friction));
    bones.Last()->Set_Frame(FRAME<TV>(TV(-(T).75,(T)3,0)));

    bones.Append(&tests.Add_Rigid_Body("miniplank25wide2",(T)1,friction));
    bones.Last()->Set_Frame(FRAME<TV>(TV((T).5,(T)3,0)));

    bones.Append(&tests.Add_Rigid_Body("miniplank25wide2",(T).7,friction));
    bones.Last()->Set_Frame(FRAME<TV>(TV((T)1.7,(T)3,0)));

    bones.Append(&tests.Add_Rigid_Body("miniplank25wide2",(T).5,friction));
    bones.Last()->Set_Frame(FRAME<TV>(TV((T)2.5,(T)3,0)));

    T density=1000;
    for(int i=1;i<=bones.m;i++){
        bones(i)->Set_Mass(density*bones(i)->Volume());}

    T joint_strengths[4]={(T)1000,(T)1000,(T)400,(T)200};
    for(int i=2;i<=bones.m;i++){
        JOINT<TV>* joint=new POINT_JOINT<TV>;
        Initialize_Joint_Between(joint,*bones(i-1),*bones(i),TV(0,0,1));
        JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
        joint_function->Set_k_p(joint_strengths[i-2]);joint_function->Set_Target_Angle(ROTATION<TV>());}

//    // make the joints
//    tests.PD_Curl(scale,center-(T).5*(extent-scale*plank_length)*TV::Axis_Vector(1),ROTATION<TV>(),25,number_of_joints,false);

    // add the fish
    RIGID_BODY_STATE<TV> fish_state(FRAME<TV>(TV(0,3,0),ROTATION<TV>((T)half_pi,TV(1,0,0))));
    if(fish_mattress){
        BOX<TV> box(-TV((T)3.93278,(T)1.07277,(T)0.384066),TV((T)2.68344,(T)1.1747,(T)0.384353));VECTOR<int,3> counts(20,15,5);
        GRID<TV> mattress_grid=GRID<TV>(counts,box);
        tests.Create_Mattress(mattress_grid,true,&fish_state);}
    else tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/fish_42K.tet",fish_state,true,false,1000,(T)1);

    // binding the deformable particles to the rigid bodies
    for(int p=1;p<=rigid_body_collection.rigid_body_particle.array_collection->Size();p++) tests.Bind_Particles_In_Rigid_Body(rigid_body_collection.Rigid_Body(p));

    ground=&tests.Add_Ground(friction,0,0);

    arb.Use_PD_Actuators();
    arb.global_post_stabilization=true;
    arb.poststabilization_projection_iterations=5;
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=5;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
}
//#####################################################################
// Function PD_Snake
//#####################################################################
void PD_Snake(const T scale,const TV shift,const ROTATION<TV> orient,const T k_p,const int number_of_joints,const T space_adjustment,const T friction=.5)
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    RIGID_BODY<TV> *parent_body=0,*child_body=0;
    T cheight=(T)0;
    initial_orientation=ROTATION<TV>((T)half_pi,TV(0,0,1));

    // Create first body
    parent_body=&tests.Add_Rigid_Body("cyllink",scale*(T).2,friction);
    parent_body->X()=shift;
    parent_body->Rotation()=orient*initial_orientation;
    parent_body->Set_Name("parent");
    parent_body->Set_Mass(50);

    // Add children and joints
    T desired_x=(T)two_pi/(T)(number_of_joints+1);
    for(int i=1;i<=number_of_joints;i++){
        cheight+=scale*((T)1.25+space_adjustment);
        child_body=&tests.Add_Rigid_Body("cyllink",scale*(T).2,friction);
        child_body->X()=orient.Rotate(TV(cheight,0,0))+shift;
        child_body->Rotation()=orient*initial_orientation;
        child_body->Set_Coefficient_Of_Restitution((T)0.5);
        child_body->Set_Name(STRING_UTILITIES::string_sprintf("child_%d",i));
        child_body->Set_Mass(50);

        ROTATION<TV> desired_rotation=ROTATION<TV>(desired_x,TV());

        JOINT<TV>* joint=new POINT_JOINT<TV>();arb.joint_mesh.Add_Articulation(child_body->particle_index-1,child_body->particle_index,joint);
        JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
        joint_function->Set_k_p(k_p);joint_function->Set_Target_Angle(desired_rotation);
        joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(0,-scale*((T).625+space_adjustment/2),0),ROTATION<TV>(-0*(T)pi/2,TV(0,1,0))));
        joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,scale*((T).625+space_adjustment/2),0),ROTATION<TV>(-0*(T)pi/2,TV(0,1,0))));

        parent_body=child_body;}
}
//#####################################################################
// Function Sidewinding
//#####################################################################
void Sidewinding()
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
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
    number_of_joints=12;
    TV center(0,(T)1.5,0);
    T scale=(T).55;
    T spacing=(T).25;
    T plank_length=(T)1;
    T extent=scale*(number_of_joints*(plank_length+spacing)+plank_length);
    ROTATION<TV> snake_rotation(-(T)pi/8,TV(0,1,0));

    joint_separation=(T)1.25*scale;
    wavelength=(T)8*joint_separation;
    period=(T)2;
    height_amplitude=(T).4;
    bend_amplitude=(T)1.4;
    start_time=(T)1.5;
    T friction=(T).8;

    // make the joints
    PD_Snake(scale,center-(T).5*(extent-scale*plank_length)*snake_rotation.Rotate(TV::Axis_Vector(1)),snake_rotation,2000,number_of_joints,0);

    tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/snake_8K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(center,snake_rotation)),false,false,1000,(T).5);

    // binding the deformable particles to the rigid bodies
    for(int p=1;p<=rigid_body_collection.rigid_body_particle.array_collection->Size();p++) tests.Bind_Particles_In_Rigid_Body(rigid_body_collection.Rigid_Body(p));

//     {RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("subdivided_box",(T).3,friction);
//     rigid_body.X()=TV((T)0,(T).3,-10);
//     rigid_body.Set_Name("obstruction");
//     rigid_body.Set_Mass(30);}

//     {RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("subdivided_box",(T).3,friction);
//     rigid_body.X()=TV((T)-2,(T).3,-10);
//     rigid_body.Set_Name("obstruction");
//     rigid_body.Set_Mass(30);}

    int steps[]={1,2,3,4,5,4,3,3,4,5,6,2,1};
    int count=sizeof steps/sizeof*steps;

    for(int i=1;i<=count;i++) for(int j=1;j<=6;j++){
        RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("plank",(T)1,friction);
        rigid_body.X()=TV((T)(-(T)35+10*j),(T).2*steps[i-1]-(T).1,(T)(-6-2*i));
        rigid_body.Rotation()=ROTATION<TV>((T)half_pi,TV(0,1,0));
        rigid_body.Set_Name("obstruction");
        rigid_body.is_static=true;
        rigid_body.Set_Mass(30);}

    ground=&tests.Add_Ground(friction,0,0);

    arb.Use_PD_Actuators();
    arb.global_post_stabilization=true;
    arb.poststabilization_projection_iterations=5;
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=5;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
}
//#####################################################################
// Function Add_Maggot
//#####################################################################
void Add_Maggot(const T scale,const RIGID_BODY_STATE<TV>& state,const std::string& resolution,T friction=(T).2)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    number_of_joints=2;
    T spacing=(T).05;
    T plank_length=(T)1;
    T extent=scale*(number_of_joints*(plank_length+spacing)+plank_length);
    ROTATION<TV> snake_rotation(-(T)pi/8,TV(0,1,0));

    joint_separation=(T)1.05*scale;

    // make the joints
    PD_Snake(scale,state.frame.t-(T).5*(extent-scale*plank_length)*state.frame.r.Rotate(TV::Axis_Vector(1)),state.frame.r,2000*scale,number_of_joints,-(T).2,friction);

    tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/maggot_"+resolution+"K.tet",state,false,false,1000,(T)1.2*scale);

    // binding the deformable particles to the rigid bodies
    for(int p=rigid_body_collection.rigid_body_particle.array_collection->Size()-number_of_joints;p<=rigid_body_collection.rigid_body_particle.array_collection->Size();p++){
        RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(p);
        rigid_body.Twist().angular=state.twist.angular;
        rigid_body.Twist().linear=state.Pointwise_Object_Velocity(rigid_body.X());
        tests.Bind_Particles_In_Rigid_Body(rigid_body);}
}
//#####################################################################
// Function Maggot
//#####################################################################
void Maggot()
{
    solids_parameters.cfl=(T)1;

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

    wavelength=(T)8*joint_separation;
    period=(T)1;
    start_time=(T)1.5;

    Add_Maggot((T).65,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)1.5,0))),"8");

    ground=&tests.Add_Ground((T).8,0,0);

    arb.Use_PD_Actuators();
    arb.global_post_stabilization=true;
    arb.poststabilization_projection_iterations=5;
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=5;
}
//#####################################################################
// Function Initialize_Dynamic_Subsampling
//#####################################################################
void Initialize_Dynamic_Subsampling()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;

    for(int i=1;TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>(i);i++){
        tetrahedralized_volume->Update_Number_Nodes();
        if(!tetrahedralized_volume->triangulated_surface) tetrahedralized_volume->Initialize_Triangulated_Surface();
        surface_elements.Append_Elements(tetrahedralized_volume->triangulated_surface->mesh.elements);}

    for(int i=1;TRIANGULATED_SURFACE<T>* triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>*>(i);i++){
        triangulated_surface->Update_Number_Nodes();
        surface_elements.Append_Elements(triangulated_surface->mesh.elements);}

    surface_elements.Flattened().Get_Unique(surface_particles);

    FREE_PARTICLES<TV>& free_particles=*FREE_PARTICLES<TV>::Create(deformable_body_collection.particles);
    deformable_body_collection.deformable_geometry.Add_Structure(&free_particles);

    old_number_particles=deformable_body_collection.particles.array_collection->Size();
    triangle_free_particles.Resize(surface_elements.m);

//    solids_parameters.write_static_variables_every_frame=true;
}
//#####################################################################
// Function Add_Subsamples
//#####################################################################
void Add_Subsamples(const int surface_triangle_index,ARRAY<BINDING<TV>*>& new_binding_list,ARRAY<VECTOR<int,2> >& new_soft_bindings)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;

    const VECTOR<int,3>& triangle=surface_elements(surface_triangle_index);
    ARRAY<int>& particle_subsamples=triangle_free_particles(surface_triangle_index);
    for(int i=1;i<=subsamples;i++){
        VECTOR<T,3> weights;do{weights=random_numbers.Get_Uniform_Vector(TV(),TV(1,1,1));}while(weights.x+weights.y>=1);weights.z=(T)1-weights.x-weights.y;
        int hard_bound_particle=particles.array_collection->Add_Element_From_Deletion_List();
        new_binding_list.Append(new LINEAR_BINDING<TV,3>(particles,hard_bound_particle,triangle,weights));
        int soft_bound_particle=particles.array_collection->Add_Element_From_Deletion_List();
        new_soft_bindings.Append(VECTOR<int,2>(soft_bound_particle,hard_bound_particle));
        particle_subsamples.Append(soft_bound_particle);}
}
//#####################################################################
// Function Delete_Subsamples
//#####################################################################
void Delete_Subsamples(const int surface_triangle_index)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    // only need to delete the particles
    ARRAY<int>& particle_subsamples=triangle_free_particles(surface_triangle_index);
    for(int i=1;i<=particle_subsamples.m;i++){
        int soft_bound_particle=particle_subsamples(i);
        int hard_bound_particle=soft_bindings.Parent(soft_bound_particle);
        particles.array_collection->Add_To_Deletion_List(soft_bound_particle);
        particles.array_collection->Add_To_Deletion_List(hard_bound_particle);}
    particle_subsamples.Remove_All();
}
//#####################################################################
// Function Persist_Subsamples
//#####################################################################
void Persist_Subsamples(const int surface_triangle_index,ARRAY<BINDING<TV>*>& new_binding_list,ARRAY<VECTOR<int,2> >& new_soft_bindings)
{
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    ARRAY<int>& particle_subsamples=triangle_free_particles(surface_triangle_index);
    if(!particle_subsamples.m) return Add_Subsamples(surface_triangle_index,new_binding_list,new_soft_bindings);

    assert(particle_subsamples.m==subsamples);
    for(int i=1;i<=particle_subsamples.m;i++){
        int soft_bound_particle=particle_subsamples(i);
        // save soft binding
        const VECTOR<int,2>& soft_binding=soft_bindings.bindings(soft_bindings.Soft_Binding(soft_bound_particle));
        new_soft_bindings.Append(soft_binding);
        // save hard binding
        int hard_binding_index=binding_list.binding_index_from_particle_index(soft_binding.y);
        new_binding_list.Append(binding_list.bindings(hard_binding_index));
        binding_list.bindings(hard_binding_index)=0;} // so that binding won't be deleted when binding_list is rebuilt
}
//#####################################################################
// Function Update_Subsamples
//#####################################################################
void Update_Subsamples()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    COLLISION_GEOMETRY_COLLECTION<TV>& collision_body_list=deformable_body_collection.collisions.collision_body_list;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    FREE_PARTICLES<TV>& free_particles=deformable_body_collection.deformable_geometry.template Find_Structure<FREE_PARTICLES<TV>&>();

    ARRAY<BINDING<TV>*> new_binding_list;
    ARRAY<VECTOR<int,2> > new_soft_bindings;

    LOG::SCOPE("Update_Subsamples");

    if(temporarily_disable_dynamic_subsampling) return;

    binding_list.Clamp_Particles_To_Embedded_Positions();

    // the minimum distance of each particle to a collision object
    particle_distances.Resize(old_number_particles);
    particle_distances.Subset(surface_particles).Fill(FLT_MAX);

    for(int i=1;i<=surface_particles.m;i++){int p=surface_particles(i);
        for(COLLISION_GEOMETRY_ID body(1);body<=collision_body_list.Size();body++)
            particle_distances(p)=PhysBAM::min(particle_distances(p),collision_body_list(body).Implicit_Geometry_Extended_Value(particles.X(p)));}

    // iterate over surface elements
    for(int t=1;t<=surface_elements.m;t++){
        const VECTOR<int,3>& triangle=surface_elements(t);
        T triangle_distance=particle_distances.Subset(triangle).Min();
        if(triangle_distance<refinement_distance) Persist_Subsamples(t,new_binding_list,new_soft_bindings);
        else Delete_Subsamples(t);}

    // rebuild free particles
    free_particles.nodes.Remove_All();
    for(int t=1;t<=surface_elements.m;t++){
        ARRAY<int>& subsamples=triangle_free_particles(t);
        if(subsamples.m) free_particles.nodes.Append_Elements(subsamples);}

    // rebuild hard bindings
    binding_list.Clean_Memory();
    for(int b=1;b<=new_binding_list.m;b++) binding_list.Add_Binding(new_binding_list(b));

    // rebuild soft bindings
    soft_bindings.Clean_Memory();
    for(int b=1;b<=new_soft_bindings.m;b++) soft_bindings.Add_Binding(new_soft_bindings(b),true);

    binding_list.Clamp_Particles_To_Embedded_Positions();binding_list.Clamp_Particles_To_Embedded_Velocities();
    soft_bindings.Clamp_Particles_To_Embedded_Positions();soft_bindings.Clamp_Particles_To_Embedded_Velocities();

    // correct number nodes
    for(int i=1;i<=deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();
    
    // correct mass
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();

    // reinitialize deformable object collisions
    deformable_body_collection.collisions.Initialize_Object_Collisions(solids_parameters.deformable_object_collision_parameters.collide_with_interior,
        solids_parameters.deformable_object_collision_parameters.collision_tolerance,
        solids_parameters.deformable_object_collision_parameters.use_spatial_partition_for_levelset_collision_objects,
        solids_parameters.deformable_object_collision_parameters.disable_multiple_levelset_collisions,
        solids_parameters.deformable_object_collision_parameters.maximum_levelset_collision_projection_velocity);

    // reinitialize fragments
    solid_body_collection.Update_Simulated_Particles();
}
//#####################################################################
};
}
#endif

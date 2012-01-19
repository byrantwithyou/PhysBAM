//#####################################################################
// Copyright 2007-2008, Jon Gretarsson, Michael Lentine, Avi Robinson-Mosher, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS_WATER
//#####################################################################
//   1. Human in water
//   2. Human in water with rigid mesh
//   9. Fish flopping in water
//#####################################################################
#ifndef __SWIMMING_TESTS_WATER__
#define __SWIMMING_TESTS_WATER__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_TET_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Dynamics/Motion/BODY_MOTION_SEQUENCE.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Standard_Tests/WATER_STANDARD_TESTS_3D.h>
namespace PhysBAM{

template<class T_input>
class SWIMMING_TESTS_WATER:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
    typedef GRID<TV> T_GRID;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename RIGID_BODY_POLICY<TV>::INERTIA_TENSOR T_INERTIA_TENSOR;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;using BASE::output_directory;using BASE::last_frame;using BASE::frame_rate;
    using BASE::stream_type;using BASE::data_directory;using BASE::solid_body_collection;using BASE::parse_args;using BASE::test_number;using BASE::resolution;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions; // silence -Woverloaded-virtual
    using BASE::Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization;using BASE::mpi_world;using BASE::Add_Thin_Shell_To_Fluid_Simulation;
    using BASE::Add_Volumetric_Body_To_Fluid_Simulation;using BASE::Add_To_Fluid_Simulation;

    WATER_STANDARD_TESTS_3D<T_GRID > water_tests;
    SOLIDS_STANDARD_TESTS<TV> solids_tests;
    T_GRID mattress_grid;
    int deformable_object_id;
    T solid_density;
    T initial_fluid_height;
    bool implicit_springs;
    MATRIX<T,4> world_to_source;
    int bodies;
    int sub_test;
    int motion_frame_rate;
    TV left_corner,right_corner;
    BODY_MOTION_SEQUENCE<T> body_motion;
    ARRAY<int,int> id_to_index;
    ARRAY<int> rigid_body_ids;
    FRAME<TV> rigid_base_transform;

    ARRAY<PAIR<int,TV> > constrained_node_positions;
    ARRAY<TRIPLE<int,T,TV> > constrained_nodes;
    ARRAY<int> rigid_bodies_to_simulate;
    ARRAY<DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>*> deformable_objects_to_simulate;
    ARRAY<int> rigid_bodies_to_collide_against;
    ORIENTED_BOX<TV> fish_bounding_box;
    LEVELSET_IMPLICIT_OBJECT<TV>* fish_levelset,*human_levelset;

    SWIMMING_TESTS_WATER(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.WATER),
        water_tests(*this,fluids_parameters,solid_body_collection.rigid_body_collection),
        solids_tests(*this,solid_body_collection),deformable_object_id(0),solid_density((T)2000),
        initial_fluid_height((T)0),implicit_springs(false),world_to_source(MATRIX<T,4>::Identity_Matrix()),bodies(5),sub_test(1),fish_levelset(0),human_levelset(0)
    {
    }

    // Unused callbacks
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Initialize_Euler_State() PHYSBAM_OVERRIDE {}
    void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Extrapolate_Phi_Into_Objects(const T time) PHYSBAM_OVERRIDE {}
    void Get_Source_Reseed_Mask(T_ARRAYS_BOOL*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE {}
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_String_Argument("-mtn","Specify the physbam motion file to use");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    water_tests.Initialize(Water_Test_Number(test_number),resolution);
    bool solid_node=mpi_world->initialized && !mpi_world->rank;
    bool mpi=mpi_world->initialized;
    LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
    last_frame=1000;

    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;

    fluids_parameters.incompressible_iterations=200;
    *fluids_parameters.grid=water_tests.grid;
    fluids_parameters.fluid_affects_solid=fluids_parameters.solid_affects_fluid=true;
    fluids_parameters.second_order_cut_cell_method=true;
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solids_parameters.verbose_dt=true;
    solid_body_collection.print_residuals=false;
    solids_parameters.rigid_body_evolution_parameters.maximum_rigid_body_time_step_fraction=(T)1;

    if(solid_node || !mpi) solids_parameters.use_rigid_deformable_contact=true;
    output_directory=STRING_UTILITIES::string_sprintf("Swimming_Tests_Water/Test_%d_Resolution_%d",test_number,resolution);
        
    if(parse_args->Is_Value_Set("-mtn")) FILE_UTILITIES::Read_From_File(stream_type,parse_args->Get_String_Value("-mtn"),body_motion);
        
    fluids_parameters.domain_walls[2][2]=fluids_parameters.domain_walls[2][1]=false;
    fluids_parameters.density=(T)1000;
    solids_parameters.implicit_solve_parameters.cg_iterations=400;

    T scale;
    switch(test_number){
        case 1:
            frame_rate=30;motion_frame_rate=120;
            fluids_parameters.reincorporate_removed_particle_velocity=true;
            fluids_parameters.removed_particle_mass_scaling=60;
            fluids_parameters.density=(T)1000;
            solids_parameters.implicit_solve_parameters.cg_iterations=800;
            fluids_parameters.domain_walls[2][2]=false;
            initial_fluid_height=(T)5.5;
            scale=2;
            fluids_parameters.grid->Initialize(32*resolution+1,(int)(36*resolution/scale)+1,24*resolution+1,(T)-8*scale,(T)8*scale,(T)0,(T)18,(T)-6*scale,(T)6*scale);
            PHYSBAM_ASSERT(fluids_parameters.grid->dX.x==fluids_parameters.grid->dX.y && fluids_parameters.grid->dX.y==fluids_parameters.grid->dX.z);
            break;
        case 9:
            fluids_parameters.reincorporate_removed_particle_velocity=true;
            fluids_parameters.removed_particle_mass_scaling=60;
            fluids_parameters.density=(T)1000;
            solids_parameters.implicit_solve_parameters.cg_iterations=800;
            fluids_parameters.domain_walls[2][2]=false;
            initial_fluid_height=(T)3.5;
            fluids_parameters.grid->Initialize(32*resolution+1,36*resolution+1,24*resolution+1,(T)-8,(T)8,(T)0,(T)18,(T)-6,(T)6);
            PHYSBAM_ASSERT(fluids_parameters.grid->dX.x==fluids_parameters.grid->dX.y && fluids_parameters.grid->dX.y==fluids_parameters.grid->dX.z);
            break;
        default:
            LOG::cerr<<"Unrecognized test number "<<test_number<<std::endl;exit(1);
    }

    switch(test_number){
        case 1:THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Rigid_Body_Walls(*this,(T).5,(T).5,&rigid_bodies_to_collide_against);break;
        case 9:THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Rigid_Body_Walls(*this,(T).5,(T).5,&rigid_bodies_to_collide_against);break;
        default:THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Rigid_Body_Walls(*this);break;}

    fluids_parameters.domain_walls[2][1]=true;

    // give mon hints
    LOG::cout<<"MONITOR begin_frame="<<this->first_frame<<std::endl;
    LOG::cout<<"MONITOR output_directory="<<(FILE_UTILITIES::Get_Working_Directory()+"/"+output_directory)<<std::endl;
    LOG::cout<<"MONITOR end_frame="<<last_frame<<std::endl;
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Water_Test_Number
//#####################################################################
static int Water_Test_Number(const int test_number)
{
    switch(test_number){
        case 1:
        case 9:
            return 1;
        default:
            return 1;}
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    int frame=(int)(time*motion_frame_rate)+1;
    T alpha=time*motion_frame_rate-frame+1;
    if(test_number==1){
        for(int i=0;i<arb.joint_mesh.joints.m;i++){JOINT<TV>& joint=*arb.joint_mesh.joints(i);
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
    if(test_number==9){
        T desired_x=(T)two_pi/16;
        ROTATION<TV> desired_rotation=ROTATION<TV>(desired_x*sin(4*time),TV(0,1,0));
        for(int i=0;i<arb.joint_mesh.joints.m;i++){JOINT<TV>& joint=*arb.joint_mesh.joints(i);
            if(joint.joint_function) joint.joint_function->Set_Target_Angle(desired_rotation);}}
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
        Initialize_Joint_Between(joint,arb.rigid_body_collection.Rigid_Body(rigid_body_ids(parent)),arb.rigid_body_collection.Rigid_Body(rigid_body_ids(child)),TV(0,1,0));
        JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
        joint_function->Set_k_p(1000);joint_function->Set_Target_Angle(ROTATION<TV>());
        Create_Joints_From_Hierarchy(child);}
}
//#####################################################################
// Function Floppy_Human
//#####################################################################
void Floppy_Human()
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;

    solids_parameters.enforce_poststabilization_in_cg=false;
    solids_parameters.cfl=(T)2;
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solid_body_collection.print_residuals=true;
    //plastic=false;use_finite_volume=true;
    //strain_limit=false;use_implicit=false;

    arb.Set_Iterative_Tolerance((T)1e-4);
    // prestabilization settings
    arb.Set_Contact_Level_Iterations(5);
    arb.Set_Shock_Propagation_Level_Iterations(5);
    arb.Set_Use_Shock_Propagation(false);
    arb.Set_Do_Final_Pass(false);
    // poststabilization settings
    arb.Set_Poststabilization_Iterations(5);
    arb.poststabilization_projection_iterations=2;

    //T friction=(T).2;
    T density=1000;

    //adding bones
    id_to_index.Resize(int(body_motion.trajectories.m));
    rigid_body_ids.Resize(body_motion.trajectories.m);
    for(int i=0;i<body_motion.trajectories.m;i++){
        T scale=(T).95,length=body_motion.trajectories(i)(1).length,height=scale*length,radius=(T).18;
        RIGID_BODY<TV>& rigid_body=solids_tests.Add_Analytic_Cylinder(height,radius);
        int id=rigid_body.particle_index;assert(id);
        if(id_to_index.Size()<id) id_to_index.Resize(id);
        id_to_index(id)=i;rigid_body_ids(i)=id;
        T volume=(T)pi*radius*radius*height;
        rigid_body.Set_Mass(density*volume);
        rigid_body.Update_Bounding_Box();
        if(i==1){
            T_INERTIA_TENSOR& inertia_tensor=arb.rigid_body_collection.rigid_body_particle.inertia_tensor(rigid_body_ids(1));
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
    TETRAHEDRALIZED_VOLUME<T>* human=&solids_tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/body_150_scaled_10.tet",human_state,false,false,1000,(T)1);

    // binding the deformable particles to the rigid bodies
    for(int p=0;p<rigid_body_collection.rigid_body_particle.array_collection->Size();p++) solids_tests.Bind_Particles_In_Rigid_Body(rigid_body_collection.Rigid_Body(p));

    arb.Use_PD_Actuators();
    arb.global_post_stabilization=true;
    arb.poststabilization_projection_iterations=5;
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=5;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;

    human->Update_Number_Nodes();human->Initialize_Triangulated_Surface();
    TRIANGULATED_SURFACE<T>& triangulated_surface=*human->triangulated_surface;
    triangulated_surface.Update_Triangle_List();triangulated_surface.Initialize_Hierarchy();
    human_levelset=solids_tests.Read_Or_Initialize_Implicit_Surface(STRING_UTILITIES::string_sprintf("%s/human_undeformed_levelset.phi",output_directory.c_str()),triangulated_surface);
}
//#####################################################################
// Function Floppy_Fish
//#####################################################################
void Floppy_Fish()
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
    solids_parameters.enforce_poststabilization_in_cg=false;
    solids_parameters.use_post_cg_constraints=false;
    arb.Set_Iterative_Tolerance((T)1e-4);
    arb.Set_Contact_Level_Iterations(5);
    arb.Set_Use_Shock_Propagation(false);
    arb.Set_Do_Final_Pass(false);
    arb.Set_Poststabilization_Iterations(5);
    arb.Use_PD_Actuators();
    arb.global_post_stabilization=true;

    T friction=(T).2;
    T scale=(T)2;
    T bone_density=200;
    T bone_unscaled_volume=.1875;
    RIGID_BODY_STATE<TV> fish_state(FRAME<TV>(TV(0,(T)3,0),ROTATION<TV>((T)half_pi,TV(1,0,0))));

    ARRAY<RIGID_BODY<TV>*> bones;
    T bone_scales[5]={(T)1,(T)1,(T)1,(T).7,(T).6};
    TV bone_positions[5]={TV(-scale*(T)2,(T)3,0),TV(-scale*(T).75,(T)3,0),TV(scale*(T).5,(T)3,0),TV(scale*(T)1.7,(T)3,0),TV(scale*(T)2.5,(T)3,0)};
    for(int i=0;i<5;i++){
        T bone_scale=scale*bone_scales[i];
        RIGID_BODY<TV>& bone=solids_tests.Add_Rigid_Body("miniplank25wide2",bone_scale,friction);
        bones.Append(&bone);
        bone.Set_Frame(FRAME<TV>(bone_positions[i]));
        bone.Set_Mass(bone_density*bone_unscaled_volume*std::pow(bone_scale,TV::dimension));}

    T joint_strengths[4]={(T)500,(T)500,(T)200,(T)100};
    for(int i=2;i<=bones.m;i++){
        JOINT<TV>* joint=new POINT_JOINT<TV>;
        Initialize_Joint_Between(joint,*bones(i-1),*bones(i),TV(0,0,1));
        JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
        joint_function->Set_k_p(joint_strengths[i-2]);joint_function->Set_Target_Angle(ROTATION<TV>());}

    // add the fish
    T flesh_density=(T)200;
    TETRAHEDRALIZED_VOLUME<T>* fish=0;
    if(!sub_test){
        BOX<TV> box(-TV((T)3.93278,(T)1.07277,(T)0.384066),TV((T)2.68344,(T)1.1747,(T)0.384353));box*=scale;
        VECTOR<int,3> counts(20,15,5);
        GRID<TV> fish_mattress_grid=GRID<TV>(counts,box);
        solids_tests.Create_Mattress(fish_mattress_grid,true,&fish_state,flesh_density);
        fish_bounding_box=ORIENTED_BOX<TV>(box,fish_state.frame);}
    else{
        fish=&solids_tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/fish_42K.tet",fish_state,true,false,flesh_density,scale);}

    // binding the deformable particles to the rigid bodies
    for(int p=0;p<rigid_body_collection.rigid_body_particle.array_collection->Size();p++) solids_tests.Bind_Particles_In_Rigid_Body(rigid_body_collection.Rigid_Body(p));

    if(fish){
        fish->Update_Number_Nodes();fish->Initialize_Triangulated_Surface();
        TRIANGULATED_SURFACE<T>& triangulated_surface=*fish->triangulated_surface;
        triangulated_surface.Update_Triangle_List();triangulated_surface.Initialize_Hierarchy();
        fish_levelset=solids_tests.Read_Or_Initialize_Implicit_Surface(STRING_UTILITIES::string_sprintf("%s/fish_undeformed_levelset.phi",output_directory.c_str()),triangulated_surface);}
}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    fluids_parameters.Use_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initial_Phi
//#####################################################################
T Initial_Phi(const TV& X) const
{
    T phi=(T)1;
    switch(test_number){
        case 1:{
            T phi_human=(*human_levelset)(X);
            phi=max(X.y-initial_fluid_height,-phi_human);
            break;}
        case 9:{
            T phi_fish=sub_test?(*fish_levelset)(X):fish_bounding_box.Signed_Distance(X);
            phi=max(X.y-initial_fluid_height,-phi_fish);
            break;}
        default:
            phi=water_tests.Initial_Phi(X);}

    for(int i=1;i<=rigid_bodies_to_simulate.m;++i)
        phi=max(phi,-solid_body_collection.rigid_body_collection.Rigid_Body(rigid_bodies_to_simulate(i)).Implicit_Geometry_Extended_Value(X));
    return phi;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    T_GRID& grid=*fluids_parameters.grid;
    ARRAY<T,VECTOR<int,3> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) phi(iterator.Cell_Index())=Initial_Phi(iterator.Location());
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() PHYSBAM_OVERRIDE
{
    for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) 
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=water_tests.Initial_Velocity(iterator.Location())[iterator.Axis()];
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    // Add rigid bodies and initialize deformable objects into solids_tests
    switch(test_number){
        case 1: Floppy_Human();break;
        case 9: Floppy_Fish();break;
        default: LOG::cerr<<"Missing implementation for test number "<<test_number<<std::endl;exit(1);}

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    // Add deformable object collision structures and forces
    switch(test_number){
        case 1:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            T stiffness=(T)2e5,damping=(T).03;
            solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>(stiffness,(T).45,damping,(T).25),true,(T).1));
            deformable_objects_to_simulate.Append(new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(tetrahedralized_volume.Get_Boundary_Object()));
            break;}
        case 9:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            T stiffness=(T)2e5,damping=(T).03;
            solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>(stiffness,(T).45,damping,(T).25),true,(T).1));
            deformable_objects_to_simulate.Append(new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(tetrahedralized_volume.Get_Boundary_Object()));
            break;}
        default:break;
    }

    // Add everything to the simulation
    solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true));
    for(int i=1;i<=deformable_objects_to_simulate.m;++i){
        DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& collision_structure=*deformable_objects_to_simulate(i);
        collision_structure.object.Initialize_Hierarchy();
        Add_To_Fluid_Simulation(collision_structure);}
    for(int i=1;i<=rigid_bodies_to_simulate.m;++i){
        RIGID_BODY<TV>& rigid_body_to_add=rigid_body_collection.Rigid_Body(rigid_bodies_to_simulate(i));
        if(rigid_body_to_add.thin_shell) Add_Thin_Shell_To_Fluid_Simulation(rigid_body_to_add); else Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_to_add);}

    // add a ground and add it to the collision body list
    RIGID_BODY<TV>* ground=0;
    switch(test_number){
        case 10:
        case 21: ground=&solids_tests.Add_Ground((T).4,(T)0,(T)1);break;
        case 12: ground=&solids_tests.Add_Ground((T).4,(T)0,(T).15);break;
        default: ground=&solids_tests.Add_Ground();break;}
    if(ground){
        if(test_number==9) rigid_bodies_to_collide_against.Append(ground->particle_index);}

    // collide structures with the ground and walls only
    if(test_number==1 || test_number==9){
        deformable_body_collection.collisions.Use_Structure_Collide_Collision_Body(true);
        for(int s=0;s<deformable_body_collection.deformable_geometry.structures.m;s++) for(int r=0;r<rigid_bodies_to_collide_against.m;r++)
            deformable_body_collection.collisions.structure_collide_collision_body(s).Set(rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(rigid_bodies_to_collide_against(r)));
        for(int i=0;i<solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->limit_time_step_by_strain_rate=false;}
}
//#####################################################################
};
}
#endif

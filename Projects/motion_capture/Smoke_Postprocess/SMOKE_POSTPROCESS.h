//#####################################################################
// Copyright 2009, Michael Lentine
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SMOKE_POSTPROCESS
//#####################################################################
//   1. Add Smoke to an existing rigid body simulation
//#####################################################################
#ifndef __SMOKE_POSTPROCESS__
#define __SMOKE_POSTPROCESS__

#include <PhysBAM_Tools/Arrays/IDENTITY_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/SMOKE_STANDARD_TESTS_3D.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class SMOKE_POSTPROCESS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    
    typedef T_input T;typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
    typedef typename GRID<TV>::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;using BASE::output_directory;using BASE::last_frame;using BASE::frame_rate;
    using BASE::data_directory;using BASE::stream_type;using BASE::solid_body_collection;using BASE::parse_args;using BASE::test_number;using BASE::resolution;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions; // silence -Woverloaded-virtual
    using BASE::Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization;using BASE::Add_Volumetric_Body_To_Fluid_Simulation;using BASE::Add_To_Fluid_Simulation;

    SMOKE_STANDARD_TESTS_3D<GRID<TV> > smoke_tests;
    SOLIDS_STANDARD_TESTS<TV> solids_tests;
    GRID<TV> mattress_grid;
    ARRAY<int> deformable_object_enslaved_nodes;
    int rigid_body_id;
    T velocity_multiplier;
    T aspect_ratio,side_length;
    T stiffness_multiplier,damping_multiplier,bending_stiffness_multiplier,bending_damping_multiplier;
    BOX<TV> source_box,smoke_box;
    TV source_velocity;
    RIGID_BODY_COLLECTION<TV> *rbl_current,*rbl_next;
    int list_frame,rigid_body_start_frame;
    std::string input_directory;

    DEFORMABLE_BODY_COLLECTION<TV> *def_current,*def_next;
    ARRAY<PAIR<int,TV> > constrained_node_positions;
    ARRAY<int> rigid_bodies_to_simulate;

    SMOKE_POSTPROCESS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.SMOKE),
        smoke_tests(*this,fluids_parameters,fluid_collection.incompressible_fluid_collection,solid_body_collection.rigid_body_collection),
        solids_tests(*this,solid_body_collection),rigid_body_id(0),velocity_multiplier((T)1),aspect_ratio((T)1.1),side_length((T).5),
        stiffness_multiplier((T)1),damping_multiplier((T)1),bending_stiffness_multiplier((T)1),bending_damping_multiplier((T)1),source_velocity((T)0,(T).5,(T)0),
        rbl_current(0),rbl_next(0),def_current(0),def_next(0)
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
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {/*this->dt=dt;*/}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_String_Argument("-input_directory","Specify the input directory");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    smoke_tests.Initialize(Smoke_Test_Number(test_number),resolution);
    LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
    last_frame=1045;
    list_frame=-1;

    if(parse_args->Is_Value_Set("-input_directory")) input_directory=parse_args->Get_String_Value("-input_directory");

    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;

    *fluids_parameters.grid=smoke_tests.grid;
    fluids_parameters.fluid_affects_solid=false;fluids_parameters.solid_affects_fluid=true;
    fluids_parameters.incompressible_iterations=100;
    fluids_parameters.domain_walls[2][2]=fluids_parameters.domain_walls[2][1]=fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][2]=true;
    fluids_parameters.domain_walls[3][1]=fluids_parameters.domain_walls[3][2]=false;
    fluids_parameters.use_vorticity_confinement=true;
    fluids_parameters.confinement_parameter=.25;
    //fluids_parameters.gravity=.98;
    //fluids_parameters.gravity_direction=TV(0,0,1);

    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solids_parameters.rigid_body_evolution_parameters.maximum_rigid_body_time_step_fraction=(T).1;
    solids_parameters.rigid_body_evolution_parameters.rigid_geometry_evolution_parameters.use_kinematic_keyframes=false;

    switch(test_number){
        case 1:
            last_frame=1000;
            rigid_body_start_frame=425;
            solids_parameters.implicit_solve_parameters.cg_iterations=400;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            fluids_parameters.density=(T)1;
            velocity_multiplier=(T)15;
            //fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[1][2]=true;fluids_parameters.domain_walls[2][1]=true;fluids_parameters.domain_walls[3][1]=true;fluids_parameters.domain_walls[3][2]=true;
            source_velocity=TV((T)0,(T)0,(T)-5);
            //source_box=BOX<TV>((T)-10,(T)10,(T)-10,(T)10,(T)42,(T)45);
            source_box=BOX<TV>((T)-3,(T)3,(T)-2,(T)4,(T)42,(T)45);
            smoke_box=BOX<TV>((T)-3,(T)3,(T)-2,(T)4,(T)42,(T)45);
            break;
        case 2:
            last_frame=600;
            fluids_parameters.grid->Initialize(10*resolution+1,15*resolution+1,10*resolution+1,-5.,5.,0,15.,-5.,5.);
            rigid_body_start_frame=50;
            solids_parameters.implicit_solve_parameters.cg_iterations=400;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            fluids_parameters.density=(T)1;
            velocity_multiplier=(T)15;
            fluids_parameters.domain_walls[2][2]=false;
            fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[1][2]=true;fluids_parameters.domain_walls[2][1]=true;fluids_parameters.domain_walls[3][1]=true;fluids_parameters.domain_walls[3][2]=true;
            smoke_tests.world_to_source=MATRIX<T,4>::Identity_Matrix();
            source_velocity=TV((T)0,(T)5,(T)0);
            source_box=BOX<TV>((T)-1,(T)1,(T)0,(T)1,(T)-1,(T)1);
            smoke_box=source_box;
            break;
        case 3:
            last_frame=2000;
            fluids_parameters.grid->Initialize(10*resolution+1,15*resolution+1,10*resolution+1,-5.,5.,0,15.,-5.,5.);
            rigid_body_start_frame=50;
            solids_parameters.implicit_solve_parameters.cg_iterations=400;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            fluids_parameters.density=(T)1;
            velocity_multiplier=(T)15;
            fluids_parameters.domain_walls[2][2]=false;
            fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[1][2]=true;fluids_parameters.domain_walls[2][1]=true;fluids_parameters.domain_walls[3][1]=true;fluids_parameters.domain_walls[3][2]=true;
            smoke_tests.world_to_source=MATRIX<T,4>::Identity_Matrix();
            source_velocity=TV((T)0,(T)5,(T)0);
            source_box=BOX<TV>((T)-1,(T)1,(T)0,(T)1,(T)-1,(T)1);
            smoke_box=source_box;
            break;
        default:
            LOG::cerr<<"Unrecognized test number "<<test_number<<std::endl;exit(1);}
        
    output_directory=STRING_UTILITIES::string_sprintf("Smoke_Postprocess/Test_%d_Resolution_%d",test_number,resolution);
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
RIGID_BODY_COLLECTION<TV>* Create_Rigid_Object(int frame)
{
    RIGID_BODY_COLLECTION<TV> *rigid_body_collection=new RIGID_BODY_COLLECTION<TV>(0,0);ARRAY<int> needs_init;
    rigid_body_collection->Read(stream_type,input_directory,frame,&needs_init);
    return rigid_body_collection;
}
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    T threshold=0;
    int current_frame=(int)(time*frame_rate+threshold),next_frame=current_frame+1;
    T alpha=time*frame_rate-current_frame;
    if(next_frame-rigid_body_start_frame>list_frame){
        if(rbl_current) delete &rbl_current;
        rbl_current=rbl_next;rbl_next=Create_Rigid_Object(next_frame-rigid_body_start_frame);
        list_frame=next_frame-rigid_body_start_frame;}
    if(!rbl_current) frame=rbl_next->Rigid_Body(id).Frame();
    else frame=FRAME<TV>::Interpolation(rbl_current->Rigid_Body(id).Frame(),rbl_next->Rigid_Body(id).Frame(),alpha);
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
DEFORMABLE_BODY_COLLECTION<TV>* Create_Deformable_Body_List(int frame)
{
    DEFORMABLE_BODY_COLLECTION<TV> *deformable_body_collection=new DEFORMABLE_BODY_COLLECTION<TV>(solid_body_collection.example_forces_and_velocities,(*new COLLISION_GEOMETRY_COLLECTION<TV>));
    bool include_static_frame=frame==0;
    deformable_body_collection->Read(stream_type,input_directory,frame,-1,include_static_frame,solids_parameters.write_from_every_process);
    return deformable_body_collection;
}
void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE
{
    if(test_number!=2) return;
    T threshold=0;
    int current_frame=(int)(time*frame_rate+threshold),next_frame=current_frame+1;
    T alpha=time*frame_rate-current_frame;
    if(next_frame-rigid_body_start_frame>list_frame){
        LOG::cout<<"Rereading def_current, next_frame="<<next_frame<<" rigid_body_start_frame="<<rigid_body_start_frame<<", list_frame="<<list_frame<<std::endl;
        if(def_current) delete def_current;
        def_current=def_next;def_next=Create_Deformable_Body_List(next_frame-rigid_body_start_frame);
        list_frame=next_frame-rigid_body_start_frame;}

    if(def_current) LOG::cout<<"def_current, so using that for position and velocity updates with alpha="<<alpha<<std::endl;
    else LOG::cout<<"NO def_current, so cannot use that for position and velocity updates"<<std::endl;

    if(!def_current) for(int i=1;i<=def_next->particles.array_collection->Size();i++) X(i)=def_next->particles.X(i);
    else for(int i=1;i<=def_next->particles.array_collection->Size();i++) X(i)=def_current->particles.X(i)*(1-alpha)+def_next->particles.X(i)*alpha;
}
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Adjust_Density_And_Temperature_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    BASE::Adjust_Density_And_Temperature_With_Sources(smoke_box,smoke_tests.world_to_source,smoke_tests.rho,fluids_parameters.temperature_products);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{
    T threshold=0;
    int current_frame=(int)(time*frame_rate+threshold),next_frame=current_frame+1;
    T alpha=time*frame_rate-current_frame;
    if(next_frame-rigid_body_start_frame>list_frame){
        if(rbl_current) delete &rbl_current;
        rbl_current=rbl_next;rbl_next=Create_Rigid_Object(next_frame-rigid_body_start_frame);
        list_frame=next_frame-rigid_body_start_frame;}
    if(!rbl_current) twist=rbl_next->rigid_body_particle.twist(id);
    else{
        twist.linear=(rbl_current->rigid_body_particle.twist(id).linear*(1-alpha)+rbl_next->rigid_body_particle.twist(id).linear*alpha);
        twist.angular=(rbl_current->rigid_body_particle.twist(id).angular*(1-alpha)+rbl_next->rigid_body_particle.twist(id).angular*alpha);}
    return true;
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(test_number!=2) return;
    T threshold=0;
    int current_frame=(int)(velocity_time*frame_rate+threshold),next_frame=current_frame+1;
    T alpha=velocity_time*frame_rate-current_frame;
    if(next_frame-rigid_body_start_frame>list_frame){
        LOG::cout<<"Rereading def_current, next_frame="<<next_frame<<" rigid_body_start_frame="<<rigid_body_start_frame<<", list_frame="<<list_frame<<std::endl;
        if(def_current) delete def_current;
        def_current=def_next;def_next=Create_Deformable_Body_List(next_frame-rigid_body_start_frame);
        list_frame=next_frame-rigid_body_start_frame;}
    if(!def_current) for(int i=1;i<=def_next->particles.array_collection->Size();i++) V(i)=def_next->particles.V(i);
    else for(int i=1;i<=def_next->particles.array_collection->Size();i++) V(i)=def_current->particles.V(i)*(1-alpha)+def_next->particles.V(i)*alpha;
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) PHYSBAM_OVERRIDE
{}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(test_number!=2) return;
    for(int i=1;i<=def_next->particles.array_collection->Size();i++) V(i)=TV();
}
//#####################################################################
// Function Smoke_Test_Number
//#####################################################################
static int Smoke_Test_Number(const int test_number)
{
    switch(test_number){
        case 1:
        case 2:
        case 3:
            return 1;
        default:
            return 1;}
}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    fluids_parameters.Use_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() PHYSBAM_OVERRIDE
{
    for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) 
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=velocity_multiplier*smoke_tests.Initial_Velocity(iterator.Location())[iterator.Axis()];
    BASE::Initialize_Velocities();
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<3> >& face_velocities,ARRAY<bool,FACE_INDEX<3> >& psi_N,const T time) PHYSBAM_OVERRIDE
{
    BASE::Get_Source_Velocities(source_box,smoke_tests.world_to_source,source_velocity);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    if(test_number==1 || test_number==3){
        rbl_next=Create_Rigid_Object(0);list_frame=0;ARRAY<int> needs_init;
        rigid_body_collection.Read(stream_type,input_directory,0,&needs_init);}
    for(int id(1);id<rigid_body_collection.rigid_body_particle.array_collection->Size();id++) if(!rigid_body_collection.Rigid_Body(id).Has_Infinite_Inertia()) rigid_body_collection.Rigid_Body(id).Is_Kinematic()=true;

    if(test_number==2){
        def_next=Create_Deformable_Body_List(0);list_frame=0;
        deformable_body_collection.Read(stream_type,input_directory,0,-1,true,solids_parameters.write_from_every_process);
        TRIANGULATED_SURFACE<T>* surface=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>*>(1);

        DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>* collision_structure=new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(*surface);
        collision_structure->object.Initialize_Hierarchy();
        Add_To_Fluid_Simulation(*collision_structure);
    }

    // add structures and rigid bodies to collisions
    //deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    //solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    //solids_parameters.collision_body_list.Add_Bodies(solid_body_collection.rigid_body_collection);

    // correct number nodes
    for(int i=1;i<=deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();
    if(test_number==1 || test_number==3){
        for(int id(1);id<rigid_body_collection.rigid_body_particle.array_collection->Size();id++)
            if(rigid_body_collection.Is_Active(id)){
                Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(id));}
        //DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(triangulated_surface);
        //deformable_collisions.object.Initialize_Hierarchy();
        //Add_To_Fluid_Simulation(deformable_collisions,true,false);
    }
    //THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Rigid_Body_Walls(*this);
}
//#####################################################################
};
}
#endif

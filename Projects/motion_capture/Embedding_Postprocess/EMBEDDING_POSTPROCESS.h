//#####################################################################
// Copyright 2009, Michael Lentine
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDING_POSTPROCESS
//#####################################################################
//   1. Add Smoke to an existing rigid body simulation
//#####################################################################
#ifndef __EMBEDDING_POSTPROCESS__
#define __EMBEDDING_POSTPROCESS__

#include <PhysBAM_Tools/Arrays/IDENTITY_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
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
class EMBEDDING_POSTPROCESS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    
    typedef T_input T;typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
    typedef typename GRID<TV>::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::output_directory;using BASE::last_frame;using BASE::frame_rate;
    using BASE::data_directory;using BASE::stream_type;using BASE::solid_body_collection;using BASE::parse_args;using BASE::test_number;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions; // silence -Woverloaded-virtual
    using BASE::Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization;using BASE::Add_Volumetric_Body_To_Fluid_Simulation;

    int read_test_number;
    SOLIDS_STANDARD_TESTS<TV> solids_tests;
    GRID<TV> mattress_grid;
    ARRAY<int> deformable_object_enslaved_nodes;
    int rigid_body_id;
    T velocity_multiplier;
    T aspect_ratio,side_length;
    T stiffness_multiplier,damping_multiplier,bending_stiffness_multiplier,bending_damping_multiplier;
    BOX<TV> source_box,smoke_box;
    TV source_velocity;
    DEFORMABLE_BODY_COLLECTION<TV> *rbl_current,*rbl_next,*rbl_prev;
    int list_frame,rigid_body_start_frame;
    std::string input_directory;
    T dt;

    ARRAY<PAIR<int,TV> > constrained_node_positions;
    ARRAY<int> rigid_bodies_to_simulate;

    EMBEDDING_POSTPROCESS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),
        solids_tests(*this,solid_body_collection),rigid_body_id(0),velocity_multiplier((T)1),aspect_ratio((T)1.1),side_length((T).5),
        stiffness_multiplier((T)1),damping_multiplier((T)1),bending_stiffness_multiplier((T)1),bending_damping_multiplier((T)1),source_velocity((T)0,(T).5,(T)0),
        rbl_current(0),rbl_next(0),rbl_prev(0),dt(0)
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
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {this->dt=dt;}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE {}
    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return false;}
    void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
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
    LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
    last_frame=1045;
    list_frame=-1;
    rigid_body_start_frame=0;
    read_test_number=0;

    if(parse_args->Is_Value_Set("-input_directory")) input_directory=parse_args->Get_String_Value("-input_directory");

    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;

    frame_rate=60;
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solids_parameters.rigid_body_evolution_parameters.maximum_rigid_body_time_step_fraction=(T).1;
    solids_parameters.rigid_body_evolution_parameters.rigid_geometry_evolution_parameters.use_kinematic_keyframes=false;

    switch(test_number){
        case 1:
            last_frame=1000;
            solids_parameters.implicit_solve_parameters.cg_iterations=400;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            break;
        case 2:
            last_frame=2000;
            read_test_number=6;
            solids_parameters.implicit_solve_parameters.cg_iterations=400;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            break;
        default:
            LOG::cerr<<"Unrecognized test number "<<test_number<<std::endl;exit(1);}
        
    output_directory=STRING_UTILITIES::string_sprintf("Embedding_Postprocess/Test_%d",test_number);
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
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
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE
{
    T threshold=0;
    int current_frame=(int)(time*frame_rate+threshold),next_frame=current_frame+1;
    T alpha=time*frame_rate-current_frame;
    if(next_frame-rigid_body_start_frame>list_frame){
        LOG::cout<<"Rereading rbl_current, next_frame="<<next_frame<<" rigid_body_start_frame="<<rigid_body_start_frame<<", list_frame="<<list_frame<<std::endl;
        if(rbl_current) delete rbl_current;
        rbl_current=rbl_next;rbl_next=Create_Deformable_Body_List(next_frame-rigid_body_start_frame);
        list_frame=next_frame-rigid_body_start_frame;}

    if(rbl_current) LOG::cout<<"rbl_current, so using that for position and velocity updates with alpha="<<alpha<<std::endl;
    else LOG::cout<<"NO rbl_current, so cannot use that for position and velocity updates"<<std::endl;

    if(!rbl_current) for(int i=1;i<=rbl_next->particles.array_collection->Size();i++) X(i)=rbl_next->particles.X(i);
    else for(int i=1;i<=rbl_next->particles.array_collection->Size();i++) X(i)=rbl_current->particles.X(i)*(1-alpha)+rbl_next->particles.X(i)*alpha;
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    T threshold=0;
    int current_frame=(int)(velocity_time*frame_rate+threshold),next_frame=current_frame+1;
    T alpha=velocity_time*frame_rate-current_frame;
    if(next_frame-rigid_body_start_frame>list_frame){
        LOG::cout<<"Rereading rbl_current, next_frame="<<next_frame<<" rigid_body_start_frame="<<rigid_body_start_frame<<", list_frame="<<list_frame<<std::endl;
        if(rbl_current) delete rbl_current;
        rbl_current=rbl_next;rbl_next=Create_Deformable_Body_List(next_frame-rigid_body_start_frame);
        list_frame=next_frame-rigid_body_start_frame;}
    if(!rbl_current) for(int i=1;i<=rbl_next->particles.array_collection->Size();i++) V(i)=rbl_next->particles.V(i);
    else for(int i=1;i<=rbl_next->particles.array_collection->Size();i++) V(i)=rbl_current->particles.V(i)*(1-alpha)+rbl_next->particles.V(i)*alpha;
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    for(int i=1;i<=rbl_next->particles.array_collection->Size();i++) V(i)=TV();
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;

    rbl_next=Create_Deformable_Body_List(0);list_frame=0;
    deformable_body_collection.Read(stream_type,input_directory,0,-1,true,solids_parameters.write_from_every_process);

    TETRAHEDRALIZED_VOLUME<T>* volume=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>(1);
    TRIANGULATED_SURFACE<T>* surface_original=TRIANGULATED_SURFACE<T>::Create(),*surface=0;
    if(FILE_UTILITIES::File_Exists(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tri",read_test_number))){
        FILE_UTILITIES::Read_From_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tri",read_test_number),*surface_original);
        if(read_test_number==6) for(int i=1;i<=surface_original->particles.array_collection->Size();++i) surface_original->particles.X(i)+=TV(0,(T)7,0);
        surface=(TRIANGULATED_SURFACE<T>*)surface_original->Append_Particles_And_Create_Copy(particles);}
    else PHYSBAM_FATAL_ERROR();
    surface->Update_Number_Nodes();
    deformable_body_collection.deformable_geometry.Add_Structure(surface);

    volume->Update_Number_Nodes();
    volume->Initialize_Hierarchy();
    SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(*volume,1,true);

    ARRAY<TRIPLE<int,int,TV> > bindings; 
    ARRAY<int> tets;const T tolerance=(T)1e-4;
    ARRAY<int> surface_particles;surface->mesh.elements.Flattened().Get_Unique(surface_particles);
    for(int i=0;i<surface_particles.m;i++){int p=surface_particles(i);
        tets.Remove_All();volume->hierarchy->Intersection_List(particles.X(p),tets,tolerance);T tol_current=tolerance;
        while(tets.m==0){tol_current*=2;volume->hierarchy->Intersection_List(particles.X(p),tets,tol_current);}
        TV bary=TETRAHEDRON<T>::First_Three_Barycentric_Coordinates(particles.X(p),particles.X.Subset(volume->mesh.elements(tets(1))));
        bindings.Append(TRIPLE<int,int,TV>(p,tets(1),bary));}
    for(int i=0;i<bindings.m;i++){
        if(bindings(i).y==0) continue;
        VECTOR<int,4> nodes=volume->mesh.elements(bindings(i).y);
        binding_list.Add_Binding(new LINEAR_BINDING<TV,4>(particles,bindings(i).x,nodes,bindings(i).z));}
    surface->Update_Number_Nodes();
    volume->Update_Number_Nodes();

    // add structures and rigid bodies to collisions
    //deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    //solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    //solids_parameters.collision_body_list.Add_Bodies(solid_body_collection.deformable_body_collection.rigid_body_particles);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();
}
//#####################################################################
};
}
#endif

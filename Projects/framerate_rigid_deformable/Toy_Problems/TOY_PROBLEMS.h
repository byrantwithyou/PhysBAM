//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TOY_PROBLEMS
//#####################################################################
// 1. Particle falling (make sure it falls at the right rate).
// 2. Oscillating spring w/o gravity (make sure the frequency is right and no energy loss)
// 3. Spring falling on the ground (handle collisions without energy loss)
// 4. Rotating spring
// 5. Contrained spring (make sure the frequency is right and no energy loss with constraints)
// 6. Two spring network
// 7. 10 spring network
//#####################################################################
#ifndef __TOY_PROBLEMS__
#define __TOY_PROBLEMS__

#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>

namespace PhysBAM{

template<class T_input>
class TOY_PROBLEMS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;
    typedef VECTOR<T_input,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    SOLIDS_STANDARD_TESTS<TV> tests;

    T stiffness_multiplier;
    T damping_multiplier;
    bool use_be,use_tr;
    
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::restart;
    using BASE::solid_body_collection;using BASE::solids_evolution;using BASE::parse_args;using BASE::test_number;using BASE::frame_rate;

    TOY_PROBLEMS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),use_be(true),use_tr(false)
    {
    }

    ~TOY_PROBLEMS()
    {}

    // Unused callbacks
    void Post_Initialization() PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > F,const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}
    void Advance_One_Time_Step_End_Callback(const T dt,const T time) {}
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) {}
    void Self_Collisions_Begin_Callback(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Set_Rigid_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Set_Deformable_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return true;}
    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE {}
    
//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Double_Argument("-stiffen",1,"","stiffness multiplier for various tests");
    parse_args->Add_Double_Argument("-dampen",1,"","damping multiplier for various tests");
    parse_args->Add_Option_Argument("-print_energy","print energy statistics");
    parse_args->Add_Option_Argument("-use_be","use backward Euler");
    parse_args->Add_Option_Argument("-use_tr","use trapezoid rule");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    output_directory=STRING_UTILITIES::string_sprintf("Toy_Problems/Test_%d",test_number);
    stiffness_multiplier=(T)parse_args->Get_Double_Value("-stiffen");
    damping_multiplier=(T)parse_args->Get_Double_Value("-dampen");
    solid_body_collection.print_energy=parse_args->Get_Option_Value("-print_energy");
    if(parse_args->Is_Value_Set("-use_tr")){use_tr=true;use_be=false;}
    if(parse_args->Is_Value_Set("-use_be")){use_tr=false;use_be=true;}
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    switch(test_number){
        case 1: Falling_Particle_Test();break;
        case 2: Oscillating_Spring_No_Gravity_Test();break;
        case 3: Spring_Collision_Test();break;
        case 4: Rotating_Spring_Test();break;
        case 5: Constrained_Spring_Test();break;
        case 6: Two_Spring_Test();break;
        case 7: Ten_Spring_Test();break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    if(use_tr) solids_parameters.use_trapezoidal_rule_for_velocities=true;
    else if(use_be) solids_parameters.use_trapezoidal_rule_for_velocities=false;
    deformable_body_collection.particles.Compute_Auxiliary_Attributes(soft_bindings);
    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // add forces
    switch(test_number){
        case 1:
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            break;
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:{
            SEGMENTED_CURVE<TV>& segmented_curve=deformable_body_collection.deformable_geometry.template Find_Structure<SEGMENTED_CURVE<TV>&>();
            T linear_stiffness=stiffness_multiplier*10/(1+sqrt((T)2)),linear_damping=damping_multiplier*15;
            LINEAR_SPRINGS<TV>* spring_force=Create_Edge_Springs(segmented_curve,linear_stiffness,linear_damping);
            solid_body_collection.Add_Force(spring_force);
            if(test_number==2 || test_number==3)
                solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            ARRAY<T> restlengths(spring_force->segment_mesh.elements.m);restlengths.Fill((T).75);
            spring_force->Set_Restlength(restlengths);
            break;}
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
}
//#####################################################################
// Function Falling_Particle_Test
//#####################################################################
void Falling_Particle_Test()
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    FREE_PARTICLES<TV>& free_particles=*FREE_PARTICLES<TV>::Create();solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&free_particles);
    int new_particle=particles.array_collection->Add_Element();
    particles.mass(new_particle)=1;particles.X(new_particle)=TV(0,1,0);
    free_particles.nodes.Append(new_particle);
}
//#####################################################################
// Function Oscillating_Spring_No_Gravity_Test
//#####################################################################
void Oscillating_Spring_No_Gravity_Test()
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    SEGMENTED_CURVE<TV>& segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve);
    int new_edge_node1=segmented_curve.particles.array_collection->Add_Element();int new_edge_node2=segmented_curve.particles.array_collection->Add_Element();
    static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node1)=static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node2)=1;
    segmented_curve.particles.X(new_edge_node1)=TV(0,(T)1,(T)-.5);segmented_curve.particles.X(new_edge_node2)=TV(0,(T)1,(T).5);
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node1,new_edge_node2));
}
//#####################################################################
// Function Spring_Collision_Test
//#####################################################################
void Spring_Collision_Test()
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    SEGMENTED_CURVE<TV>& segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve);
    int new_edge_node1=segmented_curve.particles.array_collection->Add_Element();int new_edge_node2=segmented_curve.particles.array_collection->Add_Element();
    static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node1)=static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node2)=1;
    segmented_curve.particles.X(new_edge_node1)=TV(0,(T)1,0);segmented_curve.particles.X(new_edge_node2)=TV(0,(T)2,0);
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node1,new_edge_node2));
    tests.Add_Ground();
}
//#####################################################################
// Function Rotating_Spring_Test
//#####################################################################
void Rotating_Spring_Test()
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    SEGMENTED_CURVE<TV>& segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve);
    int new_edge_node1=segmented_curve.particles.array_collection->Add_Element();int new_edge_node2=segmented_curve.particles.array_collection->Add_Element();
    static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node1)=static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node2)=1;
    segmented_curve.particles.X(new_edge_node1)=TV(0,(T)1,(T)-.5);segmented_curve.particles.X(new_edge_node2)=TV(0,(T)1,(T).5);
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node1,new_edge_node2));
    particles.V(new_edge_node1)=TV((T).5,0,0);particles.V(new_edge_node2)=TV((T)-.5,0,0);
}
//#####################################################################
// Function Constrained_Spring_Test
//#####################################################################
void Constrained_Spring_Test()
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    SEGMENTED_CURVE<TV>& segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve);
    int new_edge_node1=segmented_curve.particles.array_collection->Add_Element();int new_edge_node2=segmented_curve.particles.array_collection->Add_Element();
    static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node1)=static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node2)=1;
    segmented_curve.particles.X(new_edge_node1)=TV(0,(T)1.5,(T)0);segmented_curve.particles.X(new_edge_node2)=TV(0,(T).5,(T)0);
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node1,new_edge_node2));
    last_frame=1000;
}
//#####################################################################
// Function Two_Spring_Test
//#####################################################################
void Two_Spring_Test()
{
    Create_Spring_Chain(2);
}
//#####################################################################
// Function Ten_Spring_Test
//#####################################################################
void Ten_Spring_Test()
{
    Create_Spring_Chain(10);
}
//#####################################################################
// Function Create_Spring_Chain
//#####################################################################
void Create_Spring_Chain(int count)
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    SEGMENTED_CURVE<TV>& segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve);
    T current_endpoint_z=0;
    int current_node=segmented_curve.particles.array_collection->Add_Element();
    static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve.particles).mass(current_node)=1;
    segmented_curve.particles.X(current_node)=TV(0,(T)1,(T)current_endpoint_z);
    for(int i=0;i<count;i++){
        int new_edge_node=segmented_curve.particles.array_collection->Add_Element();
        static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node)=static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve.particles).mass(current_node);
        segmented_curve.particles.X(new_edge_node)=TV(0,(T)1,(T)current_endpoint_z+1);
        segmented_curve.mesh.elements.Append(VECTOR<int,2>(current_node,new_edge_node));
        current_endpoint_z+=1;current_node=new_edge_node;}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(test_number==5) V(1)=TV();
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(test_number==5) V(1)=TV();
}
//#####################################################################
};
}
#endif

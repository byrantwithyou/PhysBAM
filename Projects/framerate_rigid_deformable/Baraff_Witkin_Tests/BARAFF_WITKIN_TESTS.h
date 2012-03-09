//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BARAFF_WITKIN_TESTS
//#####################################################################
//  1. Cloth panel
//  2. Permuted triangle test
//  3. Permuted triangle test
//  4. Shear test
//  5. Bending test
//  6. Cloth constrained to table top
//  7. Cloth dropped onto rigid block
//  8. Cloth constrained by two corners and draped on ground
//#####################################################################
#ifndef __BARAFF_WITKIN_TESTS__
#define __BARAFF_WITKIN_TESTS__

#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BW_BENDING_FORCES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BW_GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BW_SHEAR_FORCES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BW_STRETCH_FORCES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/SOLIDS_EXAMPLE.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>

namespace PhysBAM{

template<class T_input>
class BARAFF_WITKIN_TESTS:public SOLIDS_EXAMPLE<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_EXAMPLE<TV> BASE;
    using BASE::solids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;
    using BASE::stream_type;using BASE::write_substeps_level;using BASE::solid_body_collection;using BASE::parse_args;using BASE::test_number;

    SOLIDS_STANDARD_TESTS<TV> tests;

    T cloth_cfl;
    int number_side_panels;
    T aspect_ratio;
    T side_length;
    bool use_shear,use_bend;
    int parameter;

    bool use_bridson_forces;
    T stiffness_multiplier;
    T damping_multiplier;
    T bending_stiffness_multiplier;
    T bending_damping_multiplier;
    T axial_bending_stiffness_multiplier,axial_bending_damping_multiplier;

    BARAFF_WITKIN_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type),tests(*this,solid_body_collection),number_side_panels(31),aspect_ratio((T)1.0),side_length((T)1.0),use_shear(false),use_bend(false),
        use_bridson_forces(true)
    {}

    ~BARAFF_WITKIN_TESTS()
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
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return false;}
    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Double_Argument("-clothcfl",4.,"Cloth CFL");
    parse_args->Add_Option_Argument("-print_energy","print energy statistics");
    parse_args->Add_Option_Argument("-test_system","test arb system properties");
    parse_args->Add_Option_Argument("-print_matrix","print arb system");
    parse_args->Add_Option_Argument("-use_shear","use shear forces");
    parse_args->Add_Option_Argument("-use_bend","use bending forces");
    parse_args->Add_Integer_Argument("-parameter",0,"parameter used by multiple tests to change the parameters of the test");

    parse_args->Add_Double_Argument("-stiffen",1,"","stiffness multiplier for various tests");
    parse_args->Add_Double_Argument("-dampen",1,"","damping multiplier for various tests");
    parse_args->Add_Double_Argument("-stiffen_bending",1,"","stiffness multiplier for bending springs in various cloth tests");
    parse_args->Add_Double_Argument("-dampen_bending",1,"","damping multiplier for bending springs in various cloth tests");
    parse_args->Add_Double_Argument("-dampen_axial_bending",1,"","axial damping multiplier for bending springs in various cloth tests");
    parse_args->Add_Double_Argument("-stiffen_axial_bending",1,"","axial stiffness multiplier for bending springs in various cloth tests");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    output_directory=STRING_UTILITIES::string_sprintf("Baraff_Witkin_Tests/Test_%d",test_number);
    cloth_cfl=(T)parse_args->Get_Double_Value("-clothcfl");
    solid_body_collection.print_energy=parse_args->Get_Option_Value("-print_energy");
    if(parse_args->Is_Value_Set("-test_system")) solids_parameters.implicit_solve_parameters.test_system=true;
    if(parse_args->Is_Value_Set("-print_matrix")) solids_parameters.implicit_solve_parameters.print_matrix=true;
    use_shear=parse_args->Get_Option_Value("-use_shear");
    use_bend=parse_args->Get_Option_Value("-use_bend");
    parameter=parse_args->Get_Integer_Value("-parameter");

    stiffness_multiplier=(T)parse_args->Get_Double_Value("-stiffen");
    damping_multiplier=(T)parse_args->Get_Double_Value("-dampen");
    bending_stiffness_multiplier=(T)parse_args->Get_Double_Value("-stiffen_bending");
    bending_damping_multiplier=(T)parse_args->Get_Double_Value("-dampen_bending");
    axial_bending_stiffness_multiplier=(T)parse_args->Get_Double_Value("-stiffen_axial_bending");
    axial_bending_damping_multiplier=(T)parse_args->Get_Double_Value("-dampen_axial_bending");
}
//#####################################################################
// Function Parse_Late_Options
//#####################################################################
void Parse_Late_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Late_Options();
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    switch(test_number){
        case 1: Curtain();break;
        case 2: Triangle_Permutations();break;
        case 3: Single_Triangle();break;
        case 4: Shear_Test();break;
        case 5: Bending_Test();break;
        case 6: Cloth_Table_Test();break;
        case 7: Cloth_Body_Test();break;
        case 8: Cloth_Draped_On_Ground();break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    deformable_body_collection.particles.Compute_Auxiliary_Attributes(soft_bindings);
    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // add forces
    switch(test_number){
        case 1:
        case 2:
        case 3:
            Triangulated_Surface_Forces(true,false,false);
            for(int i=0;i<deformable_body_collection.particles.array_collection->Size();i++)
                deformable_body_collection.particles.X(i)*=1.2;
            break;
        case 4:{
            Triangulated_Surface_Forces(false,true,false);
            T theta=pi/16;
            deformable_body_collection.particles.X(2).x=deformable_body_collection.particles.X(2).z*sin(theta);
            deformable_body_collection.particles.X(2).z=deformable_body_collection.particles.X(2).z*cos(theta);}
            break;
        case 5:
            Triangulated_Surface_Forces(false,false,true);
            break;
        case 6:
        case 7:
        case 8:
            if(use_bridson_forces){
                T linear_stiffness=stiffness_multiplier*10/(1+sqrt((T)2)),linear_damping=damping_multiplier*15;
                Triangulated_Surface_Forces(true,true,true,linear_stiffness,linear_damping,stiffness_multiplier*8/(1+sqrt((T)2)),damping_multiplier*4);
                solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));}
            else{
                Triangulated_Surface_Forces(true,true,true);
                solid_body_collection.Add_Force(new BW_GRAVITY<TV>(deformable_body_collection.particles,9.8,TV((T)0,(T)-1,(T)0)));}
            break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}
}
//#####################################################################
// Function Curtain
//#####################################################################
void Curtain()
{
    solids_parameters.cfl=(T)50;
    frame_rate=30;
    last_frame=(int)(7*frame_rate);
    solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;
    solids_parameters.cfl=cloth_cfl; // was 4
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-6;
    solids_parameters.implicit_solve_parameters.cg_iterations=200;
    
    tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,0);
}
//#####################################################################
// Function Triangle_Permutations
//#####################################################################
void Triangle_Permutations()
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    solids_parameters.cfl=(T)50;
    frame_rate=30;
    last_frame=(int)(7*frame_rate);
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-6;
    solids_parameters.implicit_solve_parameters.cg_iterations=200;
    
    TRIANGULATED_SURFACE<T>& surface=*TRIANGULATED_SURFACE<T>::Create(particles);solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&surface);
    particles.array_collection->Add_Elements(6);
    particles.mass.Fill((T)1);
    particles.X(0)=TV(0,0,0);particles.X(1)=TV(1,0,0);particles.X(2)=TV(0,0,1);
    particles.X(3)=TV(0,0,0);particles.X(4)=TV(1,0,0);particles.X(5)=TV(0,0,1);
    surface.mesh.elements.Append(VECTOR<int,3>(0,1,2));
    surface.mesh.elements.Append(VECTOR<int,3>(3,4,5));
}
//#####################################################################
// Function Single_Triangle
//#####################################################################
void Single_Triangle()
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    solids_parameters.cfl=(T)50;
    frame_rate=30;
    last_frame=(int)(7*frame_rate);
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-6;
    solids_parameters.implicit_solve_parameters.cg_iterations=200;
    
    TRIANGULATED_SURFACE<T>& surface=*TRIANGULATED_SURFACE<T>::Create(particles);solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&surface);
    particles.array_collection->Add_Elements(3);
    particles.mass.Fill((T)1);
    T theta=0;
    switch(parameter){
        case 0:theta=0;break;
        case 1:theta=pi/4;break;
        case 2:theta=-pi/4;break;}
    particles.X(0)=TV(0,0,0);particles.X(1)=TV(cos(theta),0,-sin(theta));particles.X(2)=TV(sin(theta),0,cos(theta));
    surface.mesh.elements.Append(VECTOR<int,3>(0,1,2));
}
//#####################################################################
// Function Shear_Test
//#####################################################################
void Shear_Test()
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    solids_parameters.cfl=(T)50;
    frame_rate=30;
    last_frame=(int)(7*frame_rate);
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-6;
    solids_parameters.implicit_solve_parameters.cg_iterations=200;
    
    TRIANGULATED_SURFACE<T>& surface=*TRIANGULATED_SURFACE<T>::Create(particles);solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&surface);
    particles.array_collection->Add_Elements(3);
    particles.mass.Fill((T)1);
    particles.X(0)=TV(0,0,0);particles.X(1)=TV(1,0,0);particles.X(2)=TV(0,0,1);
    surface.mesh.elements.Append(VECTOR<int,3>(0,1,2));
}
//#####################################################################
// Function Bending_Test
//#####################################################################
void Bending_Test()
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    solids_parameters.cfl=(T)50;
    frame_rate=30;
    last_frame=(int)(7*frame_rate);
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-6;
    solids_parameters.implicit_solve_parameters.cg_iterations=200;
    
    TRIANGULATED_SURFACE<T>& surface=*TRIANGULATED_SURFACE<T>::Create(particles);solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&surface);
    particles.array_collection->Add_Elements(4);
    particles.mass.Fill((T)1);
    T theta=pi/32;
    particles.X(0)=TV(0,0,0);particles.X(1)=TV(1,0,0);particles.X(2)=TV(0,0,1);particles.X(3)=TV(.5*(1+cos(theta)),sqrt(2)*.5*sin(theta),.5*(1+cos(theta)));
    surface.mesh.elements.Append(VECTOR<int,3>(0,1,2));
    surface.mesh.elements.Append(VECTOR<int,3>(3,2,1));
}
//#####################################################################
// Function Cloth_Body_Test
//#####################################################################
void Cloth_Table_Test()
{
    solids_parameters.cfl=(T)50;
    frame_rate=50;
    last_frame=(int)(7*frame_rate);
    solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;
    solids_parameters.cfl=cloth_cfl; // was 4
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-6;
    solids_parameters.implicit_solve_parameters.cg_iterations=200;
    T density=.1;

    TRIANGULATED_SURFACE<T>& ts=tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,0);
    for(int i=0;i<solid_body_collection.deformable_body_collection.particles.array_collection->Size();i++){
        solid_body_collection.deformable_body_collection.particles.mass(i)=0;}
    for(int i=0;i<ts.mesh.elements.m;i++){
        int node1,node2,node3;ts.mesh.elements(i).Get(node1,node2,node3);
        solid_body_collection.deformable_body_collection.particles.X(node1).y=0;
        solid_body_collection.deformable_body_collection.particles.X(node2).y=0;
        solid_body_collection.deformable_body_collection.particles.X(node3).y=0;
        T area=TRIANGLE_3D<T>(solid_body_collection.deformable_body_collection.particles.X.Subset(VECTOR<int,3>(node1,node2,node3))).Area();
        T m=density*area/(T)3.0;
        solid_body_collection.deformable_body_collection.particles.mass(node1)+=m;
        solid_body_collection.deformable_body_collection.particles.mass(node2)+=m;
        solid_body_collection.deformable_body_collection.particles.mass(node3)+=m;}
    T total_mass=0;
    for(int i=0;i<solid_body_collection.deformable_body_collection.particles.array_collection->Size();i++){
        total_mass+=solid_body_collection.deformable_body_collection.particles.mass(i);}
    LOG::cout << "total mass " << total_mass << std::endl;
}
//#####################################################################
// Function Cloth_Body_Test
//#####################################################################
void Cloth_Body_Test()
{
    solids_parameters.cfl=(T)50;
    frame_rate=30;
    last_frame=(int)(7*frame_rate);
    solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;
    solids_parameters.cfl=cloth_cfl; // was 4
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-6;
    solids_parameters.implicit_solve_parameters.cg_iterations=200;
    T density=.1;

    TRIANGULATED_SURFACE<T>& ts=tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,0);
    for(int i=0;i<ts.mesh.elements.m;i++){
        int node1,node2,node3;ts.mesh.elements(i).Get(node1,node2,node3);
        T area=TRIANGLE_3D<T>(solid_body_collection.deformable_body_collection.particles.X.Subset(VECTOR<int,3>(node1,node2,node3))).Area();
        T m=density*area/(T)3.0;
        solid_body_collection.deformable_body_collection.particles.mass(node1)+=m;
        solid_body_collection.deformable_body_collection.particles.mass(node2)+=m;
        solid_body_collection.deformable_body_collection.particles.mass(node3)+=m;}

    RIGID_BODY<TV>& body=tests.Add_Rigid_Body("box",(T).25,(T)0);
    body.Frame().t=TV((T).5,(T)-.75,(T).5);
    body.is_static=true;
}
//#####################################################################
// Function Cloth_Draped_On_Ground
//#####################################################################
void Cloth_Draped_On_Ground()
{
    solids_parameters.cfl=(T)50;
    frame_rate=30;
    last_frame=(int)(7*frame_rate);
    solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;
    solids_parameters.cfl=cloth_cfl; // was 4
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-6;
    solids_parameters.implicit_solve_parameters.cg_iterations=200;
    T density=.1;

    TRIANGULATED_SURFACE<T>& ts=tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,0);
    for(int i=0;i<ts.mesh.elements.m;i++){
        int node1,node2,node3;ts.mesh.elements(i).Get(node1,node2,node3);
        T area=TRIANGLE_3D<T>(solid_body_collection.deformable_body_collection.particles.X.Subset(VECTOR<int,3>(node1,node2,node3))).Area();
        T m=density*area/(T)3.0;
        solid_body_collection.deformable_body_collection.particles.mass(node1)+=m;
        solid_body_collection.deformable_body_collection.particles.mass(node2)+=m;
        solid_body_collection.deformable_body_collection.particles.mass(node3)+=m;}

//    tests.Add_Ground((T).5,0,1);
}
//#####################################################################
// Function Triangulated_Surface_Forces
//#####################################################################
void Triangulated_Surface_Forces(const bool use_stretch,const bool use_shear,const bool use_bend)
{
    for(int i=0;TRIANGULATED_SURFACE<T>* triangulated_surface=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>*>(i);i++){
        if(use_stretch) solid_body_collection.Add_Force(Create_BW_Stretch_Force(solid_body_collection.deformable_body_collection.particles,triangulated_surface->mesh,(T)1e6,(T)1000));
        if(use_shear) solid_body_collection.Add_Force(Create_BW_Shear_Force(solid_body_collection.deformable_body_collection.particles,triangulated_surface->mesh,(T)100,(T)100));
        if(use_bend) solid_body_collection.Add_Force(Create_BW_Bending_Force(solid_body_collection.deformable_body_collection.particles,triangulated_surface->mesh,(T).01,(T)2e-6));}

    for(int i=0;i<solid_body_collection.deformable_body_collection.deformables_forces.m;i++) 
        solid_body_collection.deformable_body_collection.deformables_forces(i)->use_implicit_velocity_independent_forces=true;
}
//#####################################################################
// Function Triangulated_Surface_Forces
//#####################################################################
void Triangulated_Surface_Forces(const bool use_edge,const bool use_bending,const bool use_altitude,const T linear_stiffness,const T linear_damping,const T altitude_stiffness,
    const T altitude_damping)
{
    for(int i=0;TRIANGULATED_SURFACE<T>* triangulated_surface=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>*>(i);i++){
        if(use_edge){
            solid_body_collection.Add_Force(Create_Edge_Springs(*triangulated_surface,linear_stiffness,linear_damping));} // were *2 and *10
        if(use_bending){
            T axial_bending_stiffness=axial_bending_stiffness_multiplier*2/(1+sqrt((T)2)),axial_bending_damping=axial_bending_damping_multiplier*8;
            solid_body_collection.Add_Force(Create_Axial_Bending_Springs(*triangulated_surface,(T).01,axial_bending_stiffness,axial_bending_damping));
            T bending_stiffness=bending_stiffness_multiplier*2/(1+sqrt((T)2)),bending_damping=bending_damping_multiplier*8;
            solid_body_collection.Add_Force(Create_Bending_Springs(*triangulated_surface,bending_stiffness,bending_damping));}
        if(use_altitude) solid_body_collection.Add_Force(Create_Altitude_Springs(*triangulated_surface,altitude_stiffness,altitude_damping));}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE 
{
    if(test_number==6){
        for(int i=3;i<9;i++) for(int j=3;j<9;j++){
            int linear_index=i*12+j;
            V(linear_index)=TV();}}
    else if(test_number==8){
        V(0)=TV();
        V((number_side_panels+1)*number_side_panels+1)=TV();}
}
//#####################################################################
};
}
#endif

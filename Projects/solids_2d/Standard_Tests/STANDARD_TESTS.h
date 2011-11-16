//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
//   1. Quasistatic mattress stretched
//   2. Quasistatic mattress stretched, refinement
//   3. Quasistatic mattress stretched, refinement, with T junctions
//   4. Dynamic mattress
//   5. Dynamic mattress stretched, refinement
//   6. Dynamic mattress stretched, refinement, with T junctions
//   7. Falling circle with bending elements.
//   8. Falling mattress.
//   9. Incompressible falling mattress.
//  10. Falling particle on segment ramp.
//  11. String of segments falling on to ground.
//  12. Curtain and ball
//  13. Falling mattress, random start
//  14. Several falling mattresses
//  15. Penalty collision test
//  16. Crush test
//  17. Matress, no gravity, random start
//  18. Matress, no gravity, point start
//  19. Matress, no gravity, line start
//  20. Skinny mattress falling on box
//  21. Stretched skinny mattress with colliding ball
//  22. Squeezing between to small boxes
//  23. Highspeed collision with ground
//  24. Big 4 sides stretch
//  25. Big 4 corners stretch
//  26. Big stretch/bend
//  27. Force inversion
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <fstream>

#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parallel_Computation/PARTITION_ID.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/COROTATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_COROTATED_BLEND.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED_REFINED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ROTATED_LINEAR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/COLLISION_AREA_PENALTY_FORCE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/INCOMPRESSIBLE_FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_BENDING_ELEMENTS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_STATE.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/QUASISTATIC_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,2> > >
{
    typedef T_input T;
    typedef VECTOR<T,2> TV;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::output_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::solid_body_collection;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions;using BASE::solids_evolution; // silence -Woverloaded-virtual
    using BASE::parse_args;using BASE::test_number;

    std::ofstream svout;

    SOLIDS_STANDARD_TESTS<TV> tests;

    ARRAY<TV> deformable_body_rest_positions;

    // test 1,2,3,4,5,6
    GRID<TV> mattress_grid;
    TV attachment_velocity;

    VECTOR<int,2> processes_per_dimension;

    ARRAY<int> segment_ramp_particles;
    bool fully_implicit;
    bool test_forces;
    bool use_extended_neohookean;
    bool use_extended_neohookean_refined;
    bool use_extended_neohookean_hyperbola;
    bool use_corotated;
    bool use_corot_blend;
    bool dump_sv;
    int kinematic_id;
    INTERPOLATION_CURVE<T,FRAME<TV> > curve;
    bool print_matrix;
    int parameter;
    T stiffness_multiplier;
    T damping_multiplier;

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),fully_implicit(false),test_forces(false),use_extended_neohookean(false),use_extended_neohookean_refined(false),use_extended_neohookean_hyperbola(false),use_corotated(false),use_corot_blend(false),dump_sv(false),
        print_matrix(false),parameter(0),stiffness_multiplier(1),damping_multiplier(1)
    {
    }

    // Unused callbacks
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {if(dump_sv)svout.close();}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Update_Fragments() PHYSBAM_OVERRIDE {}
    void Post_Initialization() PHYSBAM_OVERRIDE {}
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
    {
        if(test_forces){
            solid_body_collection.deformable_body_collection.Test_Energy(time);
            solid_body_collection.deformable_body_collection.Test_Force_Derivatives(time);}
    }
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
    {
        if (dump_sv)
        {
            FINITE_VOLUME<TV,2>& force_field = solid_body_collection.deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,2>&>();
            ARRAY<DIAGONAL_MATRIX<T,2> >& sv = force_field.Fe_hat;
            
            for (int i=1; i<=sv.m; i++)
            {
                svout << sv(i).x11 << " " << sv(i).x22 << std::endl;
            }
        }
    }
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Option_Argument("-fully_implicit","use fully implicit forces");
    parse_args->Add_Option_Argument("-test_forces","use fully implicit forces");
    parse_args->Add_Option_Argument("-use_ext_neo");
    parse_args->Add_Option_Argument("-use_ext_neo_ref");
    parse_args->Add_Option_Argument("-use_ext_neo_hyper");
    parse_args->Add_Option_Argument("-use_corotated");
    parse_args->Add_Option_Argument("-use_corot_blend");
    parse_args->Add_Option_Argument("-dump_sv");
    parse_args->Add_Integer_Argument("-parameter",0,"parameter used by multiple tests to change the parameters of the test");
    parse_args->Add_Double_Argument("-stiffen",1,"","stiffness multiplier for various tests");
    parse_args->Add_Double_Argument("-dampen",1,"","damping multiplier for various tests");
    parse_args->Add_Option_Argument("-residuals","print residuals during timestepping");
    parse_args->Add_Option_Argument("-print_energy","print energy statistics");
    parse_args->Add_Double_Argument("-cgsolids",1e-2,"CG tolerance for backward Euler");
    parse_args->Add_Option_Argument("-use_be","use backward euler");
    parse_args->Add_Option_Argument("-print_matrix");
    parse_args->Add_Option_Argument("-project_nullspace","project out nullspace");
    parse_args->Add_Integer_Argument("-projection_iterations",5,"number of iterations used for projection in cg");
    parse_args->Add_Integer_Argument("-solver_iterations",1000,"number of iterations used for solids system");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
    output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d",test_number);
    last_frame=1000;

    switch(test_number){
	case 20: case 21: case 26:
	    mattress_grid=GRID<TV>(40,8,(T)-2,(T)2,(T)-.4,(T).4);
	break;
	case 22: case 23: case 24: case 25: case 27:
	    mattress_grid=GRID<TV>(20,20,(T)-.9,(T).9,(T)-.9,(T).9);
	break;
    	default:
            mattress_grid=GRID<TV>(20,10,(T)-1,(T)1,(T)-.5,(T).5);
    }
    
    processes_per_dimension=VECTOR<int,2>(2,1);

    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    fully_implicit=parse_args->Is_Value_Set("-fully_implicit");
    test_forces=parse_args->Is_Value_Set("-test_forces");
    use_extended_neohookean=parse_args->Is_Value_Set("-use_ext_neo");
    use_extended_neohookean_refined=parse_args->Is_Value_Set("-use_ext_neo_ref"); // 
    use_corotated=parse_args->Is_Value_Set("-use_corotated");
    use_extended_neohookean_hyperbola=parse_args->Is_Value_Set("-use_ext_neo_hyper");
    use_corot_blend=parse_args->Is_Value_Set("-use_corot_blend");
    dump_sv=parse_args->Is_Value_Set("-dump_sv");
    solids_parameters.use_trapezoidal_rule_for_velocities=!parse_args->Get_Option_Value("-use_be");
    print_matrix=parse_args->Is_Value_Set("-print_matrix");
    parameter=parse_args->Get_Integer_Value("-parameter");
    stiffness_multiplier=(T)parse_args->Get_Double_Value("-stiffen");
    damping_multiplier=(T)parse_args->Get_Double_Value("-dampen");
    solid_body_collection.print_energy=parse_args->Get_Option_Value("-print_energy");
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=parse_args->Get_Integer_Value("-projection_iterations");
    if(parse_args->Is_Value_Set("-project_nullspace")) solids_parameters.implicit_solve_parameters.project_nullspace_frequency=1;
    solid_body_collection.Print_Residuals(parse_args->Get_Option_Value("-residuals"));

    switch(test_number){
        case 1:
        case 2:
        case 3:
            delete solids_evolution;
            solids_evolution=new QUASISTATIC_EVOLUTION<TV>(solids_parameters,solid_body_collection);
            solids_parameters.newton_iterations=10;
            solids_parameters.newton_tolerance=(T)1e-2;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=900;
            solids_parameters.use_partially_converged_result=true;
            solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
            attachment_velocity=TV((T).1,0);
            last_frame=120;
            break;
        case 4:
        case 5:
        case 6:
        case 24:
        case 25:
        case 26:
        case 27:
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=900;
            solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
            attachment_velocity=TV((T).1,0);
            last_frame=1500;
            break;
        case 7:
            solids_parameters.cfl=(T).5;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            last_frame=120;
            break;
        case 8: 
        case 13:
        case 14:
        case 15:
        case 16:
        case 22:
	case 21:
        case 17:
        case 18:
        case 19:
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            last_frame=200;
            break;
        case 23: 
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            last_frame=500;
            break;
        case 9:
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            solids_parameters.cfl=1;
            solids_parameters.deformable_object_collision_parameters.collide_with_interior=true;
            last_frame=120;
            break;
        case 10:
        case 11:
        case 12:
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            solids_parameters.cfl=(T)4.9;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            last_frame=500;
            frame_rate=120;
            break;
	case 20:
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=900;
            solids_parameters.deformable_object_collision_parameters.collide_with_interior=true;
            //solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
            attachment_velocity=TV((T).8,0);
	    last_frame=400;
            break;	    
        default:
            LOG::cerr<<"Unrecognized test number "<<test_number<<std::endl;exit(1);}

    output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d",test_number);
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    switch(test_number){
        case 1:
            tests.Create_Mattress(mattress_grid,true);
            break;
        case 4:
        case 24:
        case 25:
        case 26:
        case 27:
            tests.Create_Mattress(mattress_grid,true);
            break;
        case 2:
        case 3:
        case 5:
        case 6:{
            TRIANGULATED_AREA<T>& triangulated_area=*TRIANGULATED_AREA<T>::Create();
            TRIANGLE_MESH& mesh=triangulated_area.mesh;
            VECTOR<int,2> refine(3,8);int number_new_particles=(refine.y-refine.x+1)*(mattress_grid.counts.y-1);
            ARRAY<int> deleted_elements;deleted_elements.Preallocate(2*number_new_particles);
            ARRAY<VECTOR<int,3> > new_elements;new_elements.Preallocate(4*number_new_particles);
            triangulated_area.particles.array_collection->Preallocate(number_new_particles);triangulated_area.Initialize_Square_Mesh_And_Particles(mattress_grid);
            for(int t=1;t<=mesh.elements.m;t+=2){
                VECTOR<int,3> nodes1=mesh.elements(t);int i=nodes1[1]%mattress_grid.counts.x,j_minus_one=nodes1[1]/mattress_grid.counts.x;
                if(i>refine.x+j_minus_one&&i<refine.y+j_minus_one){
                    VECTOR<int,3> nodes2=mesh.elements(t+1);
                    int center_node=triangulated_area.particles.array_collection->Add_Element();
                    triangulated_area.particles.X(center_node)=mattress_grid.Center(i,j_minus_one+1);
                    deleted_elements.Append(t);deleted_elements.Append(t+1);
                    new_elements.Append(VECTOR<int,3>(nodes1[1],nodes1[2],center_node));
                    new_elements.Append(VECTOR<int,3>(nodes1[2],nodes1[3],center_node));
                    new_elements.Append(VECTOR<int,3>(nodes2[2],nodes2[3],center_node));
                    new_elements.Append(VECTOR<int,3>(nodes2[3],nodes2[1],center_node));}}
            if(test_number==3||test_number==6){ // add the bound non-dynamical particles
                for(int t=1;t<=mesh.elements.m;t+=2){
                    VECTOR<int,3> nodes1=mesh.elements(t);int i=nodes1[1]%mattress_grid.counts.x,j_minus_one=nodes1[1]/mattress_grid.counts.x;
                    if(i==refine.x+j_minus_one||i==refine.y+j_minus_one){
                        int center_node=triangulated_area.particles.array_collection->Add_Element();
                        triangulated_area.particles.X(center_node)=mattress_grid.Center(i,j_minus_one+1);
                        solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new LINEAR_BINDING<TV,2>(particles,center_node,VECTOR<int,2>(nodes1[1],nodes1[3]),TV((T).5,(T).5)));
                        if(i==refine.x+j_minus_one){
                            deleted_elements.Append(t);
                            new_elements.Append(VECTOR<int,3>(nodes1[1],nodes1[2],center_node));
                            new_elements.Append(VECTOR<int,3>(nodes1[2],nodes1[3],center_node));}
                        else{
                            VECTOR<int,3> nodes2=mesh.elements(t+1);
                            deleted_elements.Append(t+1);
                            new_elements.Append(VECTOR<int,3>(nodes2[2],nodes2[3],center_node));
                            new_elements.Append(VECTOR<int,3>(nodes2[3],nodes2[1],center_node));}}}
                mesh.Delete_Elements(deleted_elements);mesh.elements.Append_Elements(new_elements);}
            else{mesh.Delete_Sorted_Elements(deleted_elements);mesh.elements.Append_Elements(new_elements);}
            triangulated_area.Check_Signed_Areas_And_Make_Consistent(true);
            LOG::cout<<"Adding Deformable body - Total Triangles : "<<triangulated_area.mesh.elements.m<<std::endl;
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(triangulated_area,100);
            tests.Copy_And_Add_Structure(triangulated_area);
            break;}
        case 7:{
            tests.Create_Segmented_Curve(300,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,4))));
            tests.Add_Ground();
            break;}
        case 8: 
        case 13:{
            tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,4))));
            tests.Add_Ground();
            break;}
        case 23:{
            tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,200))));
            tests.Add_Ground();
            break;}
        case 9:{
            tests.Create_Mattress(mattress_grid,false,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T).5))));
            tests.Add_Ground();
            break;}
        case 10:{
            SEGMENTED_CURVE<TV>& segmented_curve=*SEGMENTED_CURVE<TV>::Create();segmented_curve.Clean_Memory();
            int num_particles=51;T scale_factor=(T).1;segmented_curve.mesh.Initialize_Straight_Mesh(num_particles,false);
            for(int p=1;p<=num_particles;p++){
                int index=segmented_curve.particles.array_collection->Add_Element();segment_ramp_particles.Append(index);
                T x=(p-num_particles/2-1)*scale_factor;segmented_curve.particles.X(index)=VECTOR<T,2>((T)p,x*x*x);}
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(segmented_curve,100);
            tests.Set_Initial_Particle_Configuration(segmented_curve.particles,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,20))),true);
            tests.Copy_And_Add_Structure(segmented_curve);
            FREE_PARTICLES<TV>& free_particles=*FREE_PARTICLES<TV>::Create();solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&free_particles);
            int new_free_node=particles.array_collection->Add_Element();particles.X(new_free_node)=TV(20,35);free_particles.nodes.Append(new_free_node);
            particles.mass(new_free_node)=particles.mass(1);
            tests.Add_Ground();
            break;}
        case 11:{
            SEGMENTED_CURVE<TV>& segmented_curve=*SEGMENTED_CURVE<TV>::Create(*new PARTICLES<TV>);segmented_curve.Clean_Memory();
            int num_particles=500;segmented_curve.mesh.Initialize_Straight_Mesh(num_particles,false);segmented_curve.particles.array_collection->Add_Elements(num_particles);
            for(int p=1;p<=num_particles;p++) segmented_curve.particles.X(p)=VECTOR<T,2>(0,(T)p/num_particles+(T).01);
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(segmented_curve,1);
            tests.Set_Initial_Particle_Configuration(segmented_curve.particles,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T).5),ROTATION<TV>::From_Angle((T).01))),true);
            tests.Copy_And_Add_Structure(segmented_curve);
            tests.Add_Ground();
            break;}
        case 14:{
            tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,4))));
            tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(4,4))));
            tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(8,4))));
            tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(12,4))));
            tests.Add_Ground();
            break;}
        case 15:{
            for(int i=1;i<=parameter;i++) tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,1.2*i))));
            tests.Add_Ground();
            break;}
        case 16: {
            tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,1))));
            RIGID_BODY<TV>& box1=tests.Add_Rigid_Body("square",10,(T)0);
            RIGID_BODY<TV>& box2=tests.Add_Rigid_Body("square",10,(T)0);
            box1.X()=TV(0,-10);
            box2.X()=TV(0,12);
            box1.is_static=true;
            box2.is_static=false;
            kinematic_id=box2.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(box2.particle_index)=true;
            curve.Add_Control_Point(0,FRAME<TV>(TV(0,12)));
            curve.Add_Control_Point(5,FRAME<TV>(TV(0,10)));
            curve.Add_Control_Point(6,FRAME<TV>(TV(0,10)));
            curve.Add_Control_Point(11,FRAME<TV>(TV(0,12)));
            last_frame=250;
            break;}
        case 22: {
            tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,0))));
            RIGID_BODY<TV>& box1=tests.Add_Rigid_Body("square",1,(T)0);
            RIGID_BODY<TV>& box2=tests.Add_Rigid_Body("square",1,(T)0);
            box1.X()=TV(0,-2);
            box2.X()=TV(0,2);
            box1.is_static=true;
            box2.is_static=false;
            kinematic_id=box2.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(box2.particle_index)=true;
            curve.Add_Control_Point(0,FRAME<TV>(TV(0,2)));
            curve.Add_Control_Point(5,FRAME<TV>(TV(0,0)));
            curve.Add_Control_Point(6,FRAME<TV>(TV(0,0)));
            curve.Add_Control_Point(11,FRAME<TV>(TV(0,2)));
            last_frame=250;
            break;}
        case 17:
        case 18:
        case 19:{
            tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,0))));
            break;}
	case 20:{
            tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,4))));
            //RIGID_BODY<TV>& box1=tests.Add_Rigid_Body("circle",4,(T)0);
            RIGID_BODY<TV>& box2=tests.Add_Rigid_Body("square",.2,(T)0);
            //box1.X()=TV(0,-10);
            box2.X()=TV(0,12);
            //box1.is_static=true;
            box2.is_static=false;
            kinematic_id=box2.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(box2.particle_index)=true;
            curve.Add_Control_Point(0,FRAME<TV>(TV(0,12)));
            curve.Add_Control_Point(5,FRAME<TV>(TV(0,0)));
            curve.Add_Control_Point(6,FRAME<TV>(TV(0,0)));
            curve.Add_Control_Point(11,FRAME<TV>(TV(0,12)));
	    break;}
	case 21:{
            tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,4))));
            RIGID_BODY<TV>& box1=tests.Add_Rigid_Body("square",1,(T)0);
            box1.X()=TV(0,-6);
            box1.is_static=true;
	    break;
}
        default:
            LOG::cerr<<"Missing implementation for test number "<<test_number<<std::endl;exit(1);}

    // add structures and rigid bodies to collisions
    solid_body_collection.deformable_body_collection.collisions.collision_structures.Append_Elements(solid_body_collection.deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(solid_body_collection.deformable_body_collection.deformable_geometry.structures);

    // correct number nodes
    for(int i=1;i<=solid_body_collection.deformable_body_collection.deformable_geometry.structures.m;i++) solid_body_collection.deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    switch(test_number){
        case 1:
        case 2:
        case 3:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true));
            FINITE_VOLUME<TV,2>* fvm=Create_Quasistatic_Finite_Volume(triangulated_area,new NEO_HOOKEAN<T,2>((T)1e5,(T).45));
            solid_body_collection.Add_Force(fvm);
            deformable_body_rest_positions=particles.X;
            break;}
        case 4:
        case 24:
        case 25:
        case 26:
        case 27:
        case 5:
        case 6:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(triangulated_area,(T)1e2,(T).45,(T).05);
            break;}
        case 7:
        case 10:
        case 11:{
            SEGMENTED_CURVE<TV>& segmented_curve=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<SEGMENTED_CURVE<TV>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true));
            solid_body_collection.Add_Force(Create_Edge_Springs(segmented_curve,(T)5e1,(T)1.5));
            solid_body_collection.Add_Force(Create_Bending_Elements(segmented_curve,(T)1e-6,(T).0001));
            break;}
        case 8: 
        case 16:
        case 23:
        case 13:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true));
            Add_Constitutive_Model(triangulated_area,(T)1e4,(T).45,(T).01);
            if(test_number==13){RANDOM_NUMBERS<T> rand;rand.Fill_Uniform(particles.X,-1,1);}
            break;}
        case 22:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(triangulated_area,(T)1e4,(T).45,(T).01);
            break;}
        case 14:{
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true));
            
            TRIANGULATED_AREA<T>& triangulated_area1=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>(1);
            solid_body_collection.Add_Force(Create_Finite_Volume(triangulated_area1,new COROTATED<T,2>((T)2e4,(T).45,(T).01)));

            TRIANGULATED_AREA<T>& triangulated_area2 = solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>(2);
            solid_body_collection.Add_Force(Create_Finite_Volume(triangulated_area2,new NEO_HOOKEAN_COROTATED_BLEND<T,2>((T)2e4,(T).45,(T).01)));

            TRIANGULATED_AREA<T>& triangulated_area3 = solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>(3);
            solid_body_collection.Add_Force(Create_Finite_Volume(triangulated_area3,new NEO_HOOKEAN_EXTRAPOLATED<T,2>((T)2e4,(T).45,(T).01)));

            TRIANGULATED_AREA<T>& triangulated_area4 = solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>(4);
            solid_body_collection.Add_Force(Create_Finite_Volume(triangulated_area4,new NEO_HOOKEAN<T,2>((T)2e4,(T).45,(T).01)));
            break;}
        case 9:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true));
            solid_body_collection.Add_Force(Create_Incompressible_Finite_Volume(triangulated_area));
            Add_Constitutive_Model(triangulated_area,(T)2e3,(T)0,(T).01);
            break;}
        case 15:{
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true));
            COLLISION_AREA_PENALTY_FORCE<TV>* penalty_force=new COLLISION_AREA_PENALTY_FORCE<TV>(particles);
            for(int i=1;i<=parameter;i++){
                TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>(i);
                penalty_force->Add_Mesh(triangulated_area);
                Add_Constitutive_Model(triangulated_area,(T)1e4,(T).45,(T).01);}
            solid_body_collection.Add_Force(penalty_force);
            break;}
        case 17:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(triangulated_area,(T)1e4,(T).45,(T).01);

            RANDOM_NUMBERS<T> rand;
            rand.Fill_Uniform(particles.X,-1,1);
            break;}
        case 18:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(triangulated_area,(T)1e4,(T).45,(T).01);
            RANDOM_NUMBERS<T> rand;
            rand.Fill_Uniform(particles.X,0,0);
            break;}
        case 19:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(triangulated_area,(T)1e4,(T).45,(T).01);
            break;}
        case 20:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>(1);
            Add_Constitutive_Model(triangulated_area,(T)1e5,(T).45,(T).01);
            break;}
        case 21:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>(1);
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true));
            Add_Constitutive_Model(triangulated_area,(T)1e4,(T).45,(T).01);
            break;}

        default:
            LOG::cerr<<"Missing implementation for test number "<<test_number<<std::endl;exit(1);}

    if(solid_body_collection.deformable_body_collection.mpi_solids){
        solid_body_collection.deformable_body_collection.mpi_solids->Simple_Partition(solid_body_collection.deformable_body_collection,solid_body_collection.rigid_body_collection.rigid_geometry_collection,particles.X,processes_per_dimension);
        switch(test_number){
            case 1:
            case 24: 
            case 25: 
            case 26: 
            case 27: 
            case 4: 
                break;
            case 3:
                for(PARTITION_ID r(1);r<=solid_body_collection.deformable_body_collection.mpi_solids->particles_of_partition.Size();r++){
                    ARRAY<int> deletion_list;
                    for(int i=1;i<=solid_body_collection.deformable_body_collection.mpi_solids->particles_of_partition(r).m;i++){int p=solid_body_collection.deformable_body_collection.mpi_solids->particles_of_partition(r)(i);
                        if(solid_body_collection.deformable_body_collection.binding_list.Binding(p)) deletion_list.Append(i);}
                    solid_body_collection.deformable_body_collection.mpi_solids->particles_of_partition(r).Remove_Sorted_Indices(deletion_list);}
                break;
            default:
                LOG::cerr<<"Missing implementation for test number "<<test_number<<std::endl;exit(1);}}
    solid_body_collection.Update_Simulated_Particles();

    if(parse_args->Is_Value_Set("-solver_iterations")) solids_parameters.implicit_solve_parameters.cg_iterations=parse_args->Get_Integer_Value("-solver_iterations");
    if(parse_args->Is_Value_Set("-cgsolids")) solids_parameters.implicit_solve_parameters.cg_tolerance=(T)parse_args->Get_Double_Value("-cgsolids");

    if(fully_implicit) for(int i=1;i<=solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->use_implicit_velocity_independent_forces=true;
    if(fully_implicit) for(int i=1;i<=solid_body_collection.rigid_body_collection.rigids_forces.m;i++)
        solid_body_collection.rigid_body_collection.rigids_forces(i)->use_implicit_velocity_independent_forces=true;
    if(fully_implicit) for(int i=1;i<=solid_body_collection.deformable_body_collection.deformables_forces.m;i++) solid_body_collection.deformable_body_collection.deformables_forces(i)->use_implicit_velocity_independent_forces=true;

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Initialize_Bodies();
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id)
{
    if(id==kinematic_id) frame=curve.Value(time);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id)
{
    if(id==kinematic_id) twist=curve.Derivative(time);
    return false;
}
//#####################################################################
// Function Set_Particle_Is_Simulated
//#####################################################################
void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE
{
    if(test_number==10){
        // make free particles simulated
        DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
        FREE_PARTICLES<TV>& free_particles=deformable_body_collection.deformable_geometry.template Find_Structure<FREE_PARTICLES<TV>&>();
        for(int i=1;i<=free_particles.nodes.Size();i++) particle_is_simulated(free_particles.nodes(i))=true;}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(test_number==4||test_number==5||test_number==6||test_number==20){
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
        TV velocity=velocity_time<5.0?attachment_velocity:TV();
        for(int j=1;j<=n;j++){V(1+m*(j-1))=-velocity;V(m+m*(j-1))=velocity;}}
    if(test_number==24){
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
        T velocity=.2;
        T final_time=50;
        TV velocity_top=velocity_time<final_time?TV((T)0,velocity):TV();
        TV velocity_bot=velocity_time<final_time?TV((T)0,-velocity):TV();
        TV velocity_rig=velocity_time<final_time?TV(velocity,(T)0):TV();
        TV velocity_lef=velocity_time<final_time?TV(-velocity,(T)0):TV();
        for(int j=n/3+1;j<=2*n/3+1;j++){V(1+m*(j-1))=velocity_lef;V(m+m*(j-1))=velocity_rig;}
        for(int i=m/3+1;i<=2*m/3+1;i++){V(i)=velocity_bot;V(m*(n-1)+i)=velocity_top;}}
    if(test_number==25){
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
        T velocity=.2;
        T final_time=50;
        TV velocity_tr=velocity_time<final_time?TV(velocity,velocity):TV();
        TV velocity_tl=velocity_time<final_time?TV(-velocity,velocity):TV();
        TV velocity_br=velocity_time<final_time?TV(velocity,-velocity):TV();
        TV velocity_bl=velocity_time<final_time?TV(-velocity,-velocity):TV();
        V(1)=velocity_bl;
        V(m)=velocity_br;
        V(m*(n-1)+1)=velocity_tl;
        V(m*n)=velocity_tr;}
    if(test_number==26){
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
        T velocity=.2;
        T final_time=50;
        TV velocity_top=velocity_time<final_time?TV((T)0,-velocity):TV();
        TV velocity_bot=velocity_time<final_time?TV((T)0,-velocity):TV();
        TV velocity_rig=velocity_time<final_time?TV(velocity,(T)0):TV();
        TV velocity_lef=velocity_time<final_time?TV(-velocity,(T)0):TV();
        for(int j=1;j<=n;j++){V(1+m*(j-1))=velocity_lef;V(m+m*(j-1))=velocity_rig;}
        for(int i=2*m/5+1;i<=3*m/5+1;i++){V(i)=velocity_bot;V(m*(n-1)+i)=velocity_top;}}
    if(test_number==27){
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
        T velocity=.1;
        T final_time=50;
        TV velocity_rig=velocity_time<final_time?TV(-velocity,(T)0):TV();
        TV velocity_lef=velocity_time<final_time?TV(velocity,(T)0):TV();
        for(int j=1;j<=n;j++){V(1+m*(j-1))=velocity_lef;V(m+m*(j-1))=velocity_rig;}}
    if(test_number==10) for(int p=1;p<=segment_ramp_particles.m;p++) V(p)=TV();
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(test_number==4||test_number==5||test_number==6||test_number==20){
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
        for(int j=1;j<=n;j++) V(1+m*(j-1))=V(m+m*(j-1))=TV();}
    if(test_number==24){
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
        for(int j=n/3+1;j<=2*n/3+1;j++){V(1+m*(j-1))=V(m+m*(j-1))=TV();}
        for(int i=m/3+1;i<=2*m/3+1;i++){V(i)=V(m*(n-1)+i)=TV();}}
    if(test_number==25){
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
        V(1)=TV();
        V(m)=TV();
        V(m*(n-1)+1)=TV();
        V(m*n)=TV();}
    if(test_number==26){
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
        for(int j=1;j<=n;j++){V(1+m*(j-1))=TV();V(m+m*(j-1))=TV();}
        for(int i=2*m/5+1;i<=3*m/5+1;i++){V(i)=TV();V(m*(n-1)+i)=TV();}}
    if(test_number==27){
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
        for(int j=1;j<=n;j++){V(1+m*(j-1))=TV();V(m+m*(j-1))=TV();}}
    if(test_number==10) for(int p=1;p<=segment_ramp_particles.m;p++) V(p)=TV();
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE
{
    if(test_number!=1 && test_number!=2 && test_number!=3) return;
    ARRAY<TV>& X_save=deformable_body_rest_positions;
    int m=mattress_grid.counts.x;
    for(int j=1;j<=mattress_grid.counts.y;j++){
        X(1+m*(j-1))=X_save(1+m*(j-1))-time*attachment_velocity;
        X(m+m*(j-1))=X_save(m+m*(j-1))+time*attachment_velocity;}
}
//#####################################################################
// Function Zero_Out_Enslaved_Position_Nodes
//#####################################################################
void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE
{
    int m=mattress_grid.counts.x;
    for(int j=1;j<=mattress_grid.counts.y;j++) X(1+m*(j-1))=X(m+m*(j-1))=TV(0,0);
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const
{
    BASE::Write_Output_Files(frame);
    tests.Write_Debug_Particles(output_directory,frame);
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame)
{
    dynamic_cast<NEWMARK_EVOLUTION<TV>&>(*solids_evolution).print_matrix=print_matrix;
    
    if (dump_sv)
    {
        std::string output_file = STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d/SV_%d",test_number,frame);
        svout.open(output_file.c_str());
    }
}
//#####################################################################
// Function Add_Constitutive_Model
//#####################################################################
void Add_Constitutive_Model(TRIANGULATED_AREA<T>& triangulated_area,T stiffness,T poissons_ratio,T damping)
{
    ISOTROPIC_CONSTITUTIVE_MODEL<T,2>* icm=0;
    if(use_extended_neohookean) icm=new NEO_HOOKEAN_EXTRAPOLATED<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,.4,20*stiffness*stiffness_multiplier);
    else if(use_extended_neohookean_refined) icm=new NEO_HOOKEAN_EXTRAPOLATED_REFINED<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,.4,.6,20*stiffness*stiffness_multiplier);
    else if(use_corotated) icm=new COROTATED<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else if(use_corot_blend) icm=new NEO_HOOKEAN_COROTATED_BLEND<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else if(use_extended_neohookean_hyperbola) icm=new NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else icm=new NEO_HOOKEAN<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    solid_body_collection.Add_Force(Create_Finite_Volume(triangulated_area,icm));
}
//#####################################################################
};
}
#endif

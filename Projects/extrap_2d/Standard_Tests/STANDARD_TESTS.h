//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
//   8. Falling mattress.
//  13. Falling mattress, random start
//  14. Several falling mattresses
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
//  270. Inverted configuration
//  28. Taffy test
//  29. Gear Test
//  30. Horizontal stretch
//  31. Triangle stretch
//  32. Triangle stretch (II)
//  33. Tangle plot
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <PhysBAM_Tools/Images/PPM_FILE.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Basic_Geometry/SMOOTH_GEAR.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/COROTATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/COROTATED_QUARTIC.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/GEN_NEO_HOOKEAN_ENERGY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/GENERAL_EXTRAPOLATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_COROTATED_BLEND.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED_REFINED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED_SMOOTH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED2.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_J_INTERP_ENERGY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/RC_EXTRAPOLATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/RC2_EXTRAPOLATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ROTATED_LINEAR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Read_Write/EPS_FILE_GEOMETRY.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <fstream>
namespace PhysBAM{
template<class TV> void Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color);
template<class TV,class ATTR> void Debug_Particle_Set_Attribute(ATTRIBUTE_ID id,const ATTR& attr);

template<class T_input>
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,2> > >
{
    typedef T_input T;
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::output_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::solid_body_collection;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions;using BASE::solids_evolution;
    using BASE::parse_args;using BASE::test_number;

    std::ofstream svout;

    SOLIDS_STANDARD_TESTS<TV> tests;

    GRID<TV> mattress_grid;
    TV attachment_velocity;
    bool semi_implicit;
    bool test_forces;
    bool use_extended_neohookean;
    bool use_extended_neohookean2;
    bool use_extended_neohookean3;
    bool use_int_j_neo;
    bool use_rc_ext;
    bool use_rc2_ext;
    bool use_extended_neohookean_refined;
    bool use_extended_neohookean_hyperbola;
    bool use_extended_neohookean_smooth;
    bool use_corotated;
    bool use_corot_blend;
    bool use_corot_quartic;
    bool dump_sv;
    int kinematic_id,kinematic_id2;
    INTERPOLATION_CURVE<T,FRAME<TV> > curve,curve2;
    bool print_matrix;
    int parameter;
    T stiffness_multiplier;
    T damping_multiplier;
    bool use_constant_ife;
    T stretch;
    T primary_contour;
    T sigma_range;
    int image_size;
    T poissons_ratio;
    bool scatter_plot;
    bool use_contrails;
    ARRAY<ARRAY<TV> > contrail;
    ARRAY<VECTOR<T,3> > contrail_colors;
    ARRAY<VECTOR<TV,2> > contour_segments;
    T input_cutoff;
    T input_efc;
    T input_poissons_ratio,input_youngs_modulus;
    bool test_model_only;
    bool plot_energy_density;

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),semi_implicit(false),test_forces(false),use_extended_neohookean(false),
        use_extended_neohookean_refined(false),use_extended_neohookean_hyperbola(false),use_extended_neohookean_smooth(false),use_corotated(false),
        use_corot_blend(false),use_corot_quartic(false),dump_sv(false),
        print_matrix(false),parameter(20),stiffness_multiplier(1),damping_multiplier(1),use_constant_ife(false),stretch(1),poissons_ratio((T).45),input_cutoff(0),input_efc(0),
        input_poissons_ratio(-1),input_youngs_modulus(0),test_model_only(false)
    {
    }

    // Unused callbacks
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE
    {
        if(dump_sv) svout.close();
        if(scatter_plot) Dump_Scatter_Plot(frame);
    }
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
        if(scatter_plot) Update_Scatter_Plot();
    }
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    //void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Option_Argument("-semi_implicit","use semi implicit forces");
    parse_args->Add_Option_Argument("-test_forces","use fully implicit forces");
    parse_args->Add_Option_Argument("-use_ext_neo");
    parse_args->Add_Option_Argument("-use_ext_neo2");
    parse_args->Add_Option_Argument("-use_ext_neo3");
    parse_args->Add_Option_Argument("-use_int_j_neo");
    parse_args->Add_Option_Argument("-use_rc_ext");
    parse_args->Add_Option_Argument("-use_rc2_ext");
    parse_args->Add_Option_Argument("-use_ext_neo_ref");
    parse_args->Add_Option_Argument("-use_ext_neo_hyper");
    parse_args->Add_Option_Argument("-use_ext_neo_smooth");
    parse_args->Add_Option_Argument("-use_corotated");
    parse_args->Add_Option_Argument("-use_corot_blend");
    parse_args->Add_Option_Argument("-use_corot_quartic");
    parse_args->Add_Option_Argument("-dump_sv");
    parse_args->Add_Integer_Argument("-parameter",0,"parameter used by multiple tests to change the parameters of the test");
    parse_args->Add_Double_Argument("-stiffen",1,"","stiffness multiplier for various tests");
    parse_args->Add_Double_Argument("-dampen",1,"","damping multiplier for various tests");
    parse_args->Add_Option_Argument("-residuals","print residuals during timestepping");
    parse_args->Add_Option_Argument("-print_energy","print energy statistics");
    parse_args->Add_Double_Argument("-cgsolids",1e-3,"CG tolerance for backward Euler");
    parse_args->Add_Option_Argument("-use_be","use backward euler");
    parse_args->Add_Option_Argument("-print_matrix");
    parse_args->Add_Option_Argument("-project_nullspace","project out nullspace");
    parse_args->Add_Integer_Argument("-projection_iterations",5,"number of iterations used for projection in cg");
    parse_args->Add_Integer_Argument("-solver_iterations",1000,"number of iterations used for solids system");
    parse_args->Add_Option_Argument("-use_constant_ife","use constant extrapolation on inverting finite element fix");
    parse_args->Add_Double_Argument("-stretch",1,"stretch");
    parse_args->Add_Option_Argument("-plot_contour","plot primary contour");
    parse_args->Add_Integer_Argument("-image_size",500,"image size for plots");
    parse_args->Add_Double_Argument("-sigma_range",5,"sigma range for plots");
    parse_args->Add_Double_Argument("-poissons_ratio",.45,"poisson's ratio");
    parse_args->Add_Option_Argument("-scatter_plot","Create contrail plot with singular values");
    parse_args->Add_Option_Argument("-use_contrails","Show contrails in plot");
    parse_args->Add_Double_Argument("-cutoff",.4,"cutoff");
    parse_args->Add_Double_Argument("-efc",20,"efc");
    parse_args->Add_Double_Argument("-poissons_ratio",-1,"poissons_ratio");
    parse_args->Add_Double_Argument("-youngs_modulus",0,"youngs modulus, only for test 41 so far");
    parse_args->Add_Option_Argument("-test_model_only");
    parse_args->Add_Option_Argument("-plot_energy_density");
    parse_args->Add_Option_Argument("-test_system");
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
    parameter=parse_args->Get_Integer_Value("-parameter");
    
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    semi_implicit=parse_args->Is_Value_Set("-semi_implicit");
    test_forces=parse_args->Is_Value_Set("-test_forces");
    use_extended_neohookean=parse_args->Is_Value_Set("-use_ext_neo");
    use_extended_neohookean2=parse_args->Is_Value_Set("-use_ext_neo2");
    use_extended_neohookean3=parse_args->Is_Value_Set("-use_ext_neo3");
    use_int_j_neo=parse_args->Is_Value_Set("-use_int_j_neo");
    use_rc_ext=parse_args->Is_Value_Set("-use_rc_ext");
    use_rc2_ext=parse_args->Is_Value_Set("-use_rc2_ext");
    use_extended_neohookean_refined=parse_args->Is_Value_Set("-use_ext_neo_ref"); //
    use_extended_neohookean_hyperbola=parse_args->Is_Value_Set("-use_ext_neo_hyper");    
    use_extended_neohookean_smooth=parse_args->Is_Value_Set("-use_ext_neo_smooth");    
    use_corotated=parse_args->Is_Value_Set("-use_corotated");
    use_corot_blend=parse_args->Is_Value_Set("-use_corot_blend");
    use_corot_quartic=parse_args->Is_Value_Set("-use_corot_quartic");
    dump_sv=parse_args->Is_Value_Set("-dump_sv");
    solids_parameters.use_trapezoidal_rule_for_velocities=!parse_args->Get_Option_Value("-use_be");
    print_matrix=parse_args->Is_Value_Set("-print_matrix");
    stiffness_multiplier=(T)parse_args->Get_Double_Value("-stiffen");
    damping_multiplier=(T)parse_args->Get_Double_Value("-dampen");
    solid_body_collection.print_energy=parse_args->Get_Option_Value("-print_energy");
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=parse_args->Get_Integer_Value("-projection_iterations");
    if(parse_args->Is_Value_Set("-project_nullspace")) solids_parameters.implicit_solve_parameters.project_nullspace_frequency=1;
    solid_body_collection.Print_Residuals(parse_args->Get_Option_Value("-residuals"));
    use_constant_ife=parse_args->Get_Option_Value("-use_constant_ife");
    stretch=(T)parse_args->Get_Double_Value("-stretch");
    primary_contour=parse_args->Get_Option_Value("-plot_contour");
    sigma_range=(T)parse_args->Get_Double_Value("-sigma_range");
    image_size=parse_args->Get_Integer_Value("-image_size");
    poissons_ratio=(T)parse_args->Get_Double_Value("-poissons_ratio");
    scatter_plot=parse_args->Get_Option_Value("-scatter_plot");
    use_contrails=parse_args->Get_Option_Value("-use_contrails");
    if(parse_args->Is_Value_Set("-cutoff")) input_cutoff=(T)parse_args->Get_Double_Value("-cutoff");
    if(parse_args->Is_Value_Set("-efc")) input_efc=(T)parse_args->Get_Double_Value("-efc");
    if(parse_args->Is_Value_Set("-poissons_ratio")) input_poissons_ratio=(T)parse_args->Get_Double_Value("-poissons_ratio");
    if(parse_args->Is_Value_Set("-youngs_modulus")) input_youngs_modulus=(T)parse_args->Get_Double_Value("-youngs_modulus");
    test_model_only=parse_args->Get_Option_Value("-test_model_only");
    plot_energy_density=parse_args->Get_Option_Value("-plot_energy_density");
    solids_parameters.implicit_solve_parameters.test_system=parse_args->Is_Value_Set("-test_system");

    switch(test_number){
    case 20: case 21: case 26: 
	    mattress_grid=GRID<TV>(40,8,(T)-2,(T)2,(T)-.4,(T).4);
	break;
        case 22: case 23: case 24: case 25: case 27: case 30:
	    mattress_grid=GRID<TV>(parameter,parameter,(T)-.9,(T).9,(T)-.9,(T).9);
	break;
        case 28: 
            mattress_grid=GRID<TV>(80,16,(T)-2,(T)2,(T)-.4,(T).4);
            break;
    	default:
            mattress_grid=GRID<TV>(20,10,(T)-1,(T)1,(T)-.5,(T).5);
            break;
    }

    switch(test_number){
        case 24:
        case 25:
        case 26:
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=900;
            solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
            last_frame=1500;
            break;
        case 8: 
        case 13:
        case 14:
        case 16:
        case 22:
	case 21:
        case 17:
        case 18:
        case 19:
        case 29:
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            last_frame=200;
            break;
        case 23: 
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            last_frame=500;
            break;
        case 20: case 28:
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=900;
            solids_parameters.deformable_object_collision_parameters.collide_with_interior=true;
            //solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
            attachment_velocity=TV((T).8,0);
	    last_frame=480;
            break;	 
        case 27: case 270: case 30: case 31: case 32: case 33:
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=900;
            solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
            last_frame=3000;
            if(test_number==33) last_frame=1;
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
        case 24:
        case 25:
        case 26:
        case 27: case 270: case 30:
            tests.Create_Mattress(mattress_grid,true);
            break;
        case 8: 
        case 13:{
            tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,4))));
            tests.Add_Ground();
            break;}
        case 23:{
            tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,200))));
            tests.Add_Ground();
            break;}
        case 14:{
            tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,4))));
            tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(4,4))));
            tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(8,4))));
            tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(12,4))));
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
	    break;}
        case 28:{
            tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,0))));
            RIGID_BODY<TV>& box1=tests.Add_Rigid_Body("circle",.4,(T)0);
            RIGID_BODY<TV>& box2=tests.Add_Rigid_Body("circle",.4,(T)0);
            box1.X()=TV(0,-5);
            box2.X()=TV(0,5);
            box1.is_static=false;
            box2.is_static=false;
            kinematic_id=box1.particle_index;
            kinematic_id2=box2.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(box1.particle_index)=true;
            rigid_body_collection.rigid_body_particle.kinematic(box2.particle_index)=true;
            for (int ind=0; ind <=20; ind++){
                curve.Add_Control_Point(ind+4,FRAME<TV>(TV(-5*sin(2.0*pi*ind/5.0),-5*cos(2.0*pi*ind/5.0))));
                curve2.Add_Control_Point(ind+4,FRAME<TV>(TV(5*sin(2.0*pi*ind/5.0),5*cos(2.0*pi*ind/5.0))));
            }
            break;}            
        case 29: last_frame=1;
            tests.Add_Analytic_Smooth_Gear(TV(1,.1),20,16);
            break;
        case 31:{
            TRIANGULATED_AREA<T>* ta=TRIANGULATED_AREA<T>::Create(particles);
            int a=particles.array_collection->Add_Element();
            int b=particles.array_collection->Add_Element();
            int c=particles.array_collection->Add_Element();
            particles.X(a)=TV(1,0);
            particles.X(b)=TV(0,.5);
            particles.X(c)=TV(0,-.5);
            particles.mass(a)=1;
            particles.mass(b)=1;
            particles.mass(c)=1;
            ta->mesh.elements.Append(VECTOR<int,3>(a,b,c));
            solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(ta);
            break;}
        case 32:{
            TRIANGULATED_AREA<T>* ta=TRIANGULATED_AREA<T>::Create(particles);
            int a=particles.array_collection->Add_Element();
            int b=particles.array_collection->Add_Element();
            int c=particles.array_collection->Add_Element();
            particles.X(a)=TV(1,0);
            particles.X(b)=TV(0,1);
            particles.X(c)=TV(0,0);
            particles.mass(a)=1;
            particles.mass(b)=1;
            particles.mass(c)=1;
            ta->mesh.elements.Append(VECTOR<int,3>(a,b,c));
            solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(ta);
            break;}
        case 33:{
            TRIANGULATED_AREA<T>* ta=TRIANGULATED_AREA<T>::Create(particles);
            ta->particles.array_collection->Add_Elements(7);
            for(int i=1;i<=7;i++) ta->particles.X(i)=TV((T).5*(i-4),(T).5*sqrt(3)*(i%2));
            for(int i=1;i<=5;i++) ta->mesh.elements.Append(VECTOR<int,3>(i,i+1,i+2));
            ta->Update_Number_Nodes();
            ta->mesh.Make_Orientations_Consistent();
            solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(ta);
            particles.mass.Fill(1);
            break;}
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
        case 24:
        case 25:
        case 26:
        case 27: case 270: case 30: case 31: case 32: case 33:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(triangulated_area,(T)1e2,poissons_ratio,(T).05);
            if(test_number==32) particles.X(1).x=stretch;
            break;}
        case 8: 
        case 16:
        case 23:
        case 13:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true));
            Add_Constitutive_Model(triangulated_area,(T)1e4,poissons_ratio,(T).01);
            if(test_number==13){RANDOM_NUMBERS<T> rand;rand.Fill_Uniform(particles.X,-1,1);}
            break;}
        case 22:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(triangulated_area,(T)1e4,poissons_ratio,(T).01);
            break;}
        case 14:{
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true));
            
            TRIANGULATED_AREA<T>& triangulated_area1=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>(1);
            solid_body_collection.Add_Force(Create_Finite_Volume(triangulated_area1,new COROTATED_QUARTIC<T,2>((T)2e4,poissons_ratio,(T).01)));

            TRIANGULATED_AREA<T>& triangulated_area2 = solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>(2);
            solid_body_collection.Add_Force(Create_Finite_Volume(triangulated_area2,new NEO_HOOKEAN_COROTATED_BLEND<T,2>((T)2e4,poissons_ratio,(T).01)));

            TRIANGULATED_AREA<T>& triangulated_area3 = solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>(3);
            solid_body_collection.Add_Force(Create_Finite_Volume(triangulated_area3,new NEO_HOOKEAN_EXTRAPOLATED<T,2>((T)2e4,poissons_ratio,(T).01)));

            TRIANGULATED_AREA<T>& triangulated_area4 = solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>(4);
            solid_body_collection.Add_Force(Create_Finite_Volume(triangulated_area4,new NEO_HOOKEAN<T,2>((T)2e4,poissons_ratio,(T).01)));
            break;}
        case 17:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(triangulated_area,(T)1e4,poissons_ratio,(T).01);

            RANDOM_NUMBERS<T> rand;
            rand.Fill_Uniform(particles.X,-1,1);
            break;}
        case 18:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(triangulated_area,(T)1e4,poissons_ratio,(T).01);
            RANDOM_NUMBERS<T> rand;
            rand.Fill_Uniform(particles.X,0,0);
            break;}
        case 19:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(triangulated_area,(T)1e4,poissons_ratio,(T).01);
            break;}
        case 20:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>(1);
            Add_Constitutive_Model(triangulated_area,(T)1e5,poissons_ratio,(T).01);
            break;}
        case 21:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>(1);
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true));
            Add_Constitutive_Model(triangulated_area,(T)1e4,poissons_ratio,(T).01);
            break;}
        case 28:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_AREA<T>&>(1);
            Add_Constitutive_Model(triangulated_area,(T)1e5,poissons_ratio,(T).01);
            break;}
        case 29: break;
        default:
            LOG::cerr<<"Missing implementation for test number "<<test_number<<std::endl;exit(1);}

    if(solid_body_collection.deformable_body_collection.mpi_solids)
        solid_body_collection.deformable_body_collection.mpi_solids->Simple_Partition(solid_body_collection.deformable_body_collection,solid_body_collection.rigid_body_collection.rigid_geometry_collection,particles.X,VECTOR<int,2>(2,1));
    solid_body_collection.Update_Simulated_Particles();

    if(parse_args->Is_Value_Set("-solver_iterations")) solids_parameters.implicit_solve_parameters.cg_iterations=parse_args->Get_Integer_Value("-solver_iterations");
    if(parse_args->Is_Value_Set("-cgsolids")) solids_parameters.implicit_solve_parameters.cg_tolerance=(T)parse_args->Get_Double_Value("-cgsolids");

    if(!semi_implicit) for(int i=1;i<=solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->use_implicit_velocity_independent_forces=true;
    if(!semi_implicit) for(int i=1;i<=solid_body_collection.rigid_body_collection.rigids_forces.m;i++)
        solid_body_collection.rigid_body_collection.rigids_forces(i)->use_implicit_velocity_independent_forces=true;
    if(!semi_implicit) for(int i=1;i<=solid_body_collection.deformable_body_collection.deformables_forces.m;i++) solid_body_collection.deformable_body_collection.deformables_forces(i)->use_implicit_velocity_independent_forces=true;

    if(scatter_plot) Init_Scatter_Plot();

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Initialize_Bodies();
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id)
{
    if(id==kinematic_id) frame=curve.Value(time);
    if(id==kinematic_id2) frame=curve2.Value(time);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id)
{
    if(id==kinematic_id) twist=curve.Derivative(time);
    if(id==kinematic_id2) twist=curve2.Derivative(time);
    return false;
}
//#####################################################################
// Function Set_Particle_Is_Simulated
//#####################################################################
void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE
{
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(test_number==20){
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
        T final_time=40;
        TV velocity_rig=velocity_time<final_time?TV(-velocity,(T)0):TV();
        TV velocity_lef=velocity_time<final_time?TV(velocity,(T)0):TV();
        for(int j=1;j<=n;j++){V(1+m*(j-1))=velocity_lef;V(m+m*(j-1))=velocity_rig;}}
    if(test_number==30){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        T velocity=.1;
        T final_time=80;
        TV velocity_rig=velocity_time<final_time?TV(velocity,(T)0):TV();
        TV velocity_lef=velocity_time<final_time?TV(-velocity,(T)0):TV();
        for(int j=1;j<=n;j++){V(1+m*(j-1))=velocity_lef;V(m+m*(j-1))=velocity_rig;}}
    if(test_number==28){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        TV velocity=velocity_time<5.0?attachment_velocity:TV();
        for(int j=1;j<=n;j++){V(1+m*(j-1))=-velocity;V(m+m*(j-1))=velocity;}}
    if(test_number==31){V(1)=TV(1,0);V(2).x=0;V(3).x=0;}
    if(test_number==32){V(1)=V(3)=TV(0,0);V(3).x=0;}
    if(test_number==33){V(1)=V(2)=V(4)=V(6)=V(7)=TV();}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(test_number==20){
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
    if(test_number==30){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        for(int j=1;j<=n;j++){V(1+m*(j-1))=TV();V(m+m*(j-1))=TV();}}
    if(test_number==28){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        for(int j=1;j<=n;j++) V(1+m*(j-1))=V(m+m*(j-1))=TV();}
    if(test_number==31){V(1)=TV();V(2).x=0;V(3).x=0;}
    if(test_number==32){V(1)=V(3)=TV();V(3).x=0;}
    if(test_number==33){V(1)=V(2)=V(4)=V(6)=V(7)=TV();}
}
void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {
    /*if(test_number==270){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;        
    }*/
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
    if(test_number==29)
    {
        RANDOM_NUMBERS<T> random;
        SMOOTH_GEAR<TV> gear(1,.1,16);
        
        if(0)
        for(int i=1;i<=100000;i++)
        {
            TV X;
            random.Fill_Uniform(X,-2,2);
//            X.Normalize();
            // T sd=gear.Signed_Distance(X);
            TV Y=gear.Surface(X);
            TV N=gear.Normal(X);
            // VECTOR<T,3> col;
            // random.Fill_Uniform(col,0,1);
            Add_Debug_Particle(X, VECTOR<T,3>(gear.Inside(X,.05),gear.Outside(X,.05),gear.Boundary(X,.05)));
//            Add_Debug_Particle(Y, VECTOR<T,3>(1,0,0));
//            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,N);
        }
    }
    if(test_number==33) Plot_Energy_Landscape();
    if(test_number==33 && frame==2)
    {
        PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
        particles.X(5)=TV(-.08,-.16);
        particles.X(3)=particles.X(5);
        particles.X(3).x*=-1;
    }
}
//#####################################################################
// Function Add_Constitutive_Model
//#####################################################################
void Add_Constitutive_Model(TRIANGULATED_AREA<T>& triangulated_area,T stiffness,T poissons_ratio,T damping, T cutoff = 0.4, T efc = 20)
{
    if(input_efc) efc=input_efc;
    if(input_cutoff) cutoff=input_cutoff;
    if(input_poissons_ratio!=-1) poissons_ratio=input_poissons_ratio;
    if(input_youngs_modulus!=0) stiffness=input_youngs_modulus;

    ISOTROPIC_CONSTITUTIVE_MODEL<T,2>* icm=0;
    if(use_extended_neohookean) icm=new NEO_HOOKEAN_EXTRAPOLATED<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_extended_neohookean2) icm=new NEO_HOOKEAN_EXTRAPOLATED2<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_extended_neohookean3) icm=new GENERAL_EXTRAPOLATED<T,2>(*new GEN_NEO_HOOKEAN_ENERGY<T>,stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_int_j_neo) icm=new GENERAL_EXTRAPOLATED<T,2>(*new NEO_J_INTERP_ENERGY<T>,stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_rc_ext) icm=new RC_EXTRAPOLATED<T,2>(*new GEN_NEO_HOOKEAN_ENERGY<T>,stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_rc2_ext) icm=new RC2_EXTRAPOLATED<T,2>(*new GEN_NEO_HOOKEAN_ENERGY<T>,stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_extended_neohookean_refined) icm=new NEO_HOOKEAN_EXTRAPOLATED_REFINED<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,.6,efc);
    else if(use_extended_neohookean_hyperbola) icm=new NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff);
    else if(use_extended_neohookean_smooth) icm=new NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else if(use_corotated) icm=new COROTATED<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else if(use_corot_blend) icm=new NEO_HOOKEAN_COROTATED_BLEND<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else if(use_corot_quartic) icm=new COROTATED_QUARTIC<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else{
        NEO_HOOKEAN<T,2>* nh=new NEO_HOOKEAN<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
        icm=nh;

        nh->use_constant_ife=use_constant_ife;}
    //std::cout << "Lambda= " << icm->constant_lambda << std::cout;
    //std::cout << "Mu    = " << icm->constant_mu << std::cout;
    solid_body_collection.Add_Force(Create_Finite_Volume(triangulated_area,icm));

    if(plot_energy_density) Plot_Energy_Density(icm,stiffness);
    if(primary_contour) Primary_Contour(*icm);
    if(scatter_plot) Add_Primary_Contour_Segments(*icm);
    if(test_model_only) Test_Model(*icm);
}
//#####################################################################
// Function Test_Model_Helper
//#####################################################################
void Test_Model_Helper(const char* str,T a0, T a1, TV da0, TV da1, TV df, T e)
{
    T av=TV::Dot_Product(da1+da0,df)/2/e;
    T dif=(a1-a0)/e;
    char buff[1000];
    sprintf(buff, "============ test ============ %s %8.5f %8.5f (%8.5f)\n", str, av, dif, fabs(av-dif));
    LOG::cout<<buff;
}
//#####################################################################
// Function Test_Model_Helper
//#####################################################################
void Test_Model_Helper(const char* str,TV a0, TV a1, const MATRIX<T,2>& da0, const MATRIX<T,2>& da1, TV df, T e)
{
    TV av=(da1+da0)*df/2/e;
    TV dif=(a1-a0)/e;
    char buff[1000];
    sprintf(buff, "============ test ============ %s %8.5f %8.5f (%8.5f)\n", str, av.Magnitude(), dif.Magnitude(), (av-dif).Magnitude());
    LOG::cout<<buff;
}
//#####################################################################
// Function Test_Model_Helper
//#####################################################################
void Test_Model_Helper(const char* str,const MATRIX<T,2>& a0, const MATRIX<T,2>& a1, const VECTOR<SYMMETRIC_MATRIX<T,2>,2>& da0, const VECTOR<SYMMETRIC_MATRIX<T,2>,2>& da1, TV df, T e)
{
    for(int i=1;i<=TV::m;i++){
        TV av=(da1(i)+da0(i))*df/2/e;
        TV dif=(a1.Transposed().Column(i)-a0.Transposed().Column(i))/e;
        char buff[1000];
        sprintf(buff, "============ test ============ %s %8.5f %8.5f (%8.5f)\n", str, av.Magnitude(), dif.Magnitude(), (av-dif).Magnitude());
        LOG::cout<<buff;}
}
//#####################################################################
// Function Test_Model_Helper
//#####################################################################
template<class RC> void 
Test_Model_Helper(ISOTROPIC_CONSTITUTIVE_MODEL<T,2>* icm, TV &f, TV &df, T e)
{
    RC* rc=dynamic_cast<RC*>(icm);
    if(!rc) return;
    if(f.Min()>0) rc->base.Test(f.x,f.y,0);
    int simplex=0;
    if(f.Product()>rc->extrapolation_cutoff) return;
    typename RC::HELPER h0;
    if(!h0.Compute_E(rc->base,rc->extra_force_coefficient*rc->youngs_modulus,rc->extrapolation_cutoff,f,simplex)) return;
    h0.Compute_dE(rc->base,rc->extra_force_coefficient*rc->youngs_modulus,f,simplex);
    h0.Compute_ddE(rc->base,rc->extra_force_coefficient*rc->youngs_modulus,f,simplex);
    typename RC::HELPER h1;
    if(!h1.Compute_E(rc->base,rc->extra_force_coefficient*rc->youngs_modulus,rc->extrapolation_cutoff,f+df,simplex)) return;
    h1.Compute_dE(rc->base,rc->extra_force_coefficient*rc->youngs_modulus,f+df,simplex);
    h1.Compute_ddE(rc->base,rc->extra_force_coefficient*rc->youngs_modulus,f+df,simplex);
#define XX(k) Test_Model_Helper(#k,h0.k, h1.k, h0.d##k, h1.d##k, df, e);Test_Model_Helper(#k,h0.d##k, h1.d##k, h0.dd##k, h1.dd##k, df, e);
//{TV av=(h1.d##k+h0.d##k)*df/2/e;TV dif=(h1.##k-h0.##k)/e;const char*va=#k;sprintf(buff, "============ test ============ %s %8.5f %8.5f (%8.5f)\n", va, av.Magnitude(), dif.Magnitude(), (av-dif).Magnitude());LOG::cout<<buff;}
    XX(m);
    XX(h);
    XX(phi);
    XX(E);
    XX(z);
    XX(xi);
    XX(s);
//#define YY(k) {for(int i=1;i<=TV::m;i++){TV av=(h1.dd##k(i)+h0.dd##k(i))*df/2/e;TV dif=(h1.d##k.Transposed().Column(i)-h0.d##k.Transposed().Column(i))/e;const char*va=#k;sprintf(buff, "============ test ============ %s %8.5f %8.5f (%8.5f)\n", va, av.Magnitude(), dif.Magnitude(), (av-dif).Magnitude());LOG::cout<<buff;}}
    XX(Q);
    XX(u);
    XX(g);
 }
//#####################################################################
// Function Test_Model
//#####################################################################
void Test_Model(ISOTROPIC_CONSTITUTIVE_MODEL<T,2>& icm)
{
    RANDOM_NUMBERS<T> random;
    for(int i=1;i<=20;i++){
        TV f;
        random.Fill_Uniform(f,0,2);
        T e=1e-5;
        TV df;
        random.Fill_Uniform(df,-e,e);
        f=f.Sorted().Reversed();
        if(random.Get_Uniform_Integer(0,1)==1) f(2)=-f(2);
        LOG::cout<<f<<std::endl;
        icm.Test(DIAGONAL_MATRIX<T,2>(f),1);
        Test_Model_Helper<RC_EXTRAPOLATED<T,2> >(&icm,f,df,e);
        Test_Model_Helper<RC2_EXTRAPOLATED<T,2> >(&icm,f,df,e);
    }
    exit(0);
}
//#####################################################################
// Function Primary_Contour
//#####################################################################
void Primary_Contour(ISOTROPIC_CONSTITUTIVE_MODEL<T,2>& icm)
{
    ARRAY<VECTOR<T,3>,VECTOR<int,2> > img(1,image_size,1,image_size);
    for(int i=1;i<=image_size;i++)
        for(int j=1;j<=image_size;j++){
            T x=(2*i-image_size)*sigma_range/image_size+1e-5;
            T y=(2*j-image_size)*sigma_range/image_size;
            TV g=icm.P_From_Strain(DIAGONAL_MATRIX<T,2>(x,y),1,1).To_Vector();
            DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2> disd;
            icm.Isotropic_Stress_Derivative(DIAGONAL_MATRIX<T,2>(x,y),disd,1);
            SYMMETRIC_MATRIX<T,2> H(disd.x1111,disd.x2211,disd.x2222);
            DIAGONAL_MATRIX<T,2> ev;
            MATRIX<T,2> eigenvectors;
            H.Fast_Solve_Eigenproblem(ev,eigenvectors);
            TV evec=fabs(ev.x11)>fabs(ev.x22)?eigenvectors.Column(1):eigenvectors.Column(2);
            if(evec.Sum()<0) evec=-evec;
            T val=TV::Dot_Product(evec,g);
            img(VECTOR<int,2>(i,j))=VECTOR<T,3>(val<0,val>=0,0);
        }

    for(int i=1;i<=image_size;i++){if(i%10>0 && i%10<5) img(VECTOR<int,2>(i,image_size/2))=img(VECTOR<int,2>(image_size/2,i))=VECTOR<T,3>(0,0,1);}

    PPM_FILE<T>::Write("primary_contour.ppm", img);
}
//#####################################################################
// Function Init_Scatter_Plot
//#####################################################################
void Init_Scatter_Plot()
{
    int s=0;
    for(int f=1;FINITE_VOLUME<TV,2>* force=solid_body_collection.deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,2>*>(f);f++){
        force->Update_Position_Based_State(0,true);
        s+=force->Fe_hat.m;}

    contrail.Resize(s);
    contrail_colors.Resize(s);
    RANDOM_NUMBERS<T> rand;
    rand.Fill_Uniform(contrail_colors,0,1);
}
//#####################################################################
// Function Update_Scatter_Plot
//#####################################################################
void Update_Scatter_Plot()
{
    for(int f=1,k=1;FINITE_VOLUME<TV,2>* force=solid_body_collection.deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,2>*>(f);f++){
        for(int i=1;i<=force->Fe_hat.m;i++){
            if(!use_contrails) contrail(k).Remove_All();
            contrail(k++).Append(force->Fe_hat(i).To_Vector());}}
}
//#####################################################################
// Function Dump_Scatter_Plot
//#####################################################################
void Dump_Scatter_Plot(int frame)
{
    RANGE<TV> box(-image_size,image_size,-image_size,image_size);
    char buff[1000];
    sprintf(buff, "svd-plot-%04d.eps", frame);
    EPS_FILE_GEOMETRY<T> eps(buff,box);
    eps.Set_Point_Size(4*sigma_range/image_size);
    eps.fixed_bounding_box=true;
    eps.bounding_box=RANGE<TV>(-sigma_range,sigma_range,-sigma_range,sigma_range);
    for(int i=1;i<=contrail.m;i++){
        eps.Line_Color(contrail_colors(i));
        for(int j=2;j<=contrail(i).m;j++)
            eps.Draw_Line(contrail(i)(j-1),contrail(i)(j));}
    for(int i=1;i<=contrail.m;i++){
        eps.Line_Color(contrail_colors(i));
        eps.Draw_Point(contrail(i).Last());}
    eps.Line_Color(VECTOR<T,3>());

    for(int i=1;i<=contour_segments.m;i++)
        eps.Draw_Line(contour_segments(i).x,contour_segments(i).y);

    eps.Draw_Line(TV(0,-image_size),TV(0,image_size));
    eps.Draw_Line(TV(-image_size,0),TV(image_size,0));
}
//#####################################################################
// Function Add_Primary_Contour_Segment
//#####################################################################
T Contour_Crossing(const TV& g0,const TV& v0,const TV& g1,const TV& v1)
{
    T a=TV::Dot_Product(g0,v0);
    T b=TV::Dot_Product(g1,v1);
    if(!a) return 0;
    if(!b) return 1;
    if(TV::Dot_Product(v0,v1)<0) b=-b;
    if((a>0) == (b>0)) return -1;
//    LOG::cout<<"CROSS "<<a<<"  "<<b<<"  "<<g0<<"  "<<g1<<"  "<<v0<<"  "<<v1<<std::endl;
    return a/(a-b);
}
//#####################################################################
// Function Add_Primary_Contour_Segments
//#####################################################################
void Add_Primary_Contour_Segments(ISOTROPIC_CONSTITUTIVE_MODEL<T,2>& icm)
{
    ARRAY<TV,VECTOR<int,2> > evec(1,image_size,1,image_size);
    ARRAY<TV,VECTOR<int,2> > grad(1,image_size,1,image_size);
    for(int i=1;i<=image_size;i++)
        for(int j=1;j<=image_size;j++){
            T x=(2*i-image_size)*sigma_range/image_size+1e-5;
            T y=(2*j-image_size)*sigma_range/image_size;
            TV g=icm.P_From_Strain(DIAGONAL_MATRIX<T,2>(x,y),1,1).To_Vector();
            DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2> disd;
            icm.Isotropic_Stress_Derivative(DIAGONAL_MATRIX<T,2>(x,y),disd,1);
            SYMMETRIC_MATRIX<T,2> H(disd.x1111,disd.x2211,disd.x2222);
            DIAGONAL_MATRIX<T,2> ev;
            MATRIX<T,2> eigenvectors;
            H.Fast_Solve_Eigenproblem(ev,eigenvectors);
            evec(VECTOR<int,2>(i,j))=fabs(ev.x11)>fabs(ev.x22)?eigenvectors.Column(1):eigenvectors.Column(2);
            grad(VECTOR<int,2>(i,j))=g;}

    for(int i=1;i<image_size;i++)
        for(int j=1;j<image_size;j++){
            TV g00=grad(VECTOR<int,2>(i,j)),g01=grad(VECTOR<int,2>(i,j+1)),g10=grad(VECTOR<int,2>(i+1,j)),g11=grad(VECTOR<int,2>(i+1,j+1));
            TV v00=evec(VECTOR<int,2>(i,j)),v01=evec(VECTOR<int,2>(i,j+1)),v10=evec(VECTOR<int,2>(i+1,j)),v11=evec(VECTOR<int,2>(i+1,j+1));
            T cx0=Contour_Crossing(g00,v00,g10,v10);
            T cx1=Contour_Crossing(g01,v01,g11,v11);
            T c0x=Contour_Crossing(g00,v00,g01,v01);
            T c1x=Contour_Crossing(g10,v10,g11,v11);
            int n=(cx0>=0)+(cx1>=0)+(c0x>=0)+(c1x>=0);
            if(n<2) continue;
            TV X00((2*i-image_size)*sigma_range/image_size+1e-5,(2*j-image_size)*sigma_range/image_size);
            TV X01((2*i-image_size)*sigma_range/image_size+1e-5,(2*(j+1)-image_size)*sigma_range/image_size);
            TV X10((2*(i+1)-image_size)*sigma_range/image_size+1e-5,(2*j-image_size)*sigma_range/image_size);
            TV X11((2*(i+1)-image_size)*sigma_range/image_size+1e-5,(2*(j+1)-image_size)*sigma_range/image_size);
            TV Yx0=X00+(X10-X00)*cx0,Yx1=X01+(X11-X01)*cx1,Y0x=X00+(X01-X00)*c0x,Y1x=X10+(X11-X10)*c1x;
            if(cx0>=0 && cx1>=0){
                contour_segments.Append(VECTOR<TV,2>(Yx0,Yx1));
                if(n==3 && c0x>=0) contour_segments.Append(VECTOR<TV,2>(Y0x,(Yx0+Yx1)/2));
                if(n==3 && c1x>=0) contour_segments.Append(VECTOR<TV,2>(Y1x,(Yx0+Yx1)/2));}
            if(c0x>=0 && c1x>=0){
                contour_segments.Append(VECTOR<TV,2>(Y0x,Y1x));
                if(n==3 && cx0>=0) contour_segments.Append(VECTOR<TV,2>(Yx0,(Y0x+Y1x)/2));
                if(n==3 && cx1>=0) contour_segments.Append(VECTOR<TV,2>(Yx1,(Y0x+Y1x)/2));}
            if(n>2) continue;
            if(c0x>=0 && cx0>=0) contour_segments.Append(VECTOR<TV,2>(Y0x,Yx0));
            if(c0x>=0 && cx1>=0) contour_segments.Append(VECTOR<TV,2>(Y0x,Yx1));
            if(c1x>=0 && cx0>=0) contour_segments.Append(VECTOR<TV,2>(Y1x,Yx0));
            if(c1x>=0 && cx1>=0) contour_segments.Append(VECTOR<TV,2>(Y1x,Yx1));}
}
//#####################################################################
// Function Plot_Energy_Landscape
//#####################################################################
void Plot_Energy_Density(ISOTROPIC_CONSTITUTIVE_MODEL<T,2>* icm,T stiffness)
{
    TRIANGULATED_AREA<T>& ta=tests.Create_Mattress(GRID<TV>(2+image_size+1,2+image_size+1,-sigma_range,sigma_range,-sigma_range,sigma_range),true,RIGID_BODY_STATE<TV>());
    TRIANGULATED_SURFACE<T> ts(ta.mesh,*new PARTICLES<VECTOR<T,3> >);
    ts.particles.array_collection->Add_Elements(ta.particles.X.m);
    for(int i=1;i<=ta.particles.X.m;i++){
        TV X=ta.particles.X(i);
        X.x+=1.1e-4;
        X.y+=2.3e-4;
        TV Z=X;
        if((fabs(Z.x)>fabs(Z.y)?Z.x:Z.y)<0) Z=-Z;
        VECTOR<T,3> Y(X.x,X.y,icm->P_From_Strain(DIAGONAL_MATRIX<T,2>(Z),1,0).x11/stiffness);
        ts.particles.X(i)=Y;}
    FILE_UTILITIES::Write_To_File(this->stream_type,"surface.tri",ts);
}
//#####################################################################
// Function Plot_Energy_Landscape
//#####################################################################
void Plot_Energy_Landscape()
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    ARRAY<TV> F(particles.X.m);
    ARRAY<TWIST<TV> > TW;
    ARRAY<TV> X;
    ARRAY<TV> N;
    ARRAY<T> E;
    ARRAY<TV> TX(particles.X);
    for(int i=-image_size/2;i<=image_size/2;i++)
        for(int j=-image_size/2;j<=image_size/2;j++){
            F.Fill(TV());
            T x=2*sigma_range*i/image_size+1.1e-5;
            T y=2*sigma_range*j/image_size+1.2e-5;
            particles.X(3)=TV(-x,y);
            particles.X(5)=TV(x,y);
            solid_body_collection.Update_Position_Based_State(0,false);
            T ke,pe;
            solid_body_collection.Compute_Energy(0,ke,pe);
            solid_body_collection.Add_Velocity_Independent_Forces(F,TW,0);
            E.Append(ke+pe);
            X.Append(particles.X(5));
            N.Append(F(5).Normalized());
        }

    INTERPOLATED_COLOR_MAP<T> mp;
    mp.Initialize_Colors(E.Min(),E.Max(),false,true,false);

    for(int i=1;i<=X.m;i++){
        Add_Debug_Particle(X(i),mp(E(i)));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,N(i));}
    particles.X=TX;
}
//#####################################################################
};
}

#endif

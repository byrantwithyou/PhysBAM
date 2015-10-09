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
//  100. Primary contour field
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <Tools/Images/PPM_FILE.h>
#include <Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/constants.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Geometry/Basic_Geometry/SMOOTH_GEAR.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Images/EPS_FILE.h>
#include <Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Deformables/Bindings/LINEAR_BINDING.h>
#include <Deformables/Bindings/SOFT_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <Deformables/Constitutive_Models/COROTATED.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/GEN_NEO_HOOKEAN_ENERGY.h>
#include <Deformables/Constitutive_Models/GENERAL_EXTRAPOLATED.h>
#include <Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED.h>
#include <Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED2.h>
#include <Deformables/Constitutive_Models/RC_EXTRAPOLATED.h>
#include <Deformables/Constitutive_Models/RC2_EXTRAPOLATED.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/FINITE_VOLUME.h>
#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
#include <Solids/Forces_And_Torques/ETHER_DRAG.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <fstream>
namespace PhysBAM{
template<class TV> void Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color);
template<class TV,class ATTR> void Debug_Particle_Set_Attribute(ATTRIBUTE_ID id,const ATTR& attr);

template<class T_input>
class STANDARD_TESTS:public SOLIDS_EXAMPLE<VECTOR<T_input,2> >
{
    typedef T_input T;
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    typedef SOLIDS_EXAMPLE<TV> BASE;
    using BASE::solids_parameters;using BASE::output_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::solid_body_collection;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions;using BASE::solids_evolution;
    using BASE::test_number;using BASE::data_directory;

    std::ofstream svout;

    SOLIDS_STANDARD_TESTS<TV> tests;

    GRID<TV> mattress_grid;
    TV attachment_velocity;
    bool semi_implicit;
    bool test_forces;
    bool use_extended_neohookean;
    bool use_extended_neohookean2;
    bool use_extended_neohookean3;
    bool use_rc_ext;
    bool use_rc2_ext;
    bool use_extended_neohookean_refined;
    bool use_corotated;
    bool use_corotated_fixed;
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
    GEOMETRY_PARTICLES<VECTOR<T,3> > energy_particles;
    TRIANGULATED_SURFACE<T>* energy_mesh;
    std::string dual_directory;
    bool energy_profile_plot;
    T energy_profile_plot_min,energy_profile_plot_range;
    T plot_scale;
    T ether_drag;
    bool project_nullspace,opt_residuals;

    STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :BASE(stream_type_input,parse_args),tests(stream_type_input,data_directory,solid_body_collection),semi_implicit(false),test_forces(false),use_extended_neohookean(false),
        use_extended_neohookean2(false),use_extended_neohookean3(false),use_rc_ext(false),use_rc2_ext(false),use_extended_neohookean_refined(false),
        use_corotated(false),use_corotated_fixed(false),
        dump_sv(false),kinematic_id(0),kinematic_id2(0),print_matrix(false),parameter(0),stiffness_multiplier(1),damping_multiplier(1),use_constant_ife(false),
        stretch(1),primary_contour(false),sigma_range(3),image_size(500),poissons_ratio((T).45),scatter_plot(false),use_contrails(false),input_cutoff(FLT_MAX),input_efc(FLT_MAX),
        input_poissons_ratio(-1),input_youngs_modulus(0),test_model_only(false),plot_energy_density(false),energy_mesh(0),energy_profile_plot(false),energy_profile_plot_min(0),
        energy_profile_plot_range(0),plot_scale(1),ether_drag(0),project_nullspace(false),opt_residuals(false)
    {
        solids_parameters.implicit_solve_parameters.cg_projection_iterations=5;
        parse_args.Add("-semi_implicit",&semi_implicit,"use semi implicit forces");
        parse_args.Add("-test_forces",&test_forces,"use fully implicit forces");
        parse_args.Add("-use_ext_neo",&use_extended_neohookean,"use_extended_neohookean");
        parse_args.Add("-use_ext_neo2",&use_extended_neohookean2,"use_extended_neohookean2");
        parse_args.Add("-use_ext_neo3",&use_extended_neohookean3,"use_extended_neohookean3");
        parse_args.Add("-use_rc_ext",&use_rc_ext,"use_rc_ext");
        parse_args.Add("-use_rc2_ext",&use_rc2_ext,"use_rc2_ext");
        parse_args.Add("-use_ext_neo_ref",&use_extended_neohookean_refined,"use_extended_neohookean_refined");
        parse_args.Add("-use_corotated",&use_corotated,"use_corotated");
        parse_args.Add("-use_corotated_fixed",&use_corotated_fixed,"use_corotated_fixed");
        parse_args.Add("-dump_sv",&dump_sv,"dump_sv");
        parse_args.Add("-parameter",&parameter,"value","parameter used by multiple tests to change the parameters of the test");
        parse_args.Add("-stiffen",&stiffness_multiplier,"","stiffness multiplier for various tests");
        parse_args.Add("-dampen",&damping_multiplier,"","damping multiplier for various tests");
        parse_args.Add("-residuals",&opt_residuals,"print residuals during timestepping");
        parse_args.Add("-print_energy",&solid_body_collection.print_energy,"print energy statistics");
        parse_args.Add("-cgsolids",&solids_parameters.implicit_solve_parameters.cg_tolerance,"value","CG tolerance for backward Euler");
        parse_args.Add_Not("-use_be",&solids_parameters.use_trapezoidal_rule_for_velocities,"use backward euler");
        parse_args.Add("-print_matrix",&print_matrix,"print_matrix");
        parse_args.Add("-project_nullspace",&project_nullspace,"project out nullspace");
        parse_args.Add("-projection_iterations",&solids_parameters.implicit_solve_parameters.cg_projection_iterations,"value","number of iterations used for projection in cg");
        parse_args.Add("-solver_iterations",&solids_parameters.implicit_solve_parameters.cg_iterations,"value","number of iterations used for solids system");
        parse_args.Add("-use_constant_ife",&use_constant_ife,"use constant extrapolation on inverting finite element fix");
        parse_args.Add("-stretch",&stretch,"value","stretch");
        parse_args.Add("-plot_contour",&primary_contour,"value","plot primary contour");
        parse_args.Add("-energy_profile_plot",&energy_profile_plot,"plot energy profiles");
        parse_args.Add("-image_size",&image_size,"value","image size for plots");
        parse_args.Add("-sigma_range",&sigma_range,"value","sigma range for plots");
        parse_args.Add("-poissons_ratio",&poissons_ratio,"value","poisson's ratio");
        parse_args.Add("-scatter_plot",&scatter_plot,"Create contrail plot with singular values");
        parse_args.Add("-use_contrails",&use_contrails,"Show contrails in plot");
        parse_args.Add("-cutoff",&input_cutoff,"value","cutoff");
        parse_args.Add("-efc",&input_efc,"value","efc");
        parse_args.Add("-poissons_ratio",&input_poissons_ratio,"value","poissons_ratio");
        parse_args.Add("-youngs_modulus",&input_youngs_modulus,"value","youngs modulus, only for test 41 so far");
        parse_args.Add("-test_model_only",&test_model_only,"test_model_only");
        parse_args.Add("-test_system",&solids_parameters.implicit_solve_parameters.test_system,"test_system");
        parse_args.Add("-plot_scale",&plot_scale,"value","Scale height of energy plot");
        parse_args.Add("-ether_drag",&ether_drag,"value","Ether drag");
        parse_args.Add("-plot_energy_density",&plot_energy_density,"plot_energy_density");
        parse_args.Parse();

        tests.data_directory=data_directory;
        LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
        output_directory=LOG::sprintf("Standard_Tests/Test_%d",test_number);
        last_frame=1000;
        if(project_nullspace) solids_parameters.implicit_solve_parameters.project_nullspace_frequency=1;
        solid_body_collection.Print_Residuals(opt_residuals);
    
        solids_parameters.triangle_collision_parameters.perform_self_collision=false;

        switch(test_number){
            case 20: case 21: case 26: 
                mattress_grid=GRID<TV>(TV_INT(40,8),RANGE<TV>(TV((T)-2,(T)-.4),TV((T)2,(T).4)));
                break;
            case 22: case 23: case 24: case 25: case 27: case 30:
                mattress_grid=GRID<TV>(TV_INT(parameter?parameter+1:11,parameter?parameter+1:11),RANGE<TV>(TV((T)-.9,(T)-.9),TV((T).9,(T).9)));
                break;
            case 28: 
                mattress_grid=GRID<TV>(TV_INT(80,16),RANGE<TV>(TV((T)-2,(T)-.4),TV((T)2,(T).4)));
                break;
            default:
                mattress_grid=GRID<TV>(TV_INT(parameter?parameter+1:11,parameter?parameter+1:11),RANGE<TV>(TV((T)-1,(T)-1),TV((T)1,(T)1)));
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
            case 34:
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
            case 27: case 270: case 30: case 31: case 32: case 33: case 100:
                solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
                solids_parameters.implicit_solve_parameters.cg_iterations=900;
                solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
                last_frame=2000;
                if(test_number==33) last_frame=1;
                if(test_number==31) last_frame=250;
                break;
            default:
                LOG::cerr<<"Unrecognized test number "<<test_number<<std::endl;exit(1);}

        output_directory=LOG::sprintf("Standard_Tests/Test_%d",test_number);
    }

    // Unused callbacks
    void Postprocess_Frame(const int frame) override
    {
        if(dump_sv) svout.close();
        if(scatter_plot) Dump_Scatter_Plot(frame);
        if(energy_profile_plot) Energy_Profile_Plot(frame);

        if(test_number==34){
            ARRAY<int> num_elements(solid_body_collection.deformable_body_collection.particles.X.m);
            ARRAY<DIAGONAL_MATRIX<T,TV::m> > F_ave(solid_body_collection.deformable_body_collection.particles.X.m);
            ARRAY<DIAGONAL_MATRIX<T,TV::m> > P_ave(solid_body_collection.deformable_body_collection.particles.X.m);
            FINITE_VOLUME<TV,2>& fv=solid_body_collection.deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,2>&>();
            for(int i=0;i<fv.strain_measure.mesh.elements.m;i++){
                num_elements.Subset(fv.strain_measure.mesh.elements(i))+=1;
                F_ave.Subset(fv.strain_measure.mesh.elements(i))+=fv.Fe_hat(i);
                P_ave.Subset(fv.strain_measure.mesh.elements(i))+=fv.isotropic_model->P_From_Strain(fv.Fe_hat(i),i)*fv.Be_scales(i);}

            ARRAY<T> F_norms(num_elements.m),P_norms(num_elements.m);
            for(int i=0;i<num_elements.m;i++)
                if(num_elements(i)){
                    F_ave(i)/=num_elements(i);
                    P_ave(i)/=num_elements(i);
                    F_norms(i)=(F_ave(i)-1).Frobenius_Norm();
                    P_norms(i)=(P_ave(i)-1).Frobenius_Norm();}

            INTERPOLATED_COLOR_MAP<T> F_map,P_map;
            F_map.Initialize_Colors(0,3,false,true,false);
            P_map.Initialize_Colors(0,500,false,true,false);
            printf("Frame MAX %g %g\n",F_norms.Max(),P_norms.Max());

            typedef VECTOR<T,3> COL;
            ARRAY<COL> F_colors(num_elements.m),P_colors(num_elements.m);
            for(int i=0;i<num_elements.m;i++){
                F_colors(i)=F_map(F_norms(i));
                P_colors(i)=P_map(P_norms(i));}

            const char* pref[3]={"frame","stress","strain"};

            for(int type=0;type<3;type++){
                char buff[100];
                sprintf(buff,"%s-%03d.eps",pref[type],frame);
                EPS_FILE<T> eps(buff);
                eps.cur_format.line_width=.05;
                eps.Use_Round_Join();
                eps.Use_Round_Cap();
                eps.Bound(TV(-6,-6));
                eps.Bound(TV(13,-21));
                for(int i=0;i<fv.strain_measure.mesh.elements.m;i++){
                    int a,b,c;
                    fv.strain_measure.mesh.elements(i).Get(a,b,c);
                    const TV& A=solid_body_collection.deformable_body_collection.particles.X(a);
                    const TV& B=solid_body_collection.deformable_body_collection.particles.X(b);
                    const TV& C=solid_body_collection.deformable_body_collection.particles.X(c);
                    if(type==1) eps.Emit_Gouraud_Triangle(A,B,C,F_colors(a),F_colors(b),F_colors(c));
                    if(type==2) eps.Emit_Gouraud_Triangle(A,B,C,P_colors(a),P_colors(b),P_colors(c));}
                
                if(!fv.strain_measure.mesh.boundary_mesh) fv.strain_measure.mesh.Initialize_Boundary_Mesh();
                for(int i=0;i<fv.strain_measure.mesh.boundary_mesh->elements.m;i++){
                    const TV& A=solid_body_collection.deformable_body_collection.particles.X(fv.strain_measure.mesh.boundary_mesh->elements(i)(0));
                    const TV& B=solid_body_collection.deformable_body_collection.particles.X(fv.strain_measure.mesh.boundary_mesh->elements(i)(1));
                    eps.Draw_Object(A,B);}}}
    }
    void Postprocess_Solids_Substep(const T time,const int substep) override {}
    void Apply_Constraints(const T dt,const T time) override {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) override {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) override {}
    void Update_Time_Varying_Material_Properties(const T time) override {}
    void Limit_Solids_Dt(T& dt,const T time) override {}
    void Post_Initialization() override {}
    void Preprocess_Substep(const T dt,const T time) override
    {
        if(test_forces){
            solid_body_collection.deformable_body_collection.Test_Energy(time);
            solid_body_collection.deformable_body_collection.Test_Force_Derivatives(time);}
    }
    void Postprocess_Substep(const T dt,const T time) override
    {
        if (dump_sv)
        {
            FINITE_VOLUME<TV,2>& force_field=solid_body_collection.deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,2>&>();
            ARRAY<DIAGONAL_MATRIX<T,2> >& sv=force_field.Fe_hat;
            
            for (int i=0;i<sv.m;i++)
            {
                svout << sv(i).x.x << " " << sv(i).x.y << std::endl;
            }
        }
        if(scatter_plot) Update_Scatter_Plot();
        if(energy_profile_plot){
            FINITE_VOLUME<TV,2>& force_field=solid_body_collection.deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,2>&>();
            for(int i=0;i<force_field.Fe_hat.m;i++) Add_Debug_Particle(force_field.Fe_hat(i).To_Vector(),VECTOR<T,3>(1,1,0));}
    }

    void Align_Deformable_Bodies_With_Rigid_Bodies() override {}
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) override {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) override {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) override {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) override {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) override {}
    //void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) override {}

void After_Initialization() override {BASE::After_Initialization();}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    switch(test_number){
        case 24:
        case 25:
        case 26:
        case 27: case 270: case 30:
            tests.Create_Mattress(mattress_grid,true,0,1000);
            if(test_number==30) contrail.Resize(2*mattress_grid.counts.Product());
            break;
        case 8: 
        case 13:{
            tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,4))),1000);
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
            box1.Frame().t=TV(0,-10);
            box2.Frame().t=TV(0,12);
            box1.is_static=true;
            box2.is_static=false;
            kinematic_id=box2.particle_index;
            rigid_body_collection.rigid_body_particles.kinematic(box2.particle_index)=true;
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
            box1.Frame().t=TV(0,-2);
            box2.Frame().t=TV(0,2);
            box1.is_static=true;
            box2.is_static=false;
            kinematic_id=box2.particle_index;
            rigid_body_collection.rigid_body_particles.kinematic(box2.particle_index)=true;
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
            //box1.Frame().t=TV(0,-10);
            box2.Frame().t=TV(0,12);
            //box1.is_static=true;
            box2.is_static=false;
            kinematic_id=box2.particle_index;
            rigid_body_collection.rigid_body_particles.kinematic(box2.particle_index)=true;
            curve.Add_Control_Point(0,FRAME<TV>(TV(0,12)));
            curve.Add_Control_Point(5,FRAME<TV>(TV(0,0)));
            curve.Add_Control_Point(6,FRAME<TV>(TV(0,0)));
            curve.Add_Control_Point(11,FRAME<TV>(TV(0,12)));
            break;}
        case 21:{
            tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,4))));
            RIGID_BODY<TV>& box1=tests.Add_Rigid_Body("square",1,(T)0);
            box1.Frame().t=TV(0,-6);
            box1.is_static=true;
            break;}
        case 28:{
            tests.Create_Mattress(mattress_grid,true,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,0))));
            RIGID_BODY<TV>& box1=tests.Add_Rigid_Body("circle",.4,(T)0);
            RIGID_BODY<TV>& box2=tests.Add_Rigid_Body("circle",.4,(T)0);
            box1.Frame().t=TV(0,-5);
            box2.Frame().t=TV(0,5);
            box1.is_static=false;
            box2.is_static=false;
            kinematic_id=box1.particle_index;
            kinematic_id2=box2.particle_index;
            rigid_body_collection.rigid_body_particles.kinematic(box1.particle_index)=true;
            rigid_body_collection.rigid_body_particles.kinematic(box2.particle_index)=true;
            for (int ind=0;ind<20;ind++){
                curve.Add_Control_Point(ind+4,FRAME<TV>(TV(-5*sin(2.0*pi*ind/5.0),-5*cos(2.0*pi*ind/5.0))));
                curve2.Add_Control_Point(ind+4,FRAME<TV>(TV(5*sin(2.0*pi*ind/5.0),5*cos(2.0*pi*ind/5.0))));
            }
            break;}            
        case 29: last_frame=1;
            tests.Add_Analytic_Smooth_Gear(TV(1,.1),20,16);
            break;
        case 31:{
            TRIANGULATED_AREA<T>* ta=TRIANGULATED_AREA<T>::Create(particles);
            int a=particles.Add_Element();
            int b=particles.Add_Element();
            int c=particles.Add_Element();
            particles.X(a)=TV(4.1,sqrt(3)/2-1.5);
            particles.X(b)=TV(4.1-.5,-1.5);
            particles.X(c)=TV(4.1+.5,-1.5);
            particles.mass(a)=1;
            particles.mass(b)=1;
            particles.mass(c)=1;
            ta->mesh.elements.Append(VECTOR<int,3>(a,b,c));
            solid_body_collection.deformable_body_collection.Add_Structure(ta);
            contrail_colors.Append(VECTOR<T,3>(1,0,0));
            contrail.Resize(contrail_colors.m);
            break;}
        case 32:{
            TRIANGULATED_AREA<T>* ta=TRIANGULATED_AREA<T>::Create(particles);
            int a=particles.Add_Element();
            int b=particles.Add_Element();
            int c=particles.Add_Element();
            particles.X(a)=TV(1,0);
            particles.X(b)=TV(0,1);
            particles.X(c)=TV(0,0);
            particles.mass(a)=1;
            particles.mass(b)=1;
            particles.mass(c)=1;
            ta->mesh.elements.Append(VECTOR<int,3>(a,b,c));
            solid_body_collection.deformable_body_collection.Add_Structure(ta);
            break;}
        case 33:{
            TRIANGULATED_AREA<T>* ta=TRIANGULATED_AREA<T>::Create(particles);
            ta->particles.Add_Elements(7);
            for(int i=0;i<7;i++) ta->particles.X(i)=TV((T).5*(i-3),(T).5*sqrt(3)*(i%2));
            for(int i=0;i<5;i++) ta->mesh.elements.Append(VECTOR<int,3>(i,i+1,i+2));
            ta->Update_Number_Nodes();
            ta->mesh.Make_Orientations_Consistent();
            solid_body_collection.deformable_body_collection.Add_Structure(ta);
            particles.mass.Fill(1);
            break;}
        case 34: {
            TRIANGULATED_AREA<T>* ta=TESSELLATION::Generate_Triangles(SPHERE<TV>(TV(),1),12);
            for(int i=0;i<ta->particles.X.m;i++){
                T t=atan2(ta->particles.X(i).y,ta->particles.X(i).x);
                T r2=cos(t)+cos(2*t)+3-.1*cos(5*t)+.5*sin(t);
                ta->particles.X(i)*=r2;
                ta->particles.X(i)=ta->particles.X(i).Orthogonal_Vector();}
            ta->particles.X+=TV(0,(T)2.4);
            tests.Copy_And_Add_Structure(*ta,0);
            particles.mass.Fill(1);
            tests.Add_Ground(.5,-20);
            break;}
        case 100:{
            TRIANGULATED_AREA<T>* ta=TRIANGULATED_AREA<T>::Create(particles);
            if(parameter==1){
                particles.Add_Elements(3);
                particles.X(0)=TV(-.5,0)*1.5;
                particles.X(1)=TV(.5,0)*1.5;
                particles.X(2)=TV(0,sqrt(3)/2)*1.5;
                ta->mesh.elements.Append(VECTOR<int,3>(0,1,2));}
            else{
                particles.Add_Elements(9);
                particles.X(0)=TV(-.5,0);
                particles.X(1)=TV(.5,0);
                particles.X(2)=TV(0,sqrt(3)/2);
                particles.X(3)=TV(-.5,0);
                particles.X(4)=TV(.5,0);
                particles.X(5)=TV(0,sqrt(3)/2);
                particles.X(6)=TV(-.5,0);
                particles.X(7)=TV(.5,0);
                particles.X(8)=TV(0,sqrt(3)/2);
                ta->mesh.elements.Append(VECTOR<int,3>(0,1,2));
                ta->mesh.elements.Append(VECTOR<int,3>(3,4,5));
                ta->mesh.elements.Append(VECTOR<int,3>(6,7,8));}
            particles.mass.Fill(1);
            solid_body_collection.deformable_body_collection.Add_Structure(ta);
            contrail_colors.Append(VECTOR<T,3>(1,0,0));
            contrail_colors.Append(VECTOR<T,3>(0,1,0));
            contrail_colors.Append(VECTOR<T,3>(0,0,1));
            contrail.Resize(contrail_colors.m);
            break;}
        default:
            LOG::cerr<<"Missing implementation for test number "<<test_number<<std::endl;exit(1);}

    // add structures and rigid bodies to collisions
    solid_body_collection.deformable_body_collection.collisions.collision_structures.Append_Elements(solid_body_collection.deformable_body_collection.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(solid_body_collection.deformable_body_collection.structures);

    // correct number nodes
    for(int i=0;i<solid_body_collection.deformable_body_collection.structures.m;i++) solid_body_collection.deformable_body_collection.structures(i)->Update_Number_Nodes();

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
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(triangulated_area,(T)1e2,poissons_ratio,(T).05);
            if(test_number==32) particles.X(0).x=stretch;
            break;}
        case 100:{
            for(int i=0;TRIANGULATED_AREA<T>* triangulated_area=solid_body_collection.deformable_body_collection.template Find_Structure<TRIANGULATED_AREA<T>*>(i);i++)
                Add_Constitutive_Model(*triangulated_area,(T)1e2,poissons_ratio,(T).05);
            if(parameter==1) Place_Triangle(1,3,.1,.2,1.2,TV(4.1,0));
            else if(parameter==2){
                Place_Triangle(1,2.6,1.8,1.5,3.34,TV(4.1,2));
                Place_Triangle(2,3,.1,.2,1.7,TV(4.1,0));
                Place_Triangle(3,1.3,.15,.2,-.3,TV(4.1,-2));}
            else{
                Place_Triangle(1,2.6,1.8,1.5,3.34,TV(4.1,2));
                Place_Triangle(2,1.3,.15,.2,-.3,TV(4.1,0));
                Place_Triangle(3,2.7,.9,1.8,.4,TV(4.1,-2));}
            break;}
        case 8:
        case 34:
        case 16:
        case 23:
        case 13:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.template Find_Structure<TRIANGULATED_AREA<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true));
            Add_Constitutive_Model(triangulated_area,(T)1e4,poissons_ratio,(T).01);
            if(test_number==13){RANDOM_NUMBERS<T> rand;rand.Set_Seed(12345);rand.Fill_Uniform(particles.X,-1,1);}
            break;}
        case 22:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(triangulated_area,(T)1e4,poissons_ratio,(T).01);
            break;}
        case 14:{
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true));

            TRIANGULATED_AREA<T>& triangulated_area3=solid_body_collection.deformable_body_collection.template Find_Structure<TRIANGULATED_AREA<T>&>(2);
            solid_body_collection.Add_Force(Create_Finite_Volume(triangulated_area3,new NEO_HOOKEAN_EXTRAPOLATED<T,2>((T)2e4,poissons_ratio,(T).01)));

            TRIANGULATED_AREA<T>& triangulated_area4=solid_body_collection.deformable_body_collection.template Find_Structure<TRIANGULATED_AREA<T>&>(3);
            solid_body_collection.Add_Force(Create_Finite_Volume(triangulated_area4,new NEO_HOOKEAN<T,2>((T)2e4,poissons_ratio,(T).01)));
            break;}
        case 17:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(triangulated_area,(T)1e4,poissons_ratio,(T).01);

            RANDOM_NUMBERS<T> rand;
            rand.Set_Seed(12345);
            rand.Fill_Uniform(particles.X,-1,1);
            break;}
        case 18:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(triangulated_area,(T)1e4,poissons_ratio,(T).01);
            RANDOM_NUMBERS<T> rand;rand.Set_Seed(12345);
            rand.Fill_Uniform(particles.X,0,0);
            break;}
        case 19:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(triangulated_area,(T)1e4,poissons_ratio,(T).01);
            break;}
        case 20:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(triangulated_area,(T)1e5,poissons_ratio,(T).01);
            break;}
        case 21:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.template Find_Structure<TRIANGULATED_AREA<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true));
            Add_Constitutive_Model(triangulated_area,(T)1e4,poissons_ratio,(T).01);
            break;}
        case 28:{
            TRIANGULATED_AREA<T>& triangulated_area=solid_body_collection.deformable_body_collection.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(triangulated_area,(T)1e5,poissons_ratio,(T).01);
            break;}
        case 29: break;
        default:
            LOG::cerr<<"Missing implementation for test number "<<test_number<<std::endl;exit(1);}

    if(ether_drag) solid_body_collection.Add_Force(new ETHER_DRAG<TV>(particles,rigid_body_collection,true,true,ether_drag,0));

    if(solid_body_collection.deformable_body_collection.mpi_solids)
        solid_body_collection.deformable_body_collection.mpi_solids->Simple_Partition(solid_body_collection.deformable_body_collection,solid_body_collection.rigid_body_collection,particles.X,VECTOR<int,2>(2,1));
    solid_body_collection.Update_Simulated_Particles();

    solids_evolution->fully_implicit=!semi_implicit;

    if(scatter_plot) Init_Scatter_Plot();

    SOLIDS_EXAMPLE<TV>::Initialize_Bodies();
}
//#####################################################################
// Function Place_Triangle
//#####################################################################
void Place_Triangle(int tri,T s1,T s2,T a,T b,TV shift)
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    ROTATION<TV> rot(ROTATION<TV>::From_Angle(a)),rot2(ROTATION<TV>::From_Angle(b+a));
    TV com(0,sqrt(3)/6);
    for(int i=0;i<3;i++) particles.X(3*tri+i)=rot2.Rotate(TV(s1,s2)*rot.Inverse_Rotate(particles.X(3*tri+i)-com))+shift;
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
    return true;
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override
{
    if(test_number==20){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        TV velocity=velocity_time<5.0?attachment_velocity:TV();
        for(int j=0;j<n;j++){V(m*j)=-velocity;V(m-1+m*j)=velocity;}}
    if(test_number==24){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        T velocity=.2;
        T final_time=50;
        TV velocity_top=velocity_time<final_time?TV((T)0,velocity):TV();
        TV velocity_bot=velocity_time<final_time?TV((T)0,-velocity):TV();
        TV velocity_rig=velocity_time<final_time?TV(velocity,(T)0):TV();
        TV velocity_lef=velocity_time<final_time?TV(-velocity,(T)0):TV();
        for(int j=n/3;j<2*n/3+1;j++){V(m*j)=velocity_lef;V(m-1+m*j)=velocity_rig;}
        for(int i=m/3;i<2*m/3+1;i++){V(i)=velocity_bot;V(m*(n-1)+i)=velocity_top;}}
    if(test_number==25){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        T velocity=.2;
        T final_time=50;
        TV velocity_tr=velocity_time<final_time?TV(velocity,velocity):TV();
        TV velocity_tl=velocity_time<final_time?TV(-velocity,velocity):TV();
        TV velocity_br=velocity_time<final_time?TV(velocity,-velocity):TV();
        TV velocity_bl=velocity_time<final_time?TV(-velocity,-velocity):TV();
        V(0)=velocity_bl;
        V(m-1)=velocity_br;
        V(m*(n-1))=velocity_tl;
        V(m*n-1)=velocity_tr;}
    if(test_number==26){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        T velocity=.2;
        T final_time=50;
        TV velocity_top=velocity_time<final_time?TV((T)0,-velocity):TV();
        TV velocity_bot=velocity_time<final_time?TV((T)0,-velocity):TV();
        TV velocity_rig=velocity_time<final_time?TV(velocity,(T)0):TV();
        TV velocity_lef=velocity_time<final_time?TV(-velocity,(T)0):TV();
        for(int j=0;j<n;j++){V(m*j)=velocity_lef;V(m-1+m*j)=velocity_rig;}
        for(int i=2*m/5;i<3*m/5+1;i++){V(i)=velocity_bot;V(m*(n-1)+i)=velocity_top;}}
    if(test_number==27){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        T velocity=.1;
        T final_time=40;
        TV velocity_rig=velocity_time<final_time?TV(-velocity,(T)0):TV();
        TV velocity_lef=velocity_time<final_time?TV(velocity,(T)0):TV();
        for(int j=0;j<n;j++){V(m*j)=velocity_lef;V(m-1+m*j)=velocity_rig;}}
    if(test_number==30){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        T velocity=.01;
        T final_time=8000;
        TV velocity_rig=velocity_time<final_time?TV(velocity,(T)0):TV();
        TV velocity_lef=velocity_time<final_time?TV(-velocity,(T)0):TV();
        for(int j=0;j<n;j++){V(m*j)=velocity_lef;V(m-1+m*j)=velocity_rig;}}
    if(test_number==28){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        TV velocity=velocity_time<5.0?attachment_velocity:TV();
        for(int j=0;j<n;j++){V(m*j)=-velocity;V(m-1+m*j)=velocity;}}
    if(test_number==31){V(0)=TV(0,.2);V(1).y=0;V(2).y=0;}
    if(test_number==32){V(0)=V(1)=TV(0,0);V(2).x=0;}
    if(test_number==33){V(0)=V(1)=V(3)=V(5)=V(6)=TV();}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override
{
    if(test_number==20){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        for(int j=0;j<n;j++) V(m*j)=V(m-1+m*j)=TV();}
    if(test_number==24){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        for(int j=n/3;j<2*n/3+1;j++){V(m*j)=V(m-1+m*j)=TV();}
        for(int i=m/3;i<2*m/3+1;i++){V(i)=V(m*(n-1)+i)=TV();}}
    if(test_number==25){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        V(0)=TV();
        V(m-1)=TV();
        V(m*(n-1))=TV();
        V(m*n-1)=TV();}
    if(test_number==26){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        for(int j=0;j<n;j++){V(m*j)=TV();V(m-1+m*j)=TV();}
        for(int i=2*m/5;i<3*m/5+1;i++){V(i)=TV();V(m*(n-1)+i)=TV();}}
    if(test_number==27){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        for(int j=0;j<n;j++){V(m*j)=TV();V(m-1+m*j)=TV();}}
    if(test_number==30){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        for(int j=0;j<n;j++){V(m*j)=TV();V(m-1+m*j)=TV();}}
    if(test_number==28){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        for(int j=0;j<n;j++) V(m*j)=V(m-1+m*j)=TV();}
    if(test_number==31){V(0)=TV();V(1).y=0;V(2).y=0;}
    if(test_number==32){V(0)=V(1)=TV();V(2).x=0;}
    if(test_number==33){V(0)=V(1)=V(3)=V(5)=V(6)=TV();}
}
void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) override {
    /*if(test_number==270){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
    }*/
}
//#####################################################################
// Function Zero_Out_Enslaved_Position_Nodes
//#####################################################################
void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) override
{
    int m=mattress_grid.counts.x;
    for(int j=0;j<mattress_grid.counts.y;j++) X(m*j)=X(m-1+m*j)=TV(0,0);
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
        std::string output_file=LOG::sprintf("Standard_Tests/Test_%d/SV_%d",test_number,frame);
        svout.open(output_file.c_str());
    }
    if(test_number==33) Plot_Energy_Landscape();
    if(test_number==100 || test_number==31 || test_number==30) Plot_Contour_Landscape(frame);
    if(test_number==33 && frame==2)
    {
        DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
        particles.X(5)=TV(-.08,-.16);
        particles.X(3)=particles.X(5);
        particles.X(3).x*=-1;
    }
}
//#####################################################################
// Function Add_Constitutive_Model
//#####################################################################
void Add_Constitutive_Model(TRIANGULATED_AREA<T>& triangulated_area,T stiffness,T poissons_ratio,T damping, T cutoff=0.4, T efc=20)
{
    if(input_efc!=FLT_MAX) efc=input_efc;
    if(input_cutoff!=FLT_MAX) cutoff=input_cutoff;
    if(input_poissons_ratio!=-1) poissons_ratio=input_poissons_ratio;
    if(input_youngs_modulus!=0) stiffness=input_youngs_modulus;

    ISOTROPIC_CONSTITUTIVE_MODEL<T,2>* icm=0;
    if(use_extended_neohookean) icm=new NEO_HOOKEAN_EXTRAPOLATED<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_extended_neohookean2) icm=new NEO_HOOKEAN_EXTRAPOLATED2<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_extended_neohookean3) icm=new GENERAL_EXTRAPOLATED<T,2>(*new GEN_NEO_HOOKEAN_ENERGY<T>,stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_rc_ext) icm=new RC_EXTRAPOLATED<T,2>(*new GEN_NEO_HOOKEAN_ENERGY<T>,stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_rc2_ext) icm=new RC2_EXTRAPOLATED<T,2>(*new GEN_NEO_HOOKEAN_ENERGY<T>,stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_corotated) icm=new COROTATED<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else if(use_corotated_fixed) icm=new COROTATED_FIXED<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else{
        NEO_HOOKEAN<T,2>* nh=new NEO_HOOKEAN<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff);
        icm=nh;

        nh->use_constant_ife=use_constant_ife;}
    //std::cout << "Lambda= " << icm->constant_lambda << std::cout;
    //std::cout << "Mu   =" << icm->constant_mu << std::cout;
    solid_body_collection.Add_Force(Create_Finite_Volume(triangulated_area,icm));

    if(plot_energy_density) Plot_Energy_Density(icm,stiffness);
    if(primary_contour) Primary_Contour(*icm);
    if(scatter_plot || test_number==100 || test_number==31 || test_number==30){image_size*=10;Add_Primary_Contour_Segments(*icm);image_size/=10;}
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
    for(int i=0;i<TV::m;i++){
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
//#define YY(k) {for(int i=0;i<TV::m;i++){TV av=(h1.dd##k(i)+h0.dd##k(i))*df/2/e;TV dif=(h1.d##k.Transposed().Column(i)-h0.d##k.Transposed().Column(i))/e;const char*va=#k;sprintf(buff, "============ test ============ %s %8.5f %8.5f (%8.5f)\n", va, av.Magnitude(), dif.Magnitude(), (av-dif).Magnitude());LOG::cout<<buff;}}
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
    for(int i=0;i<20;i++){
        TV f;
        random.Fill_Uniform(f,0,2);
        T e=1e-4;
        TV df;
        random.Fill_Uniform(df,-e,e);
        f=f.Sorted().Reversed();
        if(random.Get_Uniform_Integer(0,1)==1) f(1)=-f(1);
        LOG::cout<<f<<std::endl;
        icm.Test(DIAGONAL_MATRIX<T,2>(f),1);
        Test_Model_Helper<RC_EXTRAPOLATED<T,2> >(&icm,f,df,e);
    }
    if(RC2_EXTRAPOLATED<T,2>* rc2=dynamic_cast<RC2_EXTRAPOLATED<T,2>*>(&icm))
        rc2->Test_Model();
    exit(0);
}
//#####################################################################
// Function Primary_Contour
//#####################################################################
void Primary_Contour(ISOTROPIC_CONSTITUTIVE_MODEL<T,2>& icm)
{
    ARRAY<VECTOR<T,3>,VECTOR<int,2> > img(1,image_size,1,image_size);
    for(int i=0;i<image_size;i++)
        for(int j=0;j<image_size;j++){
            T x=(2*i-image_size)*sigma_range/image_size+1e-5;
            T y=(2*j-image_size)*sigma_range/image_size;
            TV g=icm.P_From_Strain(DIAGONAL_MATRIX<T,2>(x,y),1).To_Vector();
            DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2> disd;
            icm.Isotropic_Stress_Derivative(DIAGONAL_MATRIX<T,2>(x,y),disd,1);
            SYMMETRIC_MATRIX<T,2> H(disd.x0000,disd.x1100,disd.x1111);
            DIAGONAL_MATRIX<T,2> ev;
            MATRIX<T,2> eigenvectors;
            H.Fast_Solve_Eigenproblem(ev,eigenvectors);
            TV evec=fabs(ev.x.x)>fabs(ev.x.y)?eigenvectors.Column(0):eigenvectors.Column(1);
            if(evec.Sum()<0) evec=-evec;
            T val=TV::Dot_Product(evec,g);
            img(VECTOR<int,2>(i,j))=VECTOR<T,3>(val<0,val>=0,0);}

    for(int i=0;i<image_size;i++){if(i%10>0 && i%10<5) img(VECTOR<int,2>(i,image_size/2))=img(VECTOR<int,2>(image_size/2,i))=VECTOR<T,3>(0,0,1);}

    PPM_FILE<T>::Write("primary_contour.ppm", img);
}
//#####################################################################
// Function Init_Scatter_Plot
//#####################################################################
void Init_Scatter_Plot()
{
    int s=0;
    for(int f=0;FINITE_VOLUME<TV,2>* force=solid_body_collection.deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,2>*>(f);f++){
        force->Update_Position_Based_State(0,true,false);
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
    for(int f=0,k=0;FINITE_VOLUME<TV,2>* force=solid_body_collection.deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,2>*>(f);f++){
        for(int i=0;i<force->Fe_hat.m;i++){
            if(!use_contrails) contrail(k).Remove_All();
            contrail(k++).Append(force->Fe_hat(i).To_Vector());}}
}
//#####################################################################
// Function Dump_Scatter_Plot
//#####################################################################
void Dump_Scatter_Plot(int frame)
{
    RANGE<TV> box(RANGE<TV>::Centered_Box()*image_size);
    char buff[1000];
    sprintf(buff, "svd-plot-%04d.eps", frame);
    EPS_FILE<T> eps(buff,box);
    eps.cur_format.point_radius=4*sigma_range/image_size;
    eps.fixed_bounding_box=true;
    eps.bounding_box=RANGE<TV>(RANGE<TV>::Centered_Box()*sigma_range);
    for(int i=0;i<contrail.m;i++){
        eps.cur_format.line_color=contrail_colors(i);
        for(int j=1;j<contrail(i).m;j++)
            eps.Draw_Object(contrail(i)(j-1),contrail(i)(j));}
    for(int i=0;i<contrail.m;i++){
        eps.cur_format.line_color=contrail_colors(i);
        eps.Draw_Object(contrail(i).Last());}
    eps.cur_format.line_color=VECTOR<T,3>();

    for(int i=0;i<contour_segments.m;i++)
        eps.Draw_Object(contour_segments(i).x,contour_segments(i).y);

    eps.Draw_Object(TV(0,-image_size),TV(0,image_size));
    eps.Draw_Object(TV(-image_size,0),TV(image_size,0));
}
//#####################################################################
// Function Contour_Crossing
//#####################################################################
T Contour_Crossing(const TV& g0,const TV& v0,const TV& g1,const TV& v1)
{
    T a=TV::Dot_Product(g0,v0);
    T b=TV::Dot_Product(g1,v1);
    if(!a) return 0;
    if(!b) return 1;
    if(TV::Dot_Product(v0,v1)<0) b=-b;
    if((a>0) == (b>0)) return -1;
    return a/(a-b);
}
//#####################################################################
// Function Add_Primary_Contour_Segments
//#####################################################################
void Add_Primary_Contour_Segments(ISOTROPIC_CONSTITUTIVE_MODEL<T,2>& icm)
{
    ARRAY<TV,VECTOR<int,2> > evec(1,image_size,1,image_size);
    ARRAY<TV,VECTOR<int,2> > grad(1,image_size,1,image_size);
    for(int i=0;i<image_size;i++)
        for(int j=0;j<image_size;j++){
            T x=(2*i-image_size)*sigma_range/image_size+1e-5;
            T y=(2*j-image_size)*sigma_range/image_size;
            TV g=icm.P_From_Strain(DIAGONAL_MATRIX<T,2>(x,y),1).To_Vector();
            DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2> disd;
            icm.Isotropic_Stress_Derivative(DIAGONAL_MATRIX<T,2>(x,y),disd,1);
            SYMMETRIC_MATRIX<T,2> H(disd.x0000,disd.x1100,disd.x1111);
            DIAGONAL_MATRIX<T,2> ev;
            MATRIX<T,2> eigenvectors;
            H.Fast_Solve_Eigenproblem(ev,eigenvectors);
            evec(VECTOR<int,2>(i,j))=fabs(ev.x.x)>fabs(ev.x.y)?eigenvectors.Column(0):eigenvectors.Column(1);
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
            if(X00.Product()>2) continue;
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
// Function Plot_Energy_Density
//#####################################################################
void Plot_Energy_Density(ISOTROPIC_CONSTITUTIVE_MODEL<T,2>* icm,T stiffness)
{
    TRIANGULATED_AREA<T>& ta=tests.Create_Mattress(GRID<TV>(TV_INT()+2+image_size+1,RANGE<TV>::Centered_Box()*sigma_range),true,RIGID_BODY_STATE<TV>());
    TRIANGULATED_SURFACE<T> ts(ta.mesh,*new DEFORMABLE_PARTICLES<VECTOR<T,3> >);
    ts.particles.Add_Elements(ta.particles.X.m);
    for(int i=0;i<ta.particles.X.m;i++){
        TV X=ta.particles.X(i);
        X.x+=1.1e-4;
        X.y+=2.3e-4;
        TV Z=X;
        if((fabs(Z.x)>fabs(Z.y)?Z.x:Z.y)<0) Z=-Z;
        VECTOR<T,3> Y(X.x,X.y,icm->P_From_Strain(DIAGONAL_MATRIX<T,2>(Z),0).x.x/stiffness);
        ts.particles.X(i)=Y;}
    FILE_UTILITIES::Write_To_File(this->stream_type,"surface.tri",ts);
}
//#####################################################################
// Function Plot_Energy_Landscape
//#####################################################################
void Plot_Energy_Landscape()
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
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
            solid_body_collection.Update_Position_Based_State(0,false,false);
            T ke,pe;
            solid_body_collection.Compute_Energy(0,ke,pe);
            solid_body_collection.Add_Velocity_Independent_Forces(F,TW,0);
            E.Append(ke+pe);
            X.Append(particles.X(5));
            N.Append(F(5).Normalized());
        }

    INTERPOLATED_COLOR_MAP<T> mp;
    mp.Initialize_Colors(E.Min(),E.Max(),false,true,false);

    for(int i=0;i<X.m;i++){
        Add_Debug_Particle(X(i),mp(E(i)));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,N(i));}
    particles.X=TX;
}
//#####################################################################
// Function Plot_Contour_Landscape
//#####################################################################
void Plot_Contour_Landscape(int frame)
{
    char buff[1000];
    sprintf(buff, "%s/data", output_directory.c_str());
    FILE_UTILITIES::Create_Directory(buff);
    sprintf(buff, "%s/data/%03d.txt", output_directory.c_str(), frame);
    std::ofstream out(buff);

    if(test_number==31) out<<"tristretch"<<std::endl;

    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    FINITE_VOLUME<TV,2>& fv=solid_body_collection.deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,2>&>();
    ISOTROPIC_CONSTITUTIVE_MODEL<T,2>* icm=fv.isotropic_model;
    bool is_neo=dynamic_cast<NEO_HOOKEAN<T,2>*>(icm) || dynamic_cast<RC2_EXTRAPOLATED<T,2>*>(icm);
    T min=-image_size/2,max=image_size/2;
    for(int i=min;i<max;i++)
        for(int j=min;j<max;j++){
//            if(is_neo && (i<0 || j<0)) continue;
            TV X(2*sigma_range*(i+.5)/image_size+1.1e-5,2*sigma_range*(j+.5)/image_size+1.2e-5);
            DIAGONAL_MATRIX<T,2> F(X),P=icm->P_From_Strain(F,0)*plot_scale;
            out<<"p "<<X<<" "<<-P.To_Vector().Normalized()<<std::endl;}

    for(int i=0;i<fv.Fe_hat.m;i++){
        contrail(i).Append(fv.Fe_hat(i).To_Vector());
        out<<"c "<<contrail(i)<<std::endl;}

    for(int i=0;i<contour_segments.m;i++){
        if(is_neo && (contour_segments(i).x.Min()<.2 || contour_segments(i).y.Min()<.2 || contour_segments(i).y.Max()<.5)) continue;
        out<<"u "<<contour_segments(i)<<std::endl;}

    out<<"t "<<particles.X.Subset(fv.strain_measure.mesh.elements.Flattened())<<std::endl;
}
//#####################################################################
// Function Energy_Profile_Plot
//#####################################################################
void Energy_Profile_Plot(int frame)
{
    if(!dual_directory.length()) dual_directory=output_directory+"/dual_output";

    if(this->write_first_frame && frame==this->first_frame) FILE_UTILITIES::Write_To_Text_File(dual_directory+"/common/first_frame",frame,"\n");

    ISOTROPIC_CONSTITUTIVE_MODEL<T,2>* icm=solid_body_collection.deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,2>&>().isotropic_model;
    if(!energy_mesh){
        energy_particles.Store_Velocity(true);
        TRIANGULATED_AREA<T> ta(*new TRIANGLE_MESH,*new GEOMETRY_PARTICLES<TV>);
        ta.Initialize_Square_Mesh_And_Particles(GRID<TV>(TV_INT()+2+image_size+1,RANGE<TV>::Centered_Box()*sigma_range,false));

        energy_mesh=new TRIANGULATED_SURFACE<T>(ta.mesh,*new DEFORMABLE_PARTICLES<VECTOR<T,3> >);
        energy_mesh->particles.Add_Elements(ta.particles.X.m);
        for(int i=0;i<ta.particles.X.m;i++){
            TV X=ta.particles.X(i);
            X.x+=1.1e-4;
            X.y+=2.3e-4;
            TV Z=X;
            if((fabs(Z.x)>fabs(Z.y)?Z.x:Z.y)<0) Z=-Z;
            VECTOR<T,3> Y(X.x,X.y,icm->Energy_Density(DIAGONAL_MATRIX<T,2>(Z),0));
            energy_mesh->particles.X(i)=Y;}

        PROJECTED_ARRAY<ARRAY_VIEW<VECTOR<T,3>,int>,FIELD_PROJECTOR<VECTOR<T,3>,T,&VECTOR<T,3>::z> > pa=energy_mesh->particles.X.template Project<T,&VECTOR<T,3>::z>();
        pa-=energy_profile_plot_min=pa.Min();
        energy_profile_plot_range=pa.Max();
        pa*=plot_scale;}

    FILE_UTILITIES::Create_Directory(dual_directory);
    std::string f=LOG::sprintf("%d",frame);
    FILE_UTILITIES::Create_Directory(dual_directory+"/"+f);
    FILE_UTILITIES::Create_Directory(dual_directory+"/common");

    LOG::cout<<this->debug_particles.debug_particles.X.m<<std::endl;
    energy_particles.Add_Elements(this->debug_particles.debug_particles.X.m);
    LOG::cout<<energy_particles.X.m<<std::endl;
    for(int i=0;i<energy_particles.X.m;i++){
        TV X=this->debug_particles.debug_particles.X(i);
        DIAGONAL_MATRIX<T,2> F(X),P=icm->P_From_Strain(F,0)*plot_scale;
        energy_particles.X(i)=X.Append((icm->Energy_Density(F,0)-energy_profile_plot_min)*plot_scale+1e-2);
        energy_particles.V(i)=-P.To_Vector().Append(P.To_Vector().Magnitude_Squared()).Normalized();}

    FILE_UTILITIES::Create_Directory(LOG::sprintf("%s/%i",dual_directory.c_str(),frame));
    FILE_UTILITIES::Write_To_File(this->stream_type,LOG::sprintf("%s/%i/debug_particles",dual_directory.c_str(),frame),energy_particles);
    energy_particles.Delete_All_Elements();

    FILE_UTILITIES::Write_To_File(this->stream_type,dual_directory+"/"+FILE_UTILITIES::Number_To_String(frame)+"/deformable_object_particles",energy_mesh->particles);
    if(frame==1 || (this->restart && frame==this->first_frame)){
        FILE_UTILITIES::Create_Directory(LOG::sprintf("%s/%i",dual_directory.c_str(),0));
        FILE_UTILITIES::Write_To_File(this->stream_type,LOG::sprintf("%s/%i/debug_particles",dual_directory.c_str(),0),energy_particles);
        FILE_UTILITIES::Write_To_File(this->stream_type,dual_directory+"/0/deformable_object_particles",energy_mesh->particles);
        std::string f="common";
        std::ostream* output_raw=FILE_UTILITIES::Safe_Open_Output(dual_directory+"/"+f+"/deformable_object_structures");
        TYPED_OSTREAM output(*output_raw,this->stream_type);
        Write_Binary(output,1);
        energy_mesh->Write_Structure(output);
        delete output_raw;}

    FILE_UTILITIES::Write_To_Text_File(dual_directory+"/common/last_frame",frame,"\n");
}
//#####################################################################
};
}

#endif

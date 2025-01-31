//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
//    1. Sphere free falling to the ground
//    2. Torus free falling to the ground
//    3. Maggot free falling to the ground
//    4. Large armadillo caving in on itself
//    5. Deformable ball falling on a rigid ball
//    6. Smash test - small boxes, large sphere
//    7. Plane test - ball
//    8. Falling mattress
//    9. Deformable ball falling on a deformable ball
//   10. Tori free falling into bowl
//   11. Increasing gravity (individual)
//   12. Mattress with ether drag
//   13. Falling mattress hitting level set
//   16. Smash test - large boxes, small mattress
//   17. Matress, no gravity, random start
//   18. Matress, no gravity, point start
//   23. Big lateral stretch
//   24. Big 6 sides stretch
//   25. Big 8 corners stretch
//   26. Big stretch/bend
//   27. Force inversion
//   28. Taffy test
//   29. Armadillo collapsing and rebounding
//   30. Projectile hitting a wall
//   31. Impact with sphere
//   32  Twisting chain
//   33. Through gears
//   34. Roll-over
//   35. Dancing jello (first attempt)
//   36. Spinning jello
//   37. Two jellos falling and colliding (first high res video)
//   38. Bunch of jellos falling
//   39. One jello falling on a stationary one
//   40. Two jellos running on each other
//   41. Bunch of jellos rolling towards camera
//   42. Bunch of jellos rolling towards camera (Reloaded)
//   43. Through smooth gears
//   44. 2 jellos collision
//   47. Fish past a magnet?
//   48. Compress with sphere
//   49. See-saw?
//   50. Fish through a torus
//   51. Fish through a tube
//   52  Jello's falling one by one on each other
//   53. Stretch tet
//   54. Stretch tet (II)
//   55. Constrained tet
//   56. Size comparison - several shapes
//   57. Two-direction stretch
//   58. Various objects through gears
//   59. Random armadillo
//   60. Random armadillo
//   61. Elbow cylinder
//   62. Twisting mattress
//   63. Twisting and stretching mattress
//   64. Stretching mattress lucky vanilla newton failure
//   65. Stretching mattress with shock
//   66. Force a deformable torus through a rigid torus
//   67. Two mattresses collide
//   68. Mattress stack
//   77. Squeeze in a box
//   80. Armadillo
//  100. Primary contour field
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <Core/Arrays/PROJECTED_ARRAY.h>
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/pow.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_UTILITIES.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/SMOOTH_LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Forces_And_Torques/MOVE_RIGID_BODY_DIFF.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Deformables/Bindings/RIGID_BODY_BINDING.h>
#include <Deformables/Bindings/SOFT_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Collisions_And_Interactions/LEVELSET_VOLUME_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_PENALTY.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/ELASTIC_ETHER_DRAG.h>
#include <Deformables/Forces/FINITE_VOLUME.h>
#include <Deformables/Forces/IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION.h>
#include <Deformables/Forces/RALEIGH_DAMPING_FORCE.h>
#include <Deformables/Forces/SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION.h>
#include <Deformables/Forces/TRIANGLE_BENDING_SPRINGS.h>
#include <Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
#include <Solids/Forces_And_Torques/ETHER_DRAG.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
#include <Solids/Forces_And_Torques/RIGID_DEFORMABLE_PENALTY_WITH_FRICTION.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_EVOLUTION.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_MINIMIZATION_OBJECTIVE.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_MINIMIZATION_SYSTEM.h>
#include <Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include "STANDARD_TESTS_BASE.h"
#include <fstream>
#ifdef USE_OPENMP
#include <omp.h>
#endif
namespace PhysBAM{

template<class TV> class STANDARD_TESTS;
template<class T_input>
class STANDARD_TESTS<VECTOR<T_input,3> >:public STANDARD_TESTS_BASE<VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    typedef STANDARD_TESTS_BASE<TV> BASE;
    using BASE::solids_parameters;using BASE::viewer_dir;
    using BASE::last_frame;using BASE::frame_rate;
    using BASE::solid_body_collection;using BASE::stream_type;
    using BASE::solids_evolution;using BASE::test_number;
    using BASE::data_directory;using BASE::m;using BASE::s;using BASE::kg;
    using BASE::unit_p;using BASE::unit_J;using BASE::unit_rho;
    using BASE::backward_euler_evolution;using BASE::Add_Collision_Object;
    using BASE::rd_penalty_stiffness;using BASE::rd_penalty_friction;
    using BASE::user_last_frame;

    SOLIDS_STANDARD_TESTS<TV> tests;

    GRID<TV> mattress_grid,mattress_grid2,mattress_grid3,mattress_grid1;
    T attachment_velocity;
    bool with_bunny,with_hand,with_big_arm,gears_of_pain;
    bool override_collisions,override_no_collisions;
    ARRAY<int> kinematic_ids;
    ARRAY<INTERPOLATION_CURVE<T,FRAME<TV> > > curves;
    INTERPOLATION_CURVE<T,T> scalar_curve;
    bool print_matrix;
    int resolution;
    int fishes,jello_size,number_of_jellos;
    T stiffness_multiplier;
    T damping_multiplier;
    T degrees_wedge,degrees_incline;
    T rebound_time,rebound_stiffness,rebound_drop;
    bool forces_are_removed,self_collision_flipped,sloped_floor;
    ARRAY<int> externally_forced;
    ARRAY<int> constrained_particles;
    ARRAY<TV> constrained_velocities;
    ARRAY<TV> jello_centers;
    T stretch,plateau,repulsion_thickness;
    T hole;
    bool nobind;
    ARRAY<TV> fish_V;
    T input_cutoff;
    T input_efc;
    T input_poissons_ratio,input_youngs_modulus;
    T input_friction,stretch_cutoff;
    T hand_scale;
    T ether_drag;
    ARRAY<ARRAY<VECTOR<T,2> > > contrail;
    int image_size;
    ARRAY<VECTOR<VECTOR<T,2>,2> > contour_segments;
    ARRAY<int> stuck_particles;
    bool pin_corners;
    int tori_stack_width;
    int tori_stack_height;
    int rand_seed;
    bool use_rand_seed;
    RANDOM_NUMBERS<T> rand;
    bool use_residuals;
    bool project_nullspace;
    bool use_penalty_collisions;
    bool use_constraint_collisions;
    T penalty_collisions_stiffness,penalty_collisions_separation,penalty_collisions_length;
    bool enforce_definiteness;
    T density;
    bool use_penalty_self_collisions;
    bool use_distance_based_self_collisions;
    T rod_length;
    T rod_radius;
    T attachment_length;
    ARRAY<TV> initial_positions;
    T save_dt;
    bool self_collide_surface_only;
    T collision_height;
    T collision_speed;

    STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :BASE(stream_type_input,parse_args),tests(stream_type_input,data_directory,solid_body_collection),
        with_bunny(false),with_hand(false),with_big_arm(false),gears_of_pain(false),override_collisions(false),override_no_collisions(false),
        print_matrix(false),resolution(0),fishes(0),jello_size(20),number_of_jellos(12),stiffness_multiplier(1),damping_multiplier(1),
        degrees_wedge(2),degrees_incline(5),rebound_time((T).2),rebound_stiffness(5),rebound_drop((T)1.5),
        forces_are_removed(false),self_collision_flipped(false),sloped_floor(false),stretch(1),plateau(0),
        repulsion_thickness((T)1e-4),hole((T).5),nobind(false),input_cutoff((T).4),input_efc(20),input_poissons_ratio(-1),
        input_youngs_modulus(0),input_friction(.3),stretch_cutoff(300),hand_scale((T).8),
        ether_drag(0),image_size(500),pin_corners(false),tori_stack_width(5),tori_stack_height(5),rand_seed(1234),
        use_rand_seed(false),use_residuals(false),project_nullspace(false),
        use_penalty_collisions(false),use_constraint_collisions(true),
        penalty_collisions_stiffness((T)1e4),penalty_collisions_separation((T)1e-4),
        penalty_collisions_length(1),enforce_definiteness(false),density(pow<TV::m>(10)),
        use_penalty_self_collisions(true),use_distance_based_self_collisions(false),rod_length(4),rod_radius(.3),attachment_length(.6),save_dt(1),self_collide_surface_only(false),
        collision_height(5),collision_speed(4)
    {
        parse_args.Add("-with_bunny",&with_bunny,"use bunny");
        parse_args.Add("-with_hand",&with_hand,"use hand");
        parse_args.Add("-with_big_arm",&with_big_arm,"use big arm");
        parse_args.Add("-resolution",&resolution,"resolution","resolution used by multiple tests to change the parameters of the test");
        parse_args.Add("-stiffen",&stiffness_multiplier,"multiplier","stiffness multiplier for various tests");
        parse_args.Add("-dampen",&damping_multiplier,"multiplier","damping multiplier for various tests");
        parse_args.Add("-residuals",&use_residuals,"print residuals during timestepping");
        parse_args.Add("-print_matrix",&print_matrix,"print Krylov matrix");
        parse_args.Add("-project_nullspace",&project_nullspace,"project out nullspace");
        parse_args.Add("-seed",&rand_seed,&use_rand_seed,"seed","random seed to use");
        parse_args.Add("-collisions",&override_collisions,"Does not yet work in all sims, see code for details");
        parse_args.Add("-no_collisions",&override_no_collisions,"Does not yet work in all sims, see code for details");
        parse_args.Add("-stretch",&stretch,"stretch","stretch");
        parse_args.Add("-hole",&hole,"hole","hole size");
        parse_args.Add("-rebound_time",&rebound_time,"time","number of seconds to rebound in test 29");
        parse_args.Add("-rebound_stiffness",&rebound_stiffness,"stiffness","log10 of youngs modulus of final stiffness");
        parse_args.Add("-rebound_drop",&rebound_drop,"drop","log10 of youngs modulus of dropoff of final stiffness");
        parse_args.Add("-nobind",&nobind,"do not bind particles to rings");
        parse_args.Add("-cutoff",&input_cutoff,"cutoff","cutoff");
        parse_args.Add("-repulsion_thickness",&repulsion_thickness,"thickness","repulsion thickness");
        parse_args.Add("-efc",&input_efc,"efc","efc");
        parse_args.Add("-poissons_ratio",&input_poissons_ratio,"ratio","poissons_ratio");
        parse_args.Add("-youngs_modulus",&input_youngs_modulus,"stiffness","youngs modulus, only for test 41 so far");
        parse_args.Add("-jello_size",&jello_size,"size","resolution of each jello cube");
        parse_args.Add("-number_of_jellos",&number_of_jellos,"number","number of falling jello cubes in test 41");
        parse_args.Add("-degrees_incline",&degrees_incline,"angle","degrees of incline");
        parse_args.Add("-degrees_wedge",&degrees_wedge,"angle","degrees of side incline");
        parse_args.Add("-friction",&input_friction,"friction","amount of friction");
        parse_args.Add("-gears_of_pain",&gears_of_pain,"use gears of pain");
        parse_args.Add("-sloped_floor",&sloped_floor,"use sloped floor");
        parse_args.Add("-hand_scale",&hand_scale,"scale","hand scale on test 58");
        parse_args.Add("-stretch_cutoff",&stretch_cutoff,"cutoff","Stretch cutoff on test 58");
        parse_args.Add("-ether_drag",&ether_drag,"drag","Ether drag");
        parse_args.Add("-image_size",&image_size,"size","image size for plots");
        parse_args.Add("-pin_corners",&pin_corners,"pin corners");
        parse_args.Add("-use_penalty",&use_penalty_collisions,"use penalty collisions");
        parse_args.Add_Not("-no_constraints",&use_constraint_collisions,"disable constrained optimization for collisions");
        parse_args.Add("-penalty_stiffness",&penalty_collisions_stiffness,"tol","penalty collisions stiffness");
        parse_args.Add("-penalty_separation",&penalty_collisions_separation,"tol","penalty collisions separation");
        parse_args.Add("-penalty_length",&penalty_collisions_length,"tol","penalty collisions length scale");
        parse_args.Add("-enf_def",&enforce_definiteness,"enforce definiteness in system");
        parse_args.Add_Not("-no_self",&use_penalty_self_collisions,"disable penalty self collisions");
        parse_args.Add("-old_self",&use_distance_based_self_collisions,"use distance based penalty self collisions");
        parse_args.Add("-no_self_interior",&self_collide_surface_only,"do not process penalty self collisions against interior particles");
        parse_args.Add("-collision_height",&collision_height,"height","height of collision body in test 68");
        parse_args.Add("-collision_speed",&collision_speed,"speed","speed of collision body in test 68");
        parse_args.Parse();

        LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
        if(!this->user_output_directory)
            viewer_dir.output_directory=LOG::sprintf("Test_%d",test_number);
        override_no_collisions=override_no_collisions&&!override_collisions;
        if(use_rand_seed) rand.Set_Seed(rand_seed);
        solids_parameters.implicit_solve_parameters.project_nullspace_frequency=project_nullspace;
        hole*=m;
        rebound_time*=s;
        repulsion_thickness*=m;
        penalty_collisions_length*=m;
        penalty_collisions_separation*=m;
        input_youngs_modulus*=unit_p;
        ether_drag/=s;
        penalty_collisions_stiffness*=unit_J;
        density*=unit_rho;
        if(!this->fixed_dt && !this->max_dt && !this->min_dt) this->fixed_dt=1;

        switch(test_number){
            case 17: case 18: case 24: case 25: case 27: case 11: case 23: case 57: case 77: case 80: case 8: case 12: case 13: case 67: case 68:
                if(!resolution) resolution=10;
                mattress_grid=GRID<TV>(TV_INT(resolution+1,resolution+1,resolution+1),RANGE<TV>(TV((T)-1,(T)-1,(T)-1),TV((T)1,(T)1,(T)1))*m);
                break;
            case 34:
                mattress_grid=GRID<TV>(TV_INT(13,13,13),RANGE<TV>(TV((T)-2,(T)-2,(T)-2),TV((T)2,(T)2,(T)2))*m);
                break;
            case 26:
                mattress_grid=GRID<TV>(TV_INT(40,5,5),RANGE<TV>(TV((T)-4,(T)-.5,(T)-.5),TV((T)4,(T).5,(T).5))*m);
                break;
            case 28:
                mattress_grid=GRID<TV>(TV_INT(80,10,10),RANGE<TV>(TV((T)-8,(T)-.5,(T)-.5),TV((T)8,(T).5,(T).5))*m);
                break;
            case 16:
                mattress_grid=GRID<TV>(TV_INT(11,6,11),RANGE<TV>(TV((T)-1,(T)-.5,(T)-1),TV((T)1,(T).5,(T)1))*m);
                break;
            case 35: case 36: case 41:
                mattress_grid=GRID<TV>(TV_INT(jello_size,jello_size,jello_size),RANGE<TV>(TV((T)-0.01,(T)-0.01,(T)-0.01),TV((T)0.01,(T)0.01,(T)0.01))*m);
                mattress_grid1=GRID<TV>(TV_INT(jello_size,jello_size,jello_size),RANGE<TV>(TV((T)-0.01,(T)-0.01,(T)-0.01),TV((T)0.01,(T)0.01,(T)0.01))*m);
                mattress_grid2=GRID<TV>(TV_INT(jello_size,jello_size,jello_size),RANGE<TV>(TV((T)-0.016,(T)-0.016,(T)-0.016),TV((T)0.016,(T)0.016,(T)0.016))*m);
                mattress_grid3=GRID<TV>(TV_INT(jello_size,jello_size,jello_size),RANGE<TV>(TV((T)-0.0125,(T)-0.0125,(T)-0.0125),TV((T)0.0125,(T)0.0125,(T)0.0125))*m);
                break;
            case 37: case 39: case 40: case 38: case 44:
                mattress_grid=GRID<TV>(TV_INT(10,10,10),RANGE<TV>(TV((T)-0.01,(T)-0.01,(T)-0.01),TV((T)0.01,(T)0.01,(T)0.01))*m);
                break;
            case 42: case 52:
                mattress_grid=GRID<TV>(TV_INT(20,20,20),RANGE<TV>(TV((T)-0.01,(T)-0.01,(T)-0.01),TV((T)0.01,(T)0.01,(T)0.01))*m);
                break;
            case 62:
                if(!resolution) resolution=10;
                mattress_grid=GRID<TV>(TV_INT(5*resolution+1,resolution+1,resolution+1),RANGE<TV>(TV((T)-5,(T)-1,(T)-1),TV((T)5,(T)1,(T)1))*m);
                break;
            case 63: case 64: case 65:
                if(!resolution) resolution=10;
                mattress_grid=GRID<TV>(TV_INT(5*resolution+1,resolution+1,resolution+1),RANGE<TV>(TV((T)-1,(T)-1,(T)-1),TV((T)1,(T)1,(T)1))*m);
                break;
            case 740:
                mattress_grid=GRID<TV>(TV_INT(resolution+1,resolution+1,resolution+1),RANGE<TV>(TV((T)-.1,(T)-.1,(T)-.1),TV((T).1,(T).1,(T).1))*m);
                break;
            default:{
                if(!resolution) resolution=10;

                mattress_grid=GRID<TV>(TV_INT(2*resolution,resolution,2*resolution),RANGE<TV>(TV((T)-1,(T)-.5,(T)-1),TV((T)1,(T).5,(T)1))*m);
            }
        }

        switch(test_number){
            case 1:
            case 3:
                // case 4:
            case 7:
            case 8:
            case 12:
            case 13:
            case 80:
            case 11:
            case 16:
            case 17:
            case 18:
            case 56:
            case 67:
            case 68:
            case 69:
            case 77:
                solids_parameters.cfl=(T)5;
                /* solids_parameters.implicit_solve_parameters.cg_iterations=100000; */
                break;
            case 701:case 720:case 721:case 722:
            case 730:
            case 740:
            case 750:case 751:
            case 130:{
                // Test rigid-deformable penalty force with friction.
                // ./be_evolution 701 -no_collisions_in_solve -rd_stiffness 1e2
                break;}
            case 140:{
                // Cloth, rigid and deformable tori
                break;}
            case 2:
                if(!user_last_frame) last_frame=600;
                break;
            case 10:
                attachment_velocity=6;
                if(!user_last_frame) last_frame=330;
                break;
            case 5:
            case 6:
            case 9:
                if(!user_last_frame) last_frame=72;
                solids_parameters.cfl=(T)5.9;
                /* solids_parameters.implicit_solve_parameters.cg_iterations=100000; */
                solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
                break;
            case 23:
                attachment_velocity=1.0;
                solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
                if(!user_last_frame) last_frame=2000;
                break;
            case 24:
            case 25:
            case 26:
            case 27: case 53: case 54: case 55: case 57: case 100: case 48:
                attachment_velocity=0.2;
                solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
                if(!user_last_frame) last_frame=2000;
                break;
            case 61:
                solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
                if(!user_last_frame) last_frame=360;
                break;
            case 62: case 63:
                solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
                if(!user_last_frame) last_frame=2000;
                break;
            case 64: case 65:
                solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
                if(!user_last_frame) last_frame=120;
                break;
            case 28:
                attachment_velocity=0.4;
                solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                //solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
                if(!user_last_frame) last_frame=1000;
            case 29: case 4:
                solids_parameters.cfl=(T)5;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=true;
                solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness=1e-5;

                solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=override_collisions;
                solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=override_collisions;
                solids_parameters.triangle_collision_parameters.perform_self_collision=override_collisions;//This gets turned off later then back on
                //std::cout << "rame collisions are " << override_collisions << std::endl;
                //}
                if(!user_last_frame) last_frame=84;
                break;
            case 30:
                solids_parameters.cfl=(T)5;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=true;
                break;
            case 31:
                solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness=1e-4;
                solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                /* solids_parameters.triangle_collision_parameters.perform_self_collision=true; */
                if(!user_last_frame) last_frame=500;
                break;
            case 59:
            case 60:
                solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness=1e-4;
                solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
                if(!user_last_frame) last_frame=300;
                break;
            case 32:
                // solids_parameters.cfl=(T)5;
                solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                solids_parameters.triangle_collision_parameters.perform_self_collision=true;
                solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=true;
                solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;
                if(!user_last_frame) last_frame=4000;
                break;
            case 33:
                solids_parameters.cfl=(T)10;
                solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
                solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness=3e-3;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                solids_parameters.triangle_collision_parameters.perform_self_collision=true;
                solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=true;
                solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=true;
                solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;
                if(!user_last_frame) last_frame=600;
                break;
            case 43:
            case 58:
                solids_parameters.cfl=(T)10;
                solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness=repulsion_thickness;
                solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=true;
                solids_parameters.triangle_collision_parameters.perform_self_collision=override_collisions;
                solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=override_collisions;
                solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=override_collisions;
                if(!user_last_frame) last_frame=300;
                break;
            case 34:
                solids_parameters.cfl=(T)5;
                solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                if(!user_last_frame) last_frame=400;
                break;
            case 35: case 36: //case 41:
                solids_parameters.cfl=(T)5;
                solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                solids_parameters.triangle_collision_parameters.perform_self_collision=true;
                // solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=true;
                // solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;
                if(!user_last_frame) last_frame=1000;
                break;
            case 37:
            case 38:
            case 39:
            case 40:
            case 42:
            case 44:
                solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness=1e-4;
                solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                solids_parameters.triangle_collision_parameters.perform_self_collision=true;
                if (override_no_collisions) solids_parameters.triangle_collision_parameters.perform_self_collision=false;
                if(!user_last_frame) last_frame=250;
                break;
            case 41:
                solids_parameters.cfl=(T)5;
                solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness=2e-4;
                solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                solids_parameters.triangle_collision_parameters.perform_self_collision=true;
                solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=true;
                solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;

                if (override_no_collisions){
                    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
                    solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=false;
                    solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=false;
                }
                if(!user_last_frame) last_frame=40;
                break;
            case 52:
                solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness=2e-4;
                solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                solids_parameters.triangle_collision_parameters.perform_self_collision=true;
                solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=true;
                solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;

                if (override_no_collisions){
                    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
                    solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=false;
                    solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=false;
                }
                if(!user_last_frame) last_frame=250;
                break;
            case 47:
                if(!user_last_frame) last_frame=480;
                solids_parameters.triangle_collision_parameters.perform_self_collision=true;
                //solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=true;
                //solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;
                solids_parameters.cfl=(T)5;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                break;
            case 49:
                if(!user_last_frame) last_frame=480;
                solids_parameters.triangle_collision_parameters.perform_self_collision=true;
                //solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=true;
                //solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;
                solids_parameters.cfl=(T)5;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                break;
            case 50:
                solids_parameters.triangle_collision_parameters.perform_self_collision=true;
                //solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=true;
                //solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;
                solids_parameters.cfl=(T)5;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                break;
            case 51:
                solids_parameters.triangle_collision_parameters.perform_self_collision=true;
                //solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=true;
                //solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;
                solids_parameters.cfl=(T)5;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                break;
            case 66:
                // solids_parameters.cfl=(T)5;
                solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
                solids_parameters.implicit_solve_parameters.cg_iterations=100000;
                if(!user_last_frame) last_frame=120;
                break;
            case 200:
            case 201:
            case 202:
                // solids_parameters.cfl=(T)5;
                solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
                solids_parameters.implicit_solve_parameters.cg_iterations=10;
                if(!user_last_frame) last_frame=120;
                break;
            default:
                LOG::cerr<<"Parsing: Unrecognized test number "<<test_number<<std::endl;exit(1);}

        //solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=override_collisions;
        //solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=override_collisions;
        //solids_parameters.triangle_collision_parameters.perform_self_collision=override_collisions;


        if(use_penalty_collisions || use_constraint_collisions){
            solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=false;
            solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=false;
            solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;}

        solid_body_collection.Print_Residuals(use_residuals);
    }

    virtual ~STANDARD_TESTS()
    {}

    void Postprocess_Frame(const int frame) override
    {
        DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
        DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
        for(int i=0;i<constrained_particles.m;i++) Add_Debug_Particle(particles.X(constrained_particles(i)),TV(1,0,0));
        for(int i=0;i<externally_forced.m;i++) Add_Debug_Particle(particles.X(externally_forced(i)),TV(0,1,0));
        if(test_number==722)
            for(int i=0;i<particles.X.m;i++)
                Add_Debug_Particle(particles.X(i),TV(1,0,0));
        if(test_number==740){
            T angle=.1;
            T vel=.1;
            ROTATION<TV> Q(angle,TV(0,0,1));
            TV n,t,t1;
            Q.Get_Rotated_Frame(t,n,t1);
            TV c(.5,.5,.5);
            TV g=m/(s*s)*TV(0,-1.8,0);
            T frame_dt=1/frame_rate;
            T time=frame*frame_dt;
            T mu=rd_penalty_friction;
            T acc=g.y*(cos(angle)*mu-sin(angle));
            if(vel+time*acc>0) // sliding
            {
                T dist=(T).5*acc*sqr(time)+time*vel;
                Add_Debug_Particle(c-t*dist,VECTOR<T,3>(0,1,0));
            }
            else // stopped
            {
                T stop_time=-vel/acc;
                T dist=(T).5*acc*sqr(stop_time)+stop_time*vel;
                Add_Debug_Particle(c-t*dist,VECTOR<T,3>(1,1,0));
            }
            Add_Debug_Object(VECTOR<TV,2>(c-t,c+t),VECTOR<T,3>(1,0,0));}
    }
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    bool automatically_add_to_collision_structures=true;
    // deformable bodies
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    switch(test_number){
        case 1: case 7:{
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)5,0)*m)),true,true,density,m);
            tests.Add_Ground(0,1.99*m);
            break;}
        case 2:{
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/adaptive_torus_float.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0)*m)),true,true,density,m);
            tests.Add_Ground();
            break;}

        case 3:{
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/maggot_8K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0)*m)),true,true,density,m);
            tests.Add_Ground();
            break;}
        case 4: case 29:{
            //tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_110K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T)0.373,(T)0)*m)),true,true,density,.005*m);

            if (with_hand)
                tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/hand_30k.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T).4,(T)0)*m,ROTATION<TV>(-T(pi/2),TV(1,0,0)))),true,true,density,m);
            else if(with_bunny)
                tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/bunny.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T).4,(T)0)*m)),true,true,density,1.0*m);
            else if(with_big_arm)
                tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_380K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T).4,(T)0)*m,ROTATION<TV>(T(pi),TV(0,1,0)))),true,true,density,.005*m);
            else
                tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_110K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T).4,(T)0)*m,ROTATION<TV>(T(pi),TV(0,1,0)))),true,true,density,.005*m);
            //            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/bunny.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)2.3,0)*m)),true,true,density,5.0*m);
//            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0)*m)),true,true,density,m);
           // tests.Create_Mattress(GRID<TV>(TV_INT(13,13,13),RANGE<TV>(TV(-1,1,-1),TV(1,3,1))*m),true,0);

            tests.Add_Ground();

            break;}
        case 5:{
            RIGID_BODY<TV>& tmp_sphere=tests.Add_Rigid_Body("sphere",m,(T).5);
            //RIGID_BODY<TV>& tmp_sphere=tests.Add_Analytic_Box(TV(1,1,1)*m);
            tmp_sphere.Frame().t=TV(0,(T).25,0)*m;
            tmp_sphere.is_static=true;
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_12k.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(.3,(T)6,0)*m)),true,true,density,.5*m);
            tests.Add_Ground();
            break;}
        case 6:{
            //RIGID_BODY<TV>& tmp_sphere=tests.Add_Rigid_Body("sphere",m,(T).5);
            RIGID_BODY<TV>& bottom_box=tests.Add_Analytic_Box(TV(1,1,1)*m);
            RIGID_BODY<TV>& top_box=tests.Add_Analytic_Box(TV(1,1,1)*m);
            bottom_box.Frame().t=TV(0,(T)0,0)*m;
            top_box.Frame().t=TV(0,(T)2,0)*m;
            bottom_box.is_static=true;
            top_box.is_static=false;
            kinematic_ids.Append(top_box.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(top_box.particle_index)=true;
            curves.Resize(1);
            curves(0).Add_Control_Point(0,FRAME<TV>(TV(0,(T)2,0)));
            curves(0).Add_Control_Point(1,FRAME<TV>(TV(0,(T)1,0)));
            curves(0).Add_Control_Point(2,FRAME<TV>(TV(0,(T)1,0)));
            curves(0).Add_Control_Point(3,FRAME<TV>(TV(0,(T)2,0)));
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)1,0)*m)),true,true,density,.4*m);
            tests.Add_Ground();
            break;}
        case 8:{
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,4,0)*m));
            tests.Create_Mattress(mattress_grid,true,&initial_state,density);
            tests.Add_Ground();
            break;}
        case 9:{
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)1,0)*m)),true,true,density,.5*m);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0.1,(T)3,0)*m)),true,true,density,.5*m);
            tests.Add_Ground();
            break;}
        case 10:{
            int n=5;
            T drop_height=11*m;
            T horizontal_spacing=2.1*m;
            T vertical_spacing=3.01*m;
            T width=(n-1)*(horizontal_spacing); 
            for(int j=0;j<tori_stack_height;j++)
              for(int i=0;i<tori_stack_width;i++)
                for(int k=0;k<tori_stack_width;k++)
                  tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/adaptive_torus_float.tet",
                      RIGID_BODY_STATE<TV>(FRAME<TV>(TV((-horizontal_spacing*i+width/2),(drop_height+vertical_spacing*j),(-horizontal_spacing*k+width/2)),ROTATION<TV>(-T((i+j+k)*pi/2),TV(0,1,0)))),true,true,density,m);
            T depth=6*m;
            RIGID_BODY<TV>& bowl=tests.Add_Rigid_Body("bowl",depth,(T).5);
            bowl.coefficient_of_friction=input_friction;
            bowl.Frame().r=ROTATION<TV>((T)-pi/2,TV(1,0,0));
            bowl.Frame().t=TV();
            bowl.is_static=true;
            tests.Add_Ground();
            break;}
        case 11:{
            if(!user_last_frame) last_frame=1200;
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,1,0)*m));
            tests.Create_Mattress(mattress_grid,true,&initial_state,density);
            tests.Add_Ground();
            break;}
        case 12:{
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,4,0)*m));
            tests.Create_Mattress(mattress_grid,true,&initial_state,density);
            break;}
        case 13:{
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,4,0)*m));
            tests.Create_Mattress(mattress_grid,true,&initial_state,density);
            tests.Add_Ground();
            RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("bowl",3,(T).5);
            rigid_body.Frame().t.y=1;
            rigid_body.Frame().r=ROTATION<TV>(-pi/2,TV(1,0,0));
            rigid_body.is_static=true;
            break;}
        case 16: {
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,0,0)*m));
            tests.Create_Mattress(mattress_grid,true,&initial_state,density);
            RIGID_BODY<TV>& box1=tests.Add_Analytic_Box(TV(20,20,20)*m);
            RIGID_BODY<TV>& box2=tests.Add_Analytic_Box(TV(20,20,20)*m);
            box1.Frame().t=TV(0,-11,0)*m;
            box2.Frame().t=TV(0,11,0)*m;
            box1.is_static=true;
            box2.is_static=false;
            kinematic_ids.Append(box2.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(box2.particle_index)=true;
            curves.Resize(1);
            curves(0).Add_Control_Point(0,FRAME<TV>(TV(0,11,0)));
            curves(0).Add_Control_Point(5,FRAME<TV>(TV(0,8,0)));
            curves(0).Add_Control_Point(6,FRAME<TV>(TV(0,8,0)));
            curves(0).Add_Control_Point(11,FRAME<TV>(TV(0,11,0)));
            if(!user_last_frame) last_frame=250;
            break;}
        case 17:
        case 18:
        case 24:
        case 25:
        case 26: case 23: case 57:
        case 27:{
            tests.Create_Mattress(mattress_grid,true,0,density);
            break;}
        case 28: {
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,0,0)*m));
            tests.Create_Mattress(mattress_grid,true,&initial_state,density);
            RIGID_BODY<TV>& cylinder1=tests.Add_Analytic_Cylinder(10*m,m);
            cylinder1.Frame().t=TV(0,6,0)*m;
            cylinder1.Frame().r=ROTATION<TV>((T)pi/2.0,TV(1,0,0));
            RIGID_BODY<TV>& cylinder2=tests.Add_Analytic_Cylinder(10*m,m);
            cylinder2.Frame().t=TV(0,-6,0)*m;
            cylinder2.Frame().r=ROTATION<TV>((T)pi/2.0,TV(0,0,1));
            cylinder1.is_static=false;
            cylinder2.is_static=false;
            kinematic_ids.Append(cylinder1.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(cylinder1.particle_index)=true;
            kinematic_ids.Append(cylinder2.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(cylinder2.particle_index)=true;
            curves.Resize(2);
            T period=10.0;T start=0.0;T radius=6.0;int pts=75;T freq=.25;
            for(int ind=0;ind<pts;ind++){
                curves(0).Add_Control_Point(start+freq*ind,FRAME<TV>(TV(radius*sin(2.0*pi*freq*ind/period),radius*cos(2.0*pi*freq*ind/period),0)));
                curves(1).Add_Control_Point(start+freq*ind,FRAME<TV>(TV(-radius*sin(2.0*pi*freq*ind/period),-radius*cos(2.0*pi*freq*ind/period),0)));
            }
            if(!user_last_frame) last_frame=1000;
            break;}
        case 30: {
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)-90,(T)6,(T)11)*m)),true,true,density,1.0*m);
            RIGID_BODY<TV>& box1=tests.Add_Analytic_Box(TV(20,20,20)*m);
            box1.Frame().t=TV(19,10,10)*m;
            box1.is_static=true;
            tests.Add_Ground(1.0);
            if(!user_last_frame) last_frame=480;
            break;}
        case 31: {
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_110K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T)0.373,(T)0)*m)),true,true,density,.005*m);
            RIGID_BODY<TV>& sphere=tests.Add_Analytic_Sphere(0.125*m,density,5);
            sphere.is_static=false;
            sphere.coefficient_of_friction=input_friction;
            tests.Add_Ground();
            kinematic_ids.Append(sphere.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(sphere.particle_index)=true;
            curves.Resize(1);
            curves(0).Add_Control_Point(0,FRAME<TV>(TV(0,0.5,-10)));
            curves(0).Add_Control_Point(20,FRAME<TV>(TV(0,0.5,90)));

            RIGID_BODY<TV>& ball1=tests.Add_Analytic_Sphere(0.065*m,density,6);
            ball1.Frame().t=TV(0.27,0.645,-0.12)*m;
            RIGID_BODY<TV>& ball2=tests.Add_Analytic_Sphere(0.065*m,density,6);
            ball2.Frame().t=TV(-0.275,0.605,-0.18)*m;
            RIGID_BODY<TV>& box1=tests.Add_Analytic_Box(TV(0.4,0.06,0.4)*m);
            box1.Frame().t=TV(-0.024,0.03,0.1)*m;
            RIGID_BODY<TV>& cylinder1=tests.Add_Analytic_Cylinder(0.4*m,0.06*m,196);
            cylinder1.Frame().t=TV(-0.024,0,0.1)+TV(0.2,0,0)*m;
            RIGID_BODY<TV>& cylinder2=tests.Add_Analytic_Cylinder(0.4*m,0.06*m,196);
            cylinder2.Frame().t=TV(-0.024,0,0.1)+TV(-0.2,0,0)*m;
            RIGID_BODY<TV>& cylinder3=tests.Add_Analytic_Cylinder(0.4*m,0.06*m,196);
            cylinder3.Frame().t=TV(-0.024,0,0.1)+TV(0,0,0.2)*m;
            cylinder3.Frame().r=ROTATION<TV>((T)pi/2,TV(0,1,0));
            RIGID_BODY<TV>& cylinder4=tests.Add_Analytic_Cylinder(0.4*m,0.06*m,196);
            cylinder4.Frame().t=TV(-0.024,0,0.1)+TV(0,0,-0.2)*m;
            cylinder4.Frame().r=ROTATION<TV>((T)pi/2,TV(0,1,0));
            RIGID_BODY<TV>& sphere1=tests.Add_Analytic_Sphere(0.06*m,density,6);
            sphere1.Frame().t=TV(-0.024,0,0.1)+TV(0.2,0,0.2)*m;
            RIGID_BODY<TV>& sphere2=tests.Add_Analytic_Sphere(0.06*m,density,6);
            sphere2.Frame().t=TV(-0.024,0,0.1)+TV(0.2,0,-0.2)*m;
            RIGID_BODY<TV>& sphere3=tests.Add_Analytic_Sphere(0.06*m,density,6);
            sphere3.Frame().t=TV(-0.024,0,0.1)+TV(-0.2,0,0.2)*m;
            RIGID_BODY<TV>& sphere4=tests.Add_Analytic_Sphere(0.06*m,density,6);
            sphere4.Frame().t=TV(-0.024,0,0.1)+TV(-0.2,0,-0.2)*m;

            ball1.is_static=true;
            ball2.is_static=true;
            box1.is_static=true;
            cylinder1.is_static=true;
            cylinder2.is_static=true;
            cylinder3.is_static=true;
            cylinder4.is_static=true;
            sphere1.is_static=true;
            sphere2.is_static=true;
            sphere3.is_static=true;
            sphere4.is_static=true;

            break;}
        case 59: case 60: {
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_110K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T)0.373,(T)0)*m)),true,true,density,.005*m);
            break;}

        case 32:{
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_thin_24K.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,0,0)*m)),true,true,density,m);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_thin_24K.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(-2.8,0,0)*m,ROTATION<TV>(T(pi/2),TV(1,0,0)))),true,true,density,m);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_thin_24K.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(2.8,0,0)*m,ROTATION<TV>(T(pi/2),TV(1,0,0)))),true,true,density,m);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_thin_24K.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(-5.6,0,0)*m)),true,true,density,m);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_thin_24K.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(5.6,0,0)*m)),true,true,density,m);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_thin_24K.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(-8.4,0,0)*m,ROTATION<TV>(T(pi/2),TV(1,0,0)))),true,true,density,m);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_thin_24K.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(8.4,0,0)*m,ROTATION<TV>(T(pi/2),TV(1,0,0)))),true,true,density,m);

            RIGID_BODY<TV>& torus1=tests.Add_Analytic_Torus((T).5,(T)2,32,64);
            torus1.is_static=false;
            torus1.coefficient_of_friction=0.05;

            RIGID_BODY<TV>& torus2=tests.Add_Analytic_Torus((T).5,(T)2,32,64);
            torus2.is_static=false;
            torus2.coefficient_of_friction=0.05;

            kinematic_ids.Append(torus1.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(torus1.particle_index)=true;
            kinematic_ids.Append(torus2.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(torus2.particle_index)=true;

            T x_start=11.2;
            T x_stop=15;
            T t_stop=2;
            T t_rot=8;

            curves.Resize(2);
            curves(0).Add_Control_Point(0,FRAME<TV>(TV(x_start,0,0)*m,ROTATION<TV>((T)0,TV(1,0,0))));
            curves(0).Add_Control_Point(t_stop,FRAME<TV>(TV(x_stop,0,0)*m,ROTATION<TV>((T)0,TV(1,0,0))));

            for(int i=0;i<47;i++) curves(0).Add_Control_Point(t_rot+i*2,FRAME<TV>(TV(x_stop,0,0)*m,ROTATION<TV>((T)pi/2.0*i,TV(1,0,0))));

            curves(1).Add_Control_Point(0,FRAME<TV>(TV(-x_start,0,0)*m,ROTATION<TV>((T)0,TV(1,0,0))));
            curves(1).Add_Control_Point(t_stop,FRAME<TV>(TV(-x_stop,0,0)*m,ROTATION<TV>((T)0,TV(1,0,0))));

            for(int i=0;i<47;i++) curves(1).Add_Control_Point(t_rot+i*2,FRAME<TV>(TV(-x_stop,0,0)*m,ROTATION<TV>((T)-pi/2.0*i,TV(1,0,0))));

            break;}
        case 33:{
            T y_translate=0.1;

            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_110K.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(-0.2,2.5-y_translate,-0.04)*m,ROTATION<TV>(T(-pi/2.*0.9),TV(0,0,1))*ROTATION<TV>(T(pi/2.),TV(0,1,0)))),true,true,density,0.013*m);

            RIGID_BODY<TV>& gear1=tests.Add_Rigid_Body("gear",.375*m,1);
            RIGID_BODY<TV>& gear2=tests.Add_Rigid_Body("gear",.375*m,1);
            RIGID_BODY<TV>& cylinder=tests.Add_Analytic_Cylinder(1.5*m,.15*m,24);

            gear1.coefficient_of_friction=1;
            gear2.coefficient_of_friction=1;
            cylinder.coefficient_of_friction=0;

            kinematic_ids.Append(gear1.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(gear1.particle_index)=true;
            kinematic_ids.Append(gear2.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(gear2.particle_index)=true;
            kinematic_ids.Append(cylinder.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(cylinder.particle_index)=true;

            T angular_velocity=1;
            T scale=1;

            curves.Resize(3);
            curves(2).Add_Control_Point(0,FRAME<TV>(TV(0,4-y_translate,0)*m,ROTATION<TV>(0,TV(0,0,1))));
            curves(2).Add_Control_Point(2,FRAME<TV>(TV(0,4-y_translate,0)*m,ROTATION<TV>(0,TV(0,0,1))));
            curves(2).Add_Control_Point(4,FRAME<TV>(TV(0,2-y_translate,0)*m,ROTATION<TV>(0,TV(0,0,1))));

            for(int i=0;i<60;i++){
                curves(0).Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(-(T).4*scale,1.5*scale-y_translate,-.75*scale)*m,ROTATION<TV>(-i,TV(0,0,1))));
                curves(1).Add_Control_Point(i/angular_velocity,FRAME<TV>(TV((T).4*scale,1.5*scale-y_translate,-.75*scale)*m,ROTATION<TV>(i,TV(0,0,1))));}

            tests.Add_Ground();
            break;}

        case 58:{
            T scale=0.75;

            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/hand_30k.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3*scale,0)*m,ROTATION<TV>(T(pi/2),TV(0,1,0))*ROTATION<TV>(T(pi/2),TV(1,0,0)))),true,true,density,hand_scale*m);
            if(!gears_of_pain){ //tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/fish_42K.tet",
                //  RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,4.5*scale,-1.2*scale)*m,ROTATION<TV>((T)pi*0.525,TV(1,0,0))*ROTATION<TV>(0*(T)pi/2,TV(0,1,0)))),true,true,density,0.06*m);
                //tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/bunny.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)4.7*scale,-3.0*scale)*m)),true,true,density,.25*m);
                tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_110K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)4.6*scale,-0.0*scale)*m,ROTATION<TV>(T(-pi/2),TV(0,0,1))*ROTATION<TV>(T(pi/2),TV(0,1,0)))),true,true,density,.0065*m);
            }else{
                tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_110K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)4.6*scale,-0.0*scale)*m,ROTATION<TV>(T(-pi/2),TV(0,0,1))*ROTATION<TV>(T(pi/2),TV(0,1,0)))),true,true,density,.005*m);
            }


            RIGID_BODY<TV>& gear1=tests.Add_Rigid_Body("gear",.375*scale*m,1);
            RIGID_BODY<TV>& gear2=tests.Add_Rigid_Body("gear",.375*scale*m,1);

            gear1.coefficient_of_friction=input_friction;
            gear2.coefficient_of_friction=input_friction;
            RIGID_BODY<TV>& gear3=tests.Add_Rigid_Body("gear",.375*scale*m,1);
            RIGID_BODY<TV>& gear4=tests.Add_Rigid_Body("gear",.375*scale*m,1);

            gear3.coefficient_of_friction=input_friction;
            gear4.coefficient_of_friction=input_friction;
            RIGID_BODY<TV>& gear5=tests.Add_Rigid_Body("gear",.375*scale*m,1);
            RIGID_BODY<TV>& gear6=tests.Add_Rigid_Body("gear",.375*scale*m,1);

            gear5.coefficient_of_friction=input_friction;
            gear6.coefficient_of_friction=input_friction;

            kinematic_ids.Append(gear1.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(gear1.particle_index)=true;
            kinematic_ids.Append(gear2.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(gear2.particle_index)=true;
            kinematic_ids.Append(gear3.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(gear3.particle_index)=true;
            kinematic_ids.Append(gear4.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(gear4.particle_index)=true;
            kinematic_ids.Append(gear5.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(gear5.particle_index)=true;
            kinematic_ids.Append(gear6.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(gear6.particle_index)=true;

            T angular_velocity=1;T gear_dx;
            if (gears_of_pain) gear_dx=(T).377;else gear_dx=(T).4;

            curves.Resize(8);
            for(int i=0;i<60;i++){
                curves(0).Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(-gear_dx*scale,1.5*scale,-.75*scale)*m,ROTATION<TV>(-i,TV(0,0,1))));//.4 is default, gears barely touch at .375
                curves(1).Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(gear_dx*scale,1.5*scale,-.75*scale)*m,ROTATION<TV>(i,TV(0,0,1))));
                curves(2).Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(-gear_dx*scale,1.5*scale,-1.75*scale)*m,ROTATION<TV>(-i,TV(0,0,1))));//.4 is default, gears barely touch at .375
                curves(3).Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(gear_dx*scale,1.5*scale,-1.75*scale)*m,ROTATION<TV>(i,TV(0,0,1))));
                curves(4).Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(-gear_dx*scale,1.5*scale,.25*scale)*m,ROTATION<TV>(-i,TV(0,0,1))));//.4 is default, gears barely touch at .375
                curves(5).Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(gear_dx*scale,1.5*scale,.25*scale)*m,ROTATION<TV>(i,TV(0,0,1))));
            }

            RIGID_BODY<TV>& box0=tests.Add_Analytic_Box(TV(2.0*scale,.1*scale,2.0*scale)*m);

            //  RIGID_BODY<TV>& box2=tests.Add_Analytic_Box(TV(2.0*scale,2.0*scale,.1*scale)*m);
            //  RIGID_BODY<TV>& box3=tests.Add_Analytic_Box(TV(2.0*scale,.1*scale,2.0*scale)*m);
            RIGID_BODY<TV>& cylinder=tests.Add_Analytic_Cylinder(1.5*scale*m,.12*scale*m);
            box0.Frame().t=TV(0,4.0*scale,-0.0*scale)*m;

            /*  if(!gears_of_pain){
             RIGID_BODY<TV>& box1=tests.Add_Analytic_Box(TV(2.0*scale,2.0*scale,.1*scale)*m);
             box1.Frame().t=TV(0,1.2*scale,0.831*scale)*m;
             //box1.Frame().r=ROTATION<TV>((T)pi/4.0,TV(1,0,0));
             box1.is_static=true;
             }*/

            //box2.Frame().t=TV(0,1.2*scale,-0.831*scale)*m;
            // box2.Frame().r=ROTATION<TV>(-(T)pi/4.0,TV(1,0,0));
            // box3.Frame().t=TV(0,3.0*scale,-1.4*scale)*m;

            // box3.Frame().r=ROTATION<TV>((T)pi/4.0,TV(1,0,0));
            box0.is_static=false;//Will move later
            // box2.is_static=true;
            // box3.is_static=false;

            T drop_time;
            if(gears_of_pain) drop_time=(T)17;else drop_time=(T)17;
            cylinder.Frame().t=TV(0,5.0*scale,0*scale)*m;
            cylinder.is_static=false;
            kinematic_ids.Append(cylinder.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(cylinder.particle_index)=true;
            curves(6).Add_Control_Point(0,FRAME<TV>(TV(0,5.0*scale,0*scale)));
            curves(6).Add_Control_Point(drop_time,FRAME<TV>(TV(0,5.0*scale,0*scale)));
            curves(6).Add_Control_Point(drop_time+(T)1,FRAME<TV>(TV(0,2.0*scale,0*scale)));

            box0.coefficient_of_friction=.05;
            kinematic_ids.Append(box0.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(box0.particle_index)=true;
            curves(7).Add_Control_Point(0,FRAME<TV>(TV(0,4.0*scale,-0.0*scale)));
            curves(7).Add_Control_Point(1.5*s+hand_scale*s,FRAME<TV>(TV(0,4.0*scale,-0.0*scale)*m,ROTATION<TV>(0*(T)pi/2.0,TV(1,0,0))));
            curves(7).Add_Control_Point(1.53*s+hand_scale*s,FRAME<TV>(TV(0,3.1*scale,-1.0*scale)*m,ROTATION<TV>((T)pi/2.0,TV(1,0,0))));

            RIGID_BODY<TV>& inclined_floor=tests.Add_Ground(input_friction);
            inclined_floor.Frame().t=TV(0,.0*scale,0)*m;
            break;}

        case 34:{
            T radius=8.0;
            T velocity=7;

            RIGID_BODY<TV>& cylinder=tests.Add_Analytic_Cylinder((T)32*m,radius*m,64);
            cylinder.is_static=false;
            cylinder.coefficient_of_friction=1e8;
            kinematic_ids.Append(cylinder.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(cylinder.particle_index)=true;
            curves.Resize(1);
            for(int i=0;i<128;i++) curves(0).Add_Control_Point(i,FRAME<TV>(TV(-25+i*velocity,radius*1.05,0)*m,ROTATION<TV>(-i*velocity/radius,TV(0,0,1))));

            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,2,0)*m));
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,3,0)*m)),true,true,density,3*m);
            // tests.Create_Mattress(mattress_grid,true,&initial_state,density);
            tests.Add_Ground(1e8);
            break;}
        case 35:{
             RIGID_BODY_STATE<TV> initial_state1(FRAME<TV>(TV(8,7,.5)*m,ROTATION<TV>(T(pi/3),TV(1,0,1))));
             RIGID_BODY_STATE<TV> initial_state2(FRAME<TV>(TV(8,1.2,0)*m,ROTATION<TV>(T(pi/6),TV(0,1,0))));
            RIGID_BODY_STATE<TV> initial_state3(FRAME<TV>(TV(-.1,6,0)*m,ROTATION<TV>(T(pi/4),TV(1,1,1))));
            RIGID_BODY_STATE<TV> initial_state4(FRAME<TV>(TV(0,1.5,0)*m,ROTATION<TV>(T(0),TV(1,1,1))));

            tests.Create_Mattress(mattress_grid1,true,&initial_state1,density);
            tests.Create_Mattress(mattress_grid2,true,&initial_state2,density);
            tests.Create_Mattress(mattress_grid2,true,&initial_state3,density);
            tests.Create_Mattress(mattress_grid3,true,&initial_state4,density);
            tests.Add_Ground();
            break;
        }
        case 36:{
            RIGID_BODY_STATE<TV> initial_state1(FRAME<TV>(TV(0,8,0)*m,ROTATION<TV>(T(0),TV(1,0,1))));
            tests.Create_Mattress(mattress_grid1,true,&initial_state1,density);
            tests.Add_Ground();
            break;
        }
        case 37:
        {
            RIGID_BODY_STATE<TV> initial_state1(FRAME<TV>(TV(1.7,5.2,1.6)*m,ROTATION<TV>(T(pi/4),TV(1.3,0.3,0.7))));
            RIGID_BODY_STATE<TV> initial_state2(FRAME<TV>(TV(0,2.2,0)*m,ROTATION<TV>(T(pi/7),TV(0.5,2,0.3))));
            tests.Create_Mattress(mattress_grid,true,&initial_state1,density);
            tests.Create_Mattress(mattress_grid,true,&initial_state2,density);
            tests.Add_Ground();
            break;
        }
        case 38:
        {
            number_of_jellos=13;

            for(int i=0;i<number_of_jellos;i++)
            {
                jello_centers.Append(TV(5.3*sin(277*i),15+2*cos(123*i),5.3*cos(297*i)));
                RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(jello_centers(i)*m,ROTATION<TV>(10*sin(178*i),TV(sin(145*i),cos(345*i),cos(478*i)))));
                tests.Create_Mattress(mattress_grid,true,&initial_state,density);
            }
            tests.Add_Ground();
            break;
        }
        case 39:
        {
            RIGID_BODY_STATE<TV> initial_state1(FRAME<TV>(TV(0.07,0.3,0.09)*m,ROTATION<TV>(T(pi/0.103),TV(1.35,0.785,1.675))));
            RIGID_BODY_STATE<TV> initial_state2(FRAME<TV>(TV(0,0.01,0)*m));
            tests.Create_Mattress(mattress_grid,true,&initial_state1,density);
            tests.Create_Mattress(mattress_grid,true,&initial_state2,density);
            tests.Add_Ground();
            break;
        }
        case 40:
        {
            jello_centers.Append(TV(-0.266,0.022,0.013));
            jello_centers.Append(TV(0.266, 0.029,-0.013));
            RIGID_BODY_STATE<TV> initial_state1(FRAME<TV>(jello_centers(0)*m,ROTATION<TV>(T(pi/0.13),TV(1.3,1.5,0.7))));
            RIGID_BODY_STATE<TV> initial_state2(FRAME<TV>(jello_centers(1)*m,ROTATION<TV>(T(pi/0.076),TV(0.7,1,0.1))));
            tests.Create_Mattress(mattress_grid,true,&initial_state1,density);
            tests.Create_Mattress(mattress_grid,true,&initial_state2,density);
            tests.Add_Ground();
            break;
        }
        case 41:

        {
            //Todo: At some point, use Oriented_Box to do a tighter collision check. C-Re, C-e, R where R is the rotation matrix and e is the edge
            T max_jello_size=.036;//maximum edge length
            T bound=.2;
            TV new_center;T new_rotate;
            T board_height=.6-.1*bound;plateau=board_height;
            bool stuck=false;
            RIGID_BODY_STATE<TV> initial_state;

            //break;}
            for(int i=0;i<number_of_jellos;i++){
                do {
                    rand.Fill_Uniform(new_center,-bound,bound);
                   // new_center=TV(rand.Get_Uniform_Number(-bound,bound),rand.Get_Uniform_Number((T).5*bound,(T)1*bound),rand.Get_Uniform_Number(-bound,bound));
                    stuck=false;
                    new_center.x=.6*(new_center.x);
                    new_center.z=1.1*(new_center.z+bound)-bound;
                    if (new_center.z < -.5*bound) new_center.y=(T).5*(new_center.y+2.0*bound+3.0*(new_center.z+bound));
                    else new_center.y=.3*(new_center.y+bound)+.5*max_jello_size+board_height+.2*bound;
                    for(int j=0;j<i&&(!stuck);j++){
                        if((new_center-jello_centers(j)).Magnitude()<=(T)2*max_jello_size) stuck=true;
                }}while(stuck);
                jello_centers.Append(new_center);
                rand.Fill_Uniform(new_center,-bound,bound);
                new_rotate=rand.Get_Uniform_Number(-(T)pi,(T)pi);
                RIGID_BODY_STATE<TV> initial_state1(FRAME<TV>(jello_centers(i)*m,ROTATION<TV>(new_rotate,new_center)));
                if (i % 5 ==0) {tests.Create_Mattress(mattress_grid3,true,&initial_state1,density);}
                else if (i % 5 ==1) {tests.Create_Mattress(mattress_grid2,true,&initial_state1,density);}
                else {tests.Create_Mattress(mattress_grid1,true,&initial_state1,density);}
            }
            if(sloped_floor){
            RIGID_BODY<TV>& inclined_floor=tests.Add_Ground(input_friction);
            inclined_floor.Frame().r=ROTATION<TV>((T)pi*degrees_wedge/(T)180,TV(0,0,1))*ROTATION<TV>((T)pi*degrees_incline/(T)180,TV(1,0,0));
            RIGID_BODY<TV>& inclined_floor2=tests.Add_Ground(input_friction);
            inclined_floor2.Frame().r=ROTATION<TV>(-(T)pi*degrees_wedge/(T)180,TV(0,0,1))*ROTATION<TV>((T)pi*degrees_incline/(T)180,TV(1,0,0));
            }
            else
            {
            RIGID_BODY<TV>& inclined_floor=tests.Add_Ground(input_friction);
            inclined_floor.Frame().r=ROTATION<TV>((T)pi*degrees_incline/(T)180,TV(1,0,0));
            }

            T dy=.3*bound*sin((T)pi*degrees_incline/(T)180);
            RIGID_BODY<TV>& box1=tests.Add_Analytic_Box(TV(2.2*bound,.25*bound,.4*bound)*m);
            RIGID_BODY<TV>& box2=tests.Add_Analytic_Box(TV(2.2*bound,.25*bound,.2*bound)*m);
            RIGID_BODY<TV>& box3=tests.Add_Analytic_Box(TV(2.2*bound,.25*bound,.2*bound)*m);
            RIGID_BODY<TV>& box4=tests.Add_Analytic_Box(TV(2.2*bound,.25*bound,.2*bound)*m);
            RIGID_BODY<TV>& box5=tests.Add_Analytic_Box(TV(2.2*bound,.25*bound,.2*bound)*m);
            RIGID_BODY<TV>& box6=tests.Add_Analytic_Box(TV(2.2*bound,.25*bound,.2*bound)*m);
            RIGID_BODY<TV>& box7=tests.Add_Analytic_Box(TV(2.2*bound,.25*bound,.2*bound)*m);
            RIGID_BODY<TV>& box8=tests.Add_Analytic_Box(TV(2.2*bound,.25*bound,.4*bound)*m);
            box1.Frame().t=TV((T)0,board_height,-.4*bound);
            box2.Frame().t=TV((T)0*bound,board_height,-0.1*bound);
            box3.Frame().t=TV((T)0*bound,board_height,0.1*bound);
            box4.Frame().t=TV((T)0*bound,board_height,0.3*bound);
            box5.Frame().t=TV((T)0*bound,board_height,0.5*bound);
            box6.Frame().t=TV((T)0*bound,board_height,0.7*bound);
            box7.Frame().t=TV((T)0*bound,board_height,0.9*bound);
            box8.Frame().t=TV((T)0*bound,board_height,1.2*bound);

            curves.Resize(8);
            curves(0).Add_Control_Point(0,FRAME<TV>(box1.Frame().t));
            curves(0).Add_Control_Point(.2,FRAME<TV>(box1.Frame().t));
            curves(1).Add_Control_Point(0,FRAME<TV>(box2.Frame().t));
            curves(1).Add_Control_Point(.2,FRAME<TV>(box2.Frame().t));
            curves(1).Add_Control_Point(.27,FRAME<TV>(box2.Frame().t-TV(0,dy,0)));
            curves(2).Add_Control_Point(0,FRAME<TV>(box3.Frame().t));
            curves(2).Add_Control_Point(.2,FRAME<TV>(box3.Frame().t));
            curves(2).Add_Control_Point(.33,FRAME<TV>(box3.Frame().t-TV(0,2.0*dy,0)));
            curves(3).Add_Control_Point(0,FRAME<TV>(box4.Frame().t));
            curves(3).Add_Control_Point(.2,FRAME<TV>(box4.Frame().t));
            curves(3).Add_Control_Point(.40,FRAME<TV>(box4.Frame().t-TV(0,3.0*dy,0)));
            curves(4).Add_Control_Point(0,FRAME<TV>(box5.Frame().t));
            curves(4).Add_Control_Point(.2,FRAME<TV>(box5.Frame().t));
            curves(4).Add_Control_Point(.47,FRAME<TV>(box5.Frame().t-TV(0,4.0*dy,0)));
            curves(5).Add_Control_Point(0,FRAME<TV>(box6.Frame().t));
            curves(5).Add_Control_Point(.2,FRAME<TV>(box6.Frame().t));
            curves(5).Add_Control_Point(.53,FRAME<TV>(box6.Frame().t-TV(0,5.0*dy,0)));
            curves(6).Add_Control_Point(0,FRAME<TV>(box7.Frame().t));
            curves(6).Add_Control_Point(.2,FRAME<TV>(box7.Frame().t));
            curves(6).Add_Control_Point(.60,FRAME<TV>(box7.Frame().t-TV(0,6.0*dy,0)));
            curves(7).Add_Control_Point(0,FRAME<TV>(box8.Frame().t));
            curves(7).Add_Control_Point(.2,FRAME<TV>(box8.Frame().t));
            curves(7).Add_Control_Point(.67,FRAME<TV>(box8.Frame().t-TV(0,7.0*dy,0)));
            box1.is_static=false;
            kinematic_ids.Append(box1.particle_index);
            box2.is_static=false;
            kinematic_ids.Append(box2.particle_index);
            box3.is_static=false;
            kinematic_ids.Append(box3.particle_index);
            box4.is_static=false;
            kinematic_ids.Append(box4.particle_index);
            box5.is_static=false;
            kinematic_ids.Append(box5.particle_index);
            box6.is_static=false;
            kinematic_ids.Append(box6.particle_index);
            box7.is_static=false;
            kinematic_ids.Append(box7.particle_index);
            box8.is_static=false;
            kinematic_ids.Append(box8.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(box1.particle_index)=true;
            rigid_body_collection.rigid_body_particles.kinematic(box2.particle_index)=true;
            rigid_body_collection.rigid_body_particles.kinematic(box3.particle_index)=true;
            rigid_body_collection.rigid_body_particles.kinematic(box4.particle_index)=true;
            rigid_body_collection.rigid_body_particles.kinematic(box5.particle_index)=true;
            rigid_body_collection.rigid_body_particles.kinematic(box6.particle_index)=true;
            rigid_body_collection.rigid_body_particles.kinematic(box7.particle_index)=true;
            rigid_body_collection.rigid_body_particles.kinematic(box8.particle_index)=true;
            break;
        }
        case 42:
        {
            int count=0;
            for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
            for(int k=0;k<3;k++)
            {
                count++;
                jello_centers.Append(TV(-100+i*5,j*5+3,k*5+sin(75*count)));
                RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(jello_centers(count)*m,ROTATION<TV>(10*sin(178*count),TV(sin(145*count),cos(345*count),cos(478*count)))));
                tests.Create_Mattress(mattress_grid,true,&initial_state,density);
            }
            tests.Add_Ground();
            break;
        }
        case 52:
        {
            T flat_bottom_radius=0.02;
            T rounded_interior_radius=0.02;
            T vertical_straight_length=0.12;
            T thickness=0.01;
            T coefficient_of_friction=0.3;
            RIGID_BODY<TV>& rounding=tests.Add_Analytic_Bowl(TV(),TV(0,1,0),flat_bottom_radius,rounded_interior_radius,thickness,128,32);
            rounding.coefficient_of_friction=coefficient_of_friction;
            rounding.Frame().t=TV(0,thickness,0);
            rounding.is_static=true;

            RIGID_BODY<TV>& walls=tests.Add_Analytic_Shell(vertical_straight_length,flat_bottom_radius+rounded_interior_radius+thickness,flat_bottom_radius+rounded_interior_radius,128);
            walls.coefficient_of_friction=coefficient_of_friction;
            walls.Frame().r=ROTATION<TV>((T)pi/2,TV(1,0,0));
            walls.Frame().t=TV(0,vertical_straight_length/2,0);
            walls.is_static=true;

            RIGID_BODY<TV>& bottom=tests.Add_Analytic_Cylinder(thickness*m,(flat_bottom_radius+rounded_interior_radius+thickness)*m,128);
            bottom.coefficient_of_friction=coefficient_of_friction;
            bottom.Frame().r=ROTATION<TV>((T)pi/2,TV(1,0,0));
            bottom.Frame().t=TV(0,thickness/2,0);
            bottom.is_static=true;

            int count=0;
            for(int i=0;i<25;i++)
            {
                count++;
                jello_centers.Append(TV(i*0.05,50,0));
                RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(jello_centers(count)*m,ROTATION<TV>(10*sin(178*i),TV(sin(145*i),cos(345*i),cos(478*i)))));
                tests.Create_Mattress(mattress_grid,true,&initial_state,density);
            }

            tests.Add_Ground();
            break;
        }
        case 44:
        {
            jello_centers.Append(TV(-10,3.2,0));
            jello_centers.Append(TV(10,3.5,-1));
            RIGID_BODY_STATE<TV> initial_state1(FRAME<TV>(jello_centers(0)*m,ROTATION<TV>(T(pi/4),TV(1.3,1.5,0.7))));
            RIGID_BODY_STATE<TV> initial_state2(FRAME<TV>(jello_centers(1)*m,ROTATION<TV>(T(pi/5),TV(0.7,1,0.1))));
            tests.Create_Mattress(mattress_grid,true,&initial_state1,density);
            tests.Create_Mattress(mattress_grid,true,&initial_state2,density);
            tests.Add_Ground();
            break;
        }
        case 50:
        {
            RIGID_BODY<TV>& torus1=tests.Add_Analytic_Torus((T)1.4,(T)1.6,32,64);
            torus1.is_static=true;
            torus1.coefficient_of_friction=0.00;
            torus1.Frame().t=TV(10,10,0)*m;
            torus1.Frame().r=ROTATION<TV>((T)pi/2.0,TV(0,1,0));
            if(!user_last_frame) last_frame=240;
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/fish_42K.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,10,0)*m,ROTATION<TV>(T(pi),TV(0,1,0)))),true,true,density,m);
            tests.Add_Ground(1);
            break;
        }
        case 47:
        {
            fishes=5;
            for(int i=0;i<fishes;i++)
                tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/fish_42K.tet",
                    RIGID_BODY_STATE<TV>(FRAME<TV>(TV(-10*i,10,0)*m,ROTATION<TV>(T(pi),TV(0,1,0)))),true,true,density,m);
            T boxsize=(T)2;
            RIGID_BODY<TV>& box1=tests.Add_Analytic_Box(TV(3*boxsize,boxsize,boxsize)*m);
            RIGID_BODY<TV>& box2=tests.Add_Analytic_Box(TV(3*boxsize,boxsize,boxsize)*m);
            RIGID_BODY<TV>& box3=tests.Add_Analytic_Box(TV(boxsize,boxsize,boxsize)*m);
            RIGID_BODY<TV>& box4=tests.Add_Analytic_Box(TV(boxsize,boxsize,boxsize)*m);

            box1.Frame().t=TV(0,.5*boxsize,boxsize);
            box2.Frame().t=TV(0,.5*boxsize,-boxsize);
            box3.Frame().t=TV(-boxsize,.5*boxsize,0);
            box4.Frame().t=TV(boxsize,.5*boxsize,0);

            box1.is_static=true;
            box2.is_static=true;
            box3.is_static=true;
            box4.is_static=true;
            tests.Add_Ground();
            break;
        }
        case 49:
        {
            TV start(10*m,10*m,0);
            T outer=(T)1*m,inner=(T).5*m,length=10*m;
            RIGID_BODY<TV>& torus1=tests.Add_Analytic_Torus((outer-inner)/2,(outer+inner)/2,32,64);
            //torus1.is_static=true;
            torus1.coefficient_of_friction=0.05;
            torus1.Frame().t=start;
            torus1.Frame().r=ROTATION<TV>((T)pi/2.0,TV(0,1,0));
            if(!user_last_frame) last_frame=240;
            RIGID_BODY<TV>& torus2=tests.Add_Analytic_Torus((outer-inner)/2,(outer+inner)/2,32,64);
            //torus2.is_static=true;
            torus2.coefficient_of_friction=0.05;
            torus2.Frame().t=start+TV(length,0,0);
            torus2.Frame().r=ROTATION<TV>((T)pi/2.0,TV(0,1,0));
           // Add_Joint(1,2,new POINT_JOINT<TV>,FRAME<TV>(TV(15,10,0)));// ropes

            /*if(!user_last_frame) last_frame=240;
            RIGID_BODY<TV>& shell=tests.Add_Analytic_Shell(length,outer,inner,64);
            shell.is_static=true;
            shell.coefficient_of_friction=0.05;
            shell.Frame().t=start+TV(length/2,0,0);
            shell.Frame().r=ROTATION<TV>((T)pi/2.0,TV(0,1,0));
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/hand_30k.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,10,0)*m,ROTATION<TV>(T(pi),TV(0,1,0)))),true,true,density,m);*/

            //Add_Rigid_Body("plank",(T).125*m,(T).5,0,FRAME<TV>(TV(0,3,0)),"seesaw");
            tests.Add_Ground(1);
            break;
        }
        case 51:
        {
            TV start(10*m,10*m,0);
            T outer=(T)2.5*m,inner=hole,length=10*m;
            RIGID_BODY<TV>& torus1=tests.Add_Analytic_Torus((outer-inner)/2,(outer+inner)/2,32,64);
            torus1.is_static=true;
            torus1.coefficient_of_friction=0.05;
            torus1.Frame().t=start;
            torus1.Frame().r=ROTATION<TV>((T)pi/2.0,TV(0,1,0));
            if(!user_last_frame) last_frame=240;
            RIGID_BODY<TV>& torus2=tests.Add_Analytic_Torus((outer-inner)/2,(outer+inner)/2,32,64);
            torus2.is_static=true;
            torus2.coefficient_of_friction=0.05;
            torus2.Frame().t=start+TV(length,0,0);
            torus2.Frame().r=ROTATION<TV>((T)pi/2.0,TV(0,1,0));
            if(!user_last_frame) last_frame=240;
            RIGID_BODY<TV>& shell=tests.Add_Analytic_Shell(length,outer,inner,64);
            shell.is_static=true;
            shell.coefficient_of_friction=0.05;
            shell.Frame().t=start+TV(length/2,0,0);
            shell.Frame().r=ROTATION<TV>((T)pi/2.0,TV(0,1,0));
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/fish_42K.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,10,0)*m,ROTATION<TV>(T(pi),TV(0,1,0)))),true,true,density,m);
            tests.Add_Ground(1);
            break;
        }
        case 43:{
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/fish_42K.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,10,0)*m,ROTATION<TV>(T(0),TV(0,1,0)))),true,true,density,m);
            int cogs=40;
            T radius=.8;
            T sep=radius+.003;
            T cog_rad=.035;
            RIGID_BODY<TV>& gear1=tests.Add_Analytic_Smooth_Gear(TV(radius,cog_rad,1),cogs,8);
            RIGID_BODY<TV>& gear2=tests.Add_Analytic_Smooth_Gear(TV(radius,cog_rad,1),cogs,8);
            gear1.coefficient_of_friction=2;
            gear2.coefficient_of_friction=2;
            kinematic_ids.Append(gear1.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(gear1.particle_index)=true;
            kinematic_ids.Append(gear2.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(gear2.particle_index)=true;
            T angular_velocity=1;

            curves.Resize(2);
            for(int i=0;i<60;i++){
                curves(0).Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(-sep,1.5,0)*m,ROTATION<TV>(-i,TV(0,0,1))*ROTATION<TV>(pi/(2*cogs),TV(0,0,1))));
                curves(1).Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(sep,1.5,0)*m,ROTATION<TV>(i,TV(0,0,1))*ROTATION<TV>(pi/(2*cogs),TV(0,0,1))));}

            tests.Add_Ground();
            break;}
        case 48:{
            tests.Create_Mattress(GRID<TV>(TV_INT()+(resolution?resolution:10)+1,RANGE<TV>::Centered_Box()*m),true,0,density);
            break;}
        case 53:{
            TETRAHEDRALIZED_VOLUME<T>* tv=TETRAHEDRALIZED_VOLUME<T>::Create(particles);
            int a=particles.Add_Element();
            int b=particles.Add_Element();
            int c=particles.Add_Element();
            int d=particles.Add_Element();
            particles.X(a)=TV(1,0,0);
            particles.X(b)=TV(0,-1,0);
            particles.X(c)=TV(0,.5,.5);
            particles.X(d)=TV(0,.5,-.5);
            particles.mass(a)=1;
            particles.mass(b)=1;
            particles.mass(c)=1;
            particles.mass(d)=1;
            tv->mesh.elements.Append(VECTOR<int,4>(a,b,c,d));
            solid_body_collection.deformable_body_collection.Add_Structure(tv);
            break;}
        case 54: case 55:{
            TETRAHEDRALIZED_VOLUME<T>* tv=TETRAHEDRALIZED_VOLUME<T>::Create(particles);
            int a=particles.Add_Element();
            int b=particles.Add_Element();
            int c=particles.Add_Element();
            int d=particles.Add_Element();
            particles.X(a)=TV(1,0,0);
            particles.X(b)=TV(0,0,0);
            particles.X(c)=TV(0,0,1);
            particles.X(d)=TV(0,1,0);
            particles.mass(a)=1;
            particles.mass(b)=1;
            particles.mass(c)=1;
            particles.mass(d)=1;
            tv->mesh.elements.Append(VECTOR<int,4>(a,b,c,d));
            solid_body_collection.deformable_body_collection.Add_Structure(tv);
            break;}
        case 56:{
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)5,0)*m)),true,true,density,.5*m);
            tests.Add_Ground();
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/fish_42K.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(5,5,0)*m,ROTATION<TV>(T(pi),TV(0,1,0)))),true,true,density,.2*m);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_110K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(10,5,0)*m,ROTATION<TV>(T(pi),TV(0,1,0)))),true,true,density,.005*m);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/hand_30k.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(15,5,0)*m,ROTATION<TV>(-T(pi/2),TV(1,0,0)))),true,true,density,m);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/bunny.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(20,(T)5,0)*m)),true,true,density,m);
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(25,5,0)*m));
            tests.Create_Mattress(mattress_grid,true,&initial_state,density);
            break;}
        case 61:{
            int num_sides=resolution?resolution:3;
            tests.Create_Cylinder(CYLINDER<T>(TV(-rod_length/2,0,0),TV(rod_length/2,0,0),rod_radius),(int)ceil(num_sides/rod_radius*rod_length),num_sides,true,0,1000);
            for(int i=0;i<particles.X.m;i++){
                T x=particles.X(i).x;
                ARRAY<int>& ar=(x>=0)?constrained_particles:externally_forced;
                x=abs(x);
                if(x>rod_length/2-attachment_length) ar.Append(i);
                if(x>attachment_length && particles.X(i).Projected_Orthogonal_To_Unit_Direction(TV(1,0,0)).Magnitude() < .1) ar.Append(i);}
            for(int i=0;i<constrained_particles.m;i++) Add_Debug_Particle(particles.X(constrained_particles(i)),TV(1,0,0));
            for(int i=0;i<externally_forced.m;i++) Add_Debug_Particle(particles.X(externally_forced(i)),TV(0,1,0));
            initial_positions=particles.X;
            constrained_velocities.Resize(constrained_particles.m);
            scalar_curve.Add_Control_Point(0,0);
            scalar_curve.Add_Control_Point(4,-.85*pi);
            scalar_curve.Add_Control_Point(5,-.85*pi);
            scalar_curve.Add_Control_Point(9,0);
            break;}
        case 62:{
            tests.Create_Mattress(mattress_grid,true,0,density);
            for(int i=0;i<particles.X.m;i++){
                T x=particles.X(i).x;
                ARRAY<int>& ar=(x>=0)?constrained_particles:externally_forced;
                x=abs(x);
                if(x>4.99) ar.Append(i);}
            for(int i=0;i<constrained_particles.m;i++) Add_Debug_Particle(particles.X(constrained_particles(i)),TV(1,0,0));
            for(int i=0;i<externally_forced.m;i++) Add_Debug_Particle(particles.X(externally_forced(i)),TV(0,1,0));
            initial_positions=particles.X;
            constrained_velocities.Resize(constrained_particles.m);
            scalar_curve.Add_Control_Point(0,0);
            scalar_curve.Add_Control_Point(9,4*pi);
            break;}
        case 63:{
            tests.Create_Mattress(mattress_grid,true,0,density);
            for(int i=0;i<particles.X.m;i++){
                T x=particles.X(i).x;
                ARRAY<int>& ar=(x>=0)?constrained_particles:externally_forced;
                x=abs(x);
                if(x>0.99) ar.Append(i);}
            for(int i=0;i<constrained_particles.m;i++) Add_Debug_Particle(particles.X(constrained_particles(i)),TV(1,0,0));
            for(int i=0;i<externally_forced.m;i++) Add_Debug_Particle(particles.X(externally_forced(i)),TV(0,1,0));
            initial_positions=particles.X;
            constrained_velocities.Resize(constrained_particles.m);
            scalar_curve.Add_Control_Point(0,0);
            scalar_curve.Add_Control_Point(9,4*pi);
            break;}
        case 64:{
            tests.Create_Mattress(mattress_grid,true,0,density);
            for(int i=0;i<particles.X.m;i++){
                T x=particles.X(i).x;
                ARRAY<int>& ar=(x>=0)?constrained_particles:externally_forced;
                x=abs(x);
                if(x>0.99) ar.Append(i);}
            for(int i=0;i<constrained_particles.m;i++) Add_Debug_Particle(particles.X(constrained_particles(i)),TV(1,0,0));
            for(int i=0;i<externally_forced.m;i++) Add_Debug_Particle(particles.X(externally_forced(i)),TV(0,1,0));
            initial_positions=particles.X;
            constrained_velocities.Resize(constrained_particles.m);
            scalar_curve.Add_Control_Point(0,0);
            scalar_curve.Add_Control_Point(1,2);
            break;}
        case 65:{
            tests.Create_Mattress(mattress_grid,true,0,density);
            for(int i=0;i<particles.X.m;i++){
                T x=particles.X(i).x;
                ARRAY<int>& ar=(x>=0)?constrained_particles:externally_forced;
                x=abs(x);
                if(x>0.99) ar.Append(i);}
            for(int i=0;i<constrained_particles.m;i++) Add_Debug_Particle(particles.X(constrained_particles(i)),TV(1,0,0));
            for(int i=0;i<externally_forced.m;i++) Add_Debug_Particle(particles.X(externally_forced(i)),TV(0,1,0));
            initial_positions=particles.X;
            constrained_velocities.Resize(constrained_particles.m);
            scalar_curve.Add_Control_Point(0,0);
            scalar_curve.Add_Control_Point(3,3);
            scalar_curve.Add_Control_Point(3.1,2.5);
            scalar_curve.Add_Control_Point(3.2,3);
            break;}
        case 66:{
//            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_115K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)2.5,0)*m,ROTATION<TV>(pi/2,TV(1,0,0)))),true,true,density,m);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/adaptive_torus_float.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)2.5,0)*m,ROTATION<TV>(pi/2,TV(1,0,0)))),true,true,density,m);
            tests.Add_Ground();
            T sphere_radius=0.75*m;
            T torus_radius=1.6*m;
            T shell_outer_radius=1.1*m;
            T shell_height=1*m;
            T spacing=0.1*m;

            RIGID_BODY<TV>& torus=tests.Add_Analytic_Torus(torus_radius-shell_outer_radius-spacing,torus_radius,32,64);
            torus.Frame()=FRAME<TV>(TV(0,(T)2,0)*m,ROTATION<TV>(pi/2,TV(1,0,0)));
            torus.is_static=true;

            RIGID_BODY<TV>& sphere=tests.Add_Rigid_Body("sphere",sphere_radius,(T)0.3);
            sphere.Frame().t=TV(0,(T)2,0)*m;
            sphere.is_static=true;

            RIGID_BODY<TV>& shell=tests.Add_Analytic_Shell(shell_height,shell_outer_radius,sphere_radius+spacing,32);
            shell.is_static=false;

            kinematic_ids.Append(shell.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(shell.particle_index)=true;

            T y_start=3.75*m;
            T y_stop=0.75*m;
            T t_stop=2;

            curves.Resize(1);
            curves(0).Add_Control_Point(0,FRAME<TV>(TV(0,y_start,0)*m,ROTATION<TV>(pi/2,TV(1,0,0))));
            curves(0).Add_Control_Point(t_stop,FRAME<TV>(TV(0,y_stop,0)*m,ROTATION<TV>(pi/2,TV(1,0,0))));
            break;}
        case 67:{
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,1.01,0)*m));
            tests.Create_Mattress(mattress_grid,true,&initial_state,density);
            RIGID_BODY_STATE<TV> initial_state2(FRAME<TV>(TV(1.1,3.1,1.1)*m));
            tests.Create_Mattress(mattress_grid,true,&initial_state2,density);
            tests.Add_Ground();
            break;}
        case 68:{
            for(int i=0;i<5;i++){
              RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,1.001+2.001*i,0)*m));
              tests.Create_Mattress(mattress_grid,true,&initial_state,density);
            }
            RIGID_BODY<TV>& spike=tests.Add_Analytic_Box(TV(2.5,2.5,2.5)*m);
            ROTATION<TV> r(MATRIX<T,3>(sqrt(2)/2,sqrt(6)/6,sqrt(3)/3,-sqrt(2)/2,sqrt(6)/6,sqrt(3)/3,0,-sqrt(6)/3,sqrt(3)/3));
            TV start(0,collision_height*m,collision_speed*3*m);
            TV end(0,collision_height*m,-collision_speed*30*m);
            spike.Frame().t=start;
            spike.Frame().r=r;
            spike.is_static=false;
            kinematic_ids.Append(spike.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(spike.particle_index)=true;
            spike.coefficient_of_friction=input_friction;
            curves.Resize(1);
            curves(0).Add_Control_Point(0,FRAME<TV>(start,r));
            curves(0).Add_Control_Point(60,FRAME<TV>(end,r));
            tests.Add_Ground(input_friction);
            break;}
        case 69:{
            int spikes_per_side=4;
            T spike_scale=0.18;
            T spike_spacing=sqrt(3)*spike_scale;
            RIGID_BODY<TV>& box1=tests.Add_Analytic_Box(TV(0.0001,0.0001,0.0001)*m);
            RIGID_BODY<TV>& box2=tests.Add_Analytic_Box(TV(0.0001,0.0001,0.0001)*m);
            box1.Frame().t=TV(0,0.5,-1)*m;
            box2.Frame().t=TV(0,0.5,1)*m;
            box1.is_static=false;
            box2.is_static=false;
            kinematic_ids.Append(box1.particle_index);
            kinematic_ids.Append(box2.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(box1.particle_index)=true;
            rigid_body_collection.rigid_body_particles.kinematic(box2.particle_index)=true;
            box1.coefficient_of_friction=input_friction;
            box2.coefficient_of_friction=input_friction;
            curves.Resize(2+2*spikes_per_side*spikes_per_side);
            int c=0;
            T final_x = 1;
            T y = 0.78;
            T final_z = 0.55;
            curves(c).Add_Control_Point(0,FRAME<TV>(TV(0,y,-1)*m));
            curves(c).Add_Control_Point(1,FRAME<TV>(TV(0,y,-final_z)*m));
            curves(c).Add_Control_Point(2,FRAME<TV>(TV(0,y+spike_spacing,-final_z)*m));
            curves(c).Add_Control_Point(3,FRAME<TV>(TV(spike_spacing,y+spike_spacing,-final_z)*m));
            curves(c).Add_Control_Point(4,FRAME<TV>(TV(spike_spacing,y,-final_z)*m));
            curves(c).Add_Control_Point(8,FRAME<TV>(TV(final_x,y,-final_z)*m));
            c++;
            curves(c).Add_Control_Point(0,FRAME<TV>(TV(0,y,1)*m));
            curves(c).Add_Control_Point(1,FRAME<TV>(TV(0,y,final_z)*m));
            curves(c).Add_Control_Point(4,FRAME<TV>(TV(0,y,final_z)*m));
            curves(c).Add_Control_Point(8,FRAME<TV>(TV(-final_x,y,final_z)*m));
            c++;
            for(int i=0;i<spikes_per_side;i++){
              for(int j=0;j<spikes_per_side;j++){
                RIGID_BODY<TV>& spike=tests.Add_Analytic_Box(TV(1,1,1)*spike_scale*m);
                TV d((i-spikes_per_side/2)*spike_spacing*m,(j-spikes_per_side/2)*spike_spacing*m,(0.5-spike_spacing/6)*m);
                ROTATION<TV> r(MATRIX<T,3>(sqrt(2)/2,sqrt(6)/6,sqrt(3)/3,-sqrt(2)/2,sqrt(6)/6,sqrt(3)/3,0,-sqrt(6)/3,sqrt(3)/3));
                spike.Frame().t=box1.Frame().t+d;
                spike.Frame().r=r;
                spike.is_static=false;
                kinematic_ids.Append(spike.particle_index);
                rigid_body_collection.rigid_body_particles.kinematic(spike.particle_index)=true;
                spike.coefficient_of_friction=input_friction;
                const ARRAY<typename INTERPOLATION_CURVE<T,FRAME<TV> >::CONTROL_POINT>& cp=curves(0).control_points;
                for(int k=0;k<cp.m;k++)
                  curves(c).Add_Control_Point(cp(k).t,FRAME<TV>(cp(k).value.t+d,r));
                c++;}}
            for(int i=0;i<spikes_per_side;i++){
              for(int j=0;j<spikes_per_side;j++){
                RIGID_BODY<TV>& spike=tests.Add_Analytic_Box(TV(1,1,1)*spike_scale*m);
                TV d((i-spikes_per_side/2+0.5)*spike_spacing*m,(j-spikes_per_side/2+0.5)*spike_spacing*m,(-0.5+spike_spacing/6)*m);
                ROTATION<TV> r(MATRIX<T,3>(sqrt(2)/2,sqrt(6)/6,sqrt(3)/3,-sqrt(2)/2,sqrt(6)/6,sqrt(3)/3,0,-sqrt(6)/3,sqrt(3)/3));
                spike.Frame().t=box2.Frame().t+d;
                spike.Frame().r=r;
                spike.is_static=false;
                kinematic_ids.Append(spike.particle_index);
                rigid_body_collection.rigid_body_particles.kinematic(spike.particle_index)=true;
                spike.coefficient_of_friction=input_friction;
                const ARRAY<typename INTERPOLATION_CURVE<T,FRAME<TV> >::CONTROL_POINT>& cp=curves(1).control_points;
                for(int k=0;k<cp.m;k++)
                  curves(c).Add_Control_Point(cp(k).t,FRAME<TV>(cp(k).value.t+d,r));
                c++;}}
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_110K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T)0.522,(T)0)*m)),true,true,density,.007*m);
            tests.Add_Ground(10000);
            break;}
        case 77: {
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,0,0)*m));
            tests.Create_Mattress(mattress_grid,true,&initial_state,density);

            RIGID_BODY<TV>& box_bottom=tests.Add_Analytic_Box(TV(6,2,6)*m);
            RIGID_BODY<TV>& box_side_1=tests.Add_Analytic_Box(TV(2,6,6)*m);
            RIGID_BODY<TV>& box_side_2=tests.Add_Analytic_Box(TV(2,6,6)*m);
            RIGID_BODY<TV>& box_side_3=tests.Add_Analytic_Box(TV(6,6,2)*m);
            RIGID_BODY<TV>& box_side_4=tests.Add_Analytic_Box(TV(6,6,2)*m);
            RIGID_BODY<TV>& box_top=tests.Add_Analytic_Box(TV(6,2,6)*m);

            box_bottom.Frame().t=TV(0,-2,0)*m;
            box_side_1.Frame().t=TV(-2,0,0)*m;
            box_side_2.Frame().t=TV(2,0,0)*m;
            box_side_3.Frame().t=TV(0,0,-2)*m;
            box_side_4.Frame().t=TV(0,0,2)*m;

            box_bottom.is_static=true;
            box_side_1.is_static=true;
            box_side_2.is_static=true;
            box_side_3.is_static=true;
            box_side_4.is_static=true;
            box_top.is_static=false;

            box_bottom.coefficient_of_friction=0;
            box_side_1.coefficient_of_friction=0;
            box_side_2.coefficient_of_friction=0;
            box_side_3.coefficient_of_friction=0;
            box_side_4.coefficient_of_friction=0;
            box_top.coefficient_of_friction=0;

            kinematic_ids.Append(box_top.particle_index);
            rigid_body_collection.rigid_body_particles.kinematic(box_top.particle_index)=true;
            curves.Resize(1);
            curves(0).Add_Control_Point(0,FRAME<TV>(TV(0,2,0)));
            curves(0).Add_Control_Point(10,FRAME<TV>(TV(0,0,0)));
            if(!user_last_frame) last_frame=300;
            break;}
        case 80:{
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,4,0)*m));
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_20K.tet.gz", // adaptive_torus_float.tet
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(15,5,0)*m,ROTATION<TV>(-T(pi/2),TV(1,0,0)))),true,true,density,m);
            break;}
        case 100:{
            TETRAHEDRALIZED_VOLUME<T>* tv=TETRAHEDRALIZED_VOLUME<T>::Create(particles);
            particles.Add_Elements(4);
            particles.X(0)=TV(0,1,0);
/*            particles.X(1)=TV(1/sqrt(3.),0,0);
            particles.X(2)=TV(-1/sqrt(12),0,.5);
            particles.X(3)=TV(-1/sqrt(12),0,-.5);*/
            particles.X(1)=ROTATION<TV>(pi*1./12.,TV(0,1,0)).Rotate(TV(1/sqrt(3.),0,0));
            particles.X(2)=TV(particles.X(2).z,0,particles.X(2).x);
            particles.X(3)=TV(-1/sqrt(6.),0,-1/sqrt(6.)+1e-2);
            tv->mesh.elements.Append(VECTOR<int,4>(1,2,3,4));
            particles.mass.Fill(1);
            solid_body_collection.deformable_body_collection.Add_Structure(tv);
            contrail.Resize(1);
            break;}
        case 201:
        case 200:
        case 202:
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/liver_153k_tets_transformed.tet",RIGID_BODY_STATE<TV>(FRAME<TV>()),true,true,density,m);
            if(test_number==202){
                std::ifstream fin(data_directory+"/Tetrahedralized_Volumes/liver_153k_tets_transformed_real_hook_dirichlet_nodes.txt");
                int n;
                while(fin>>n) stuck_particles.Append(n);
                LOG::printf("stuck_particles %P\n",stuck_particles.m);}
            break;
        case 720:{
            GRID<TV> box_grid(TV_INT()+(resolution+1),RANGE<TV>::Centered_Box()*m);
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,3.02,0)*m));
            tests.Create_Mattress(box_grid,true,&initial_state,density);
            initial_state.frame.t=TV(.01,5.03,.02)*m;
            tests.Create_Mattress(box_grid,true,&initial_state,density);
            tests.Add_Analytic_Box(TV()+2,TV_INT()+resolution,density).Frame().t.y=1.01*m;
            tests.Add_Analytic_Box(TV()+2,TV_INT()+resolution,density).Frame().t.y=7.04*m;
            tests.Add_Analytic_Box(TV()+2,TV_INT()+resolution,density).Frame().t.y=9.05*m;
            tests.Add_Ground(0);
            break;}
        case 721:{
            GRID<TV> box_grid(TV_INT()+(resolution+1),RANGE<TV>::Centered_Box()*m);
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,1.1,0)*m));
            tests.Create_Mattress(box_grid,true,&initial_state,density);
            initial_state.frame.t.y=3.2*m;
            tests.Create_Mattress(box_grid,true,&initial_state,density);
            tests.Add_Ground(0);
            break;}
            // ./be_evolution 722 -last_frame 40 -no_self -rd_stiffness 5e6 -kry_tol 1e-12 -no_collisions_in_solve -newton_it 20 -dt .01 -rd_friction .1 -debug_newton -test_forces
        case 722:{
            int p=particles.Add_Element();
            particles.X(p)=TV(0,3,0);
            particles.mass.Fill(1);
//            tests.Add_Ground(0);
            Add_Collision_Object(new ANALYTIC_IMPLICIT_OBJECT<PLANE<T> >(PLANE<T>(TV(0,1,0),TV())));
            break;}
        case 701:{
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)5,0)*m)),true,true,density,m);
            tests.Add_Ground(0,1.99*m);
            break;}
        case 130:{
            tests.Create_Tetrahedralized_Volume(
                data_directory+"/Tetrahedralized_Volumes/adaptive_torus_float.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(4,0.5,4)*m,ROTATION<TV>::From_Euler_Angles((T)-0.1,(T)0,(T)0))),true,true,density,0.1*m);
            RIGID_BODY<TV>& torus=tests.Add_Analytic_Torus((T)0.05,(T)0.1,16,32,2*density);
            torus.Frame().t=TV(4,0.45,3.7);
            RIGID_BODY<TV>& ground=tests.Add_Analytic_Box(TV(8,0.1,8));
            ground.Frame().t=TV(4,0.05,4);
            ground.is_static=true;
            break;}
        case 140:{
            int number_side_panels=32;
            TRIANGULATED_SURFACE<T>& cloth=tests.Create_Cloth_Panel(number_side_panels,2,1,
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(),ROTATION<TV>::From_Euler_Angles(pi,(T)0,(T)0))));
            cloth.particles.X*=(T).5;
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(cloth,1);
            int n=number_side_panels+1;
            for(int i=0;i<n;i++)
                stuck_particles.Append(i);
            for(int i=0;i<n;i++)
                stuck_particles.Append(n*(n-1)+i);
            for(int i=1;i<n-1;i++){
                stuck_particles.Append(i*n);
                stuck_particles.Append(i*n+n-1);}
            tests.Create_Tetrahedralized_Volume(
                data_directory+"/Tetrahedralized_Volumes/torus.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(-0.3,0.16,0)*m,ROTATION<TV>())),true,true,density,0.4*m);
            tests.Create_Tetrahedralized_Volume(
                data_directory+"/Tetrahedralized_Volumes/torus.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0.3,0.16,0.2)*m,ROTATION<TV>())),true,true,density,0.4*m);
            RIGID_BODY<TV>& rt0=tests.Add_Analytic_Torus((T)0.05,(T)0.1,16,32,2*density);
            rt0.Frame().t=TV(0,0.16,0.3);
            RIGID_BODY<TV>& rt1=tests.Add_Analytic_Torus((T)0.05,(T)0.1,16,32,2*density);
            rt1.Frame().t=TV(0.2,0.16,0);
            RIGID_BODY<TV>& ground=tests.Add_Analytic_Box(TV(8,0.1,8));
            ground.Frame().t=TV(0,-(T)0.5,0);
            ground.is_static=true;
            break;}
        case 730:{
            T gap=0.005;
            T r=0.15;
            T R=r+gap;
            T k=0.8;
            TV center(0.5-R,r,0.5+tan(pi/6)*R);
            tests.Create_Tetrahedralized_Volume(
                data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(center)),true,true,density,r*m);
            center+=TV(2*R,0,0);
            tests.Create_Tetrahedralized_Volume(
                data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(center)),true,true,density,r*m);
            center+=TV(-R,0,-sin(pi/3)*2*R);
            tests.Create_Tetrahedralized_Volume(
                data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(center)),true,true,density,r*m);
            particles.X.template Project<T,&TV::y>()*=k;
            particles.X.template Project<T,&TV::y>()+=0.1;
            RIGID_BODY<TV>& sphere=tests.Add_Analytic_Sphere(r,density);
            sphere.Frame().t=TV(0.5,r*k+sqrt(6.0)/3*2*r+0.1-0.03,0.5);
            RIGID_BODY<TV>& ground=tests.Add_Analytic_Box(TV(1,0.1,1));
            ground.Frame().t=TV(0.5,0.05,0.5);
            ground.is_static=true;
            break;}
        case 740:{
            T angle=.1;
            T vel=.1;
            T angle_t=.2;
            ROTATION<TV> Q(angle,TV(0,0,1)),Qn(angle_t,TV(0,1,0));
            ROTATION<TV> R=Q*Qn;
            TV n,t,t1;
            Q.Get_Rotated_Frame(t,n,t1);
            TV c(.5,.5,.5);

            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(R.Rotate(TV(0,0.1,0))+c,R));
            initial_state.twist.linear=-t*vel;
            tests.Create_Mattress(mattress_grid,true,&initial_state,1);
            TV g=m/(s*s)*TV(0,-1.8,0);
            RIGID_BODY<TV>& rigid_body=tests.Add_Analytic_Box(TV(2,2,2));
            rigid_body.Frame().t=(c-TV(0,1e-4,0)-Q.Rotate(TV(0,1,0)))*m;
            rigid_body.Frame().r=Q;
            rigid_body.is_static=true;
            break;}
        case 750:case 751:{
            RIGID_BODY<TV>& lw=tests.Add_Analytic_Box(TV(0.4,8,8));
            lw.Frame().t=TV(-0.2,0.5,0.5);
            lw.is_static=true;
            RIGID_BODY<TV>& rw=tests.Add_Analytic_Box(TV(0.4,8,8));
            rw.Frame().t=TV(1.2,0.5,0.5);
            rw.is_static=true;

            TV g=m/(s*s)*TV(0,-1.8,0);
            RIGID_BODY<TV>& lcube=tests.Add_Analytic_Box(TV(0.52,0.5,0.5),TV_INT()+5,(T)1);
            lcube.Frame().t=TV(0.25,0.5,0.5)*m;
            tests.Create_Tetrahedralized_Volume(
                data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0.75,0.5,0.5))),true,true,(T)1,0.26*m);

            if(test_number==751){
                RIGID_BODY<TV>& sphere=tests.Add_Analytic_Sphere(.2,(T)10);
                sphere.Frame().t=TV(0.25,2,0.5)*m;}
            break;}
        default:
            LOG::cerr<<"Initial Data: Unrecognized test number "<<test_number<<std::endl;exit(1);}

    this->After_Get_Initial_Data(automatically_add_to_collision_structures);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    Get_Initial_Data();

    switch(test_number){
        case 1:{
                 use_penalty_self_collisions=false;}
               // Fallthrough
        case 2:
        case 3:
        case 5:
        case 6:
        case 8:
        case 13:
        case 16:
        case 66:
        case 80:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            if(test_number!=80) Add_Gravity();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e5*unit_p,(T).45,(T).01*s);
            if(test_number==80) particles.X.template Project<T,&TV::x>()*=-(T).97;
            if(test_number==80) particles.X.template Project<T,&TV::y>()*=(T).98;
            break;}
        case 67:
        case 9:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(0);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume2=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
            Add_Gravity();
            Add_Constitutive_Model(tetrahedralized_volume1,(T)1e5*unit_p,(T).45,(T).01*s);
            Add_Constitutive_Model(tetrahedralized_volume2,(T)1e5*unit_p,(T).45,(T).01*s);
            break;}
        case 10:{
            for(int i=0; i<tori_stack_height*sqr(tori_stack_width); i++){
              TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(i);
              Add_Constitutive_Model(tetrahedralized_volume,(T)5e5*unit_p,(T).45,(T).01*s);}
            deformable_body_collection.particles.V.Fill(TV(0,-attachment_velocity,0));
            Add_Gravity();
            break;}
        case 77:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e5*unit_p,(T).45,(T).01*s);
            break;}
        case 4: {
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Gravity();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e6*unit_p,(T).45,(T).01*s);
            break;}
        case 7:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Gravity();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e5*unit_p,(T).45,(T).01*s);
            for(int i=0;i<deformable_body_collection.particles.X.m;i++) deformable_body_collection.particles.X(i).y=3;
            break;}
        case 56:{break;}
        case 11:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)5e4*unit_p,(T).45,(T).01*s);
            Add_Gravity();
            break;}
        case 12:{
            if(!ether_drag) ether_drag=(T).6;
            deformable_body_collection.particles.V.Fill(TV(1,2,-3));
            break;}
        case 17:{
            use_penalty_self_collisions=false;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e5*unit_p,(T).45,(T).01*s);
            int s=resolution+1,s2=s*s,s3=s*s2;
            if(pin_corners){
                stuck_particles.Append(0);
                stuck_particles.Append(s-1);
                stuck_particles.Append(s2-1);
                stuck_particles.Append(s3-1);
                stuck_particles.Append(s2-s);
                stuck_particles.Append(s3-s2);
                stuck_particles.Append(s3-s);
                stuck_particles.Append(s3-s2+s-1);}
            ARRAY<TV> OX(particles.X.Subset(stuck_particles));
            rand.Fill_Uniform(particles.X,-1,1);
            particles.X.Subset(stuck_particles)=OX;
            break;}
        case 18:{
            use_penalty_self_collisions=false;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e5*unit_p,(T).45,(T).01*s);
            rand.Fill_Uniform(particles.X,-0,0);
            break;}
        case 23:
        case 24:
        case 25:
        case 26:
        case 27:
        case 57:
        case 28:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e5*unit_p,(T).45,(T).01*s);
            if(test_number==57){
                int m=mattress_grid.counts.x;
                int n=mattress_grid.counts.y;
                int p=mattress_grid.counts.z;
                for(int i=0;i<m;i++) for(int j=0;j<n;j++){constrained_particles.Append(i+m*j);constrained_particles.Append(i+m*j+(p-1)*m*n);}
                for(int i=0;i<m;i++) for(int k=0;k<p;k++){constrained_particles.Append(i+m*n*k);constrained_particles.Append(i+m*(n-1)+m*n*k);}
                constrained_velocities=particles.X.Subset(constrained_particles)*attachment_velocity;
                constrained_velocities.template Project<T,&TV::x>().Fill(0);}
            break;}
        case 29:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Gravity();
            Add_Constitutive_Model(tetrahedralized_volume,(T)0e2*unit_p,(T).45,(T).01*s);
            forces_are_removed=true;
            break;}

        case 30:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e6*unit_p,(T).40,(T).01*s);
            Add_Gravity();
            for(int i=0;i<deformable_body_collection.particles.X.m;i++){ deformable_body_collection.particles.V(i).x=(T)60;deformable_body_collection.particles.V(i).y=0;deformable_body_collection.particles.V(i).z=0;}
            break;}
        case 31:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e4*unit_p,(T).4,(T).001*s);

            for(int i=0;i<particles.X.m;i++){
                T y  =particles.X(i).y;
                TV v1=particles.X(i)-TV(0.27,0.645,-0.12);
                TV v2=particles.X(i)-TV(-0.275,0.605,-0.18);
                if ((y<=0.06) || (v1.Magnitude_Squared()<=sqr(0.065)) || (v2.Magnitude_Squared()<=sqr(0.065)))
                    constrained_particles.Append(i);}
            constrained_velocities.Resize(constrained_particles.m);
            break;}
        case 59:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e4*unit_p,(T).3,(T).001*s);

            ARRAY<TV> rand_particles(particles.X.m);
            rand.Fill_Uniform(rand_particles,-0.3,.3);

            for(int i=0;i<particles.X.m;i++){
                T y  =particles.X(i).y;
                TV v1=particles.X(i)-TV(0.27,0.645,-0.12);
                TV v2=particles.X(i)-TV(-0.275,0.605,-0.18);

                if ((y<=0.06) || (v1.Magnitude_Squared()<=sqr(0.065)) || (v2.Magnitude_Squared()<=sqr(0.065)))
                    constrained_particles.Append(i);
                else particles.X(i)=rand_particles(i)+TV(0,0.3,0)+TV(-0.024,0.06,0.1);}
            constrained_velocities.Resize(constrained_particles.m);

            break;}
        case 60:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e4*unit_p,(T).3,(T).001*s);
            rand.Fill_Uniform(particles.X,-0.3,.3);
            break;}
        case 32:{
            T youngs_modulus=7e5*unit_p;
            T poissons_ratio=.4;
            T damping=0.05*s;
            T g=0.8*m/(s*s);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(0);
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume2=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
            Add_Constitutive_Model(tetrahedralized_volume2,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume3=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(2);
            Add_Constitutive_Model(tetrahedralized_volume3,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume4=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(3);
            Add_Constitutive_Model(tetrahedralized_volume4,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume5=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(4);
            Add_Constitutive_Model(tetrahedralized_volume5,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume6=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(5);
            Add_Constitutive_Model(tetrahedralized_volume6,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume7=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(6);
            Add_Constitutive_Model(tetrahedralized_volume7,youngs_modulus,poissons_ratio,damping);
            Add_Gravity().gravity.y=-g;
            break;}
        case 43:
        case 33:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,1e4*unit_p,0.35,0.005*s);
            Add_Gravity().gravity=TV()*m/(s*s);
            break;}
        case 58:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(0);
            Add_Constitutive_Model(tetrahedralized_volume1,2e4*unit_p,0.4,0.005*s);
            if(!gears_of_pain){
                TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume2=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
                Add_Constitutive_Model(tetrahedralized_volume2,3e4*unit_p,0.4,0.005*s);
                //TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume3=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(2);
                //Add_Constitutive_Model(tetrahedralized_volume3,1e4*unit_p,0.4,0.005*s);
            }

            Add_Gravity();
            break;}
        case 34:{
            T youngs_modulus=1e5*unit_p;
            T poissons_ratio=.4;
            T damping=0.1*s;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            Add_Gravity();
            break;}
        case 35:{
            T youngs_modulus=2e5*unit_p;
            T poissons_ratio=.45;
            T damping=0.01*s;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(0);
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume2=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
            Add_Constitutive_Model(tetrahedralized_volume2,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume3=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(2);
            Add_Constitutive_Model(tetrahedralized_volume3,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume4=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(3);
            Add_Constitutive_Model(tetrahedralized_volume4,youngs_modulus,poissons_ratio,damping);
            Add_Gravity();
            break;}
        case 36:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e6*unit_p,(T).45,(T).01*s);
            //Add_Gravity();
            for(int i=0;i<deformable_body_collection.particles.X.m/2;i++){ deformable_body_collection.particles.V(i).x=(T)1;deformable_body_collection.particles.V(i).y=0;deformable_body_collection.particles.V(i).z=0;}
            break;}
        case 37:{
            T youngs_modulus=2.25e5*unit_p;
            T poissons_ratio=.4;
            T damping=0*s;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(0);
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume2=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
            Add_Constitutive_Model(tetrahedralized_volume2,youngs_modulus,poissons_ratio,damping);
            Add_Gravity();
            int m=mattress_grid.counts.x;
            int n=mattress_grid.counts.y;
            int mn=mattress_grid.counts.z;
            for(int i=0;i<m;i++)
            for(int j=0;j<n;j++)
            for(int k=0;k<mn;k++)
            {
                particles.V(i+m*j+m*n*k)=TV(-2*sin(j/(T)n),-cos(2*k/(T)mn)*2,-2*sin(3*i/(T)m));
                particles.V(i+m*j+m*n*k+m*n*mn)=TV(1*sin(2*k/(T)mn),0.5*cos(3*i/(T)m)*2,1*sin(j/(T)n));
            }
            break;}
        case 38:{
            T youngs_modulus=4e5*unit_p;
            T poissons_ratio=.4;
            T damping=0.001*s;
            int m=mattress_grid.counts.x;
            int n=mattress_grid.counts.y;
            int mn=mattress_grid.counts.z;
            for(int p=0;p<jello_centers.m;p++)
            {
                TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(p);
                Add_Constitutive_Model(tetrahedralized_volume,youngs_modulus,poissons_ratio,damping);
                for(int i=0;i<m;i++)
                for(int j=0;j<n;j++)
                for(int ij=0;ij<mn;ij++)
                {
                    particles.V(i+m*j+m*n*ij+(p-1)*m*n*mn)=TV(-3*sin(sin(127*p)*j/(T)n),-5*cos(sin(384*p)*4*ij/(T)mn)*2+5*sin(385*p),-6*sin(3*i*sin(457*p)/(T)m));
                }
            }
            Add_Gravity();
            break;}
        case 39:{
            T youngs_modulus=1e4*unit_p;
            T poissons_ratio=.4;
            T damping=0.001*s;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(0);
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume2=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
            Add_Constitutive_Model(tetrahedralized_volume2,youngs_modulus,poissons_ratio,damping);
            Add_Gravity();
            int m=mattress_grid.counts.x;
            int n=mattress_grid.counts.y;
            int mn=mattress_grid.counts.z;
            for(int i=0;i<m;i++)
                for(int j=0;j<n;j++)
                    for(int ij=0;ij<mn;ij++)
                    {
                        particles.V(i+m*j+m*n*ij)=TV(-0.3*sin(j/(T)n)-0.15,-cos(2*ij/(T)mn)*0.3,-0.15*sin(3*i/(T)m)-0.3);
                    }
            break;}
        case 40:{
            T youngs_modulus=1e4*unit_p;
            T poissons_ratio=.4;
            T damping=0.001*s;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(0);
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume2=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
            Add_Constitutive_Model(tetrahedralized_volume2,youngs_modulus,poissons_ratio,damping);
            Add_Gravity();
            int m=mattress_grid.counts.x;
            int n=mattress_grid.counts.y;
            int mn=mattress_grid.counts.z;
            for(int i=0;i<m;i++)
                for(int j=0;j<n;j++)
                    for(int ij=0;ij<mn;ij++)
                    {
                        particles.V(i+m*j+m*n*ij)=TV(0.32*cos(6*(ij+1)/(T)mn)+2.15,-0.28*cos(4*(i+1)/(T)m)+0.55,-0.31*sin(4*(j+1)/(T)n));
                        particles.V(i+m*j+m*n*ij+m*n*mn)=TV(0.27*sin(6*(j+1)/(T)n)-2.17,0.32*cos(6*(ij+1)/(T)mn)+0.53,0.25*sin(4*(i+1)/(T)m));
                    }
            for(int i=0;i<m*n*mn;i++)
            {
                int index=i;
                particles.V(index)+=TV(particles.X(index).y-jello_centers(0).y,-(particles.X(index).x-jello_centers(0).x),0)*51;
                particles.V(index)+=TV(0,particles.X(index).z-jello_centers(0).z,-(particles.X(index).y-jello_centers(0).y))*(-13);
                index+=m*n*mn;
                particles.V(index)+=TV(particles.X(index).y-jello_centers(1).y,-(particles.X(index).x-jello_centers(1).x),0)*(-39);
                particles.V(index)+=TV(0,particles.X(index).z-jello_centers(1).z,-(particles.X(index).y-jello_centers(1).y))*23;
            }
            break;}
        case 41:{
            T youngs_modulus=1e3*unit_p;
            T poissons_ratio=.4;
            T damping=0.001*s;

            for(int k=0;k<number_of_jellos;k++){
                TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(k);
                Add_Constitutive_Model(tetrahedralized_volume,youngs_modulus,poissons_ratio,damping);}
            Add_Gravity();
            break;
        }
        case 42:{
            T youngs_modulus=3e5*unit_p;
            T poissons_ratio=.4;
            T damping=0.001*s;
            int m=mattress_grid.counts.x;
            int n=mattress_grid.counts.y;
            int mn=mattress_grid.counts.z;
            for(int p=0;p<jello_centers.m;p++)
            {
                TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(p);
                Add_Constitutive_Model(tetrahedralized_volume,youngs_modulus,poissons_ratio,damping);
                for(int i=0;i<m;i++)
                for(int j=0;j<n;j++)
                for(int ij=0;ij<mn;ij++){
                    particles.V(i+m*j+m*n*ij+(p-1)*m*n*mn)=TV(-2*sin(sin(127*p)*j/(T)n)+17+4*sin(181*p),-3*cos(sin(384*p)*4*ij/(T)mn)*2+(2+4*sin(461*p)),1.5*sin(3*i*sin(457*p)/(T)m));}
                for(int i=0;i<m*n*mn;i++){
                    int index=i+(p-1)*m*n*mn;
                    particles.V(index)+=TV(particles.X(index).y-jello_centers(p).y,-(particles.X(index).x-jello_centers(p).x),0)*(6+2*sin(453*p));
                    particles.V(index)+=TV(0,particles.X(index).z-jello_centers(p).z,-(particles.X(index).y-jello_centers(p).y))*(cos(413*p));}
            }
            Add_Gravity();
            break;}
        case 52:{
            T youngs_modulus=1e4*unit_p;
            T poissons_ratio=.4;
            T damping=0.001*s;
            for(int k=0;k<jello_centers.m;k++)
            {
                TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(k);
                Add_Constitutive_Model(tetrahedralized_volume,youngs_modulus,poissons_ratio,damping);
            }
            Add_Gravity();
            break;}
        case 44:{
            T youngs_modulus=4e5*unit_p;
            T poissons_ratio=.4;
            T damping=0.001*s;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(0);
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume2=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
            Add_Constitutive_Model(tetrahedralized_volume2,youngs_modulus,poissons_ratio,damping);
            Add_Gravity();
            int m=mattress_grid.counts.x;
            int n=mattress_grid.counts.y;
            int mn=mattress_grid.counts.z;
            for(int i=0;i<m;i++)
                for(int j=0;j<n;j++)
                    for(int ij=0;ij<mn;ij++)
                    {
                        particles.V(i+m*j+m*n*ij)=TV(-2*sin(j/(T)n)+19,-cos(2*ij/(T)mn)*2+1,-2*sin(3*i/(T)m));
                        particles.V(i+m*j+m*n*ij+m*n*mn)=TV(1*sin(2*ij/(T)mn)-21,0.5*cos(3*i/(T)m)*2+1.2,1*sin(j/(T)n));
                    }
            for(int i=0;i<m*n*mn;i++)
            {
                int index=i;
                particles.V(index)+=TV(particles.X(index).y-jello_centers(0).y,-(particles.X(index).x-jello_centers(0).x),0)*9;
                particles.V(index)+=TV(0,particles.X(index).z-jello_centers(0).z,-(particles.X(index).y-jello_centers(0).y))*2;
                index+=m*n*mn;
                particles.V(index)+=TV(particles.X(index).y-jello_centers(1).y,-(particles.X(index).x-jello_centers(1).x),0)*(-12);
                particles.V(index)+=TV(0,particles.X(index).z-jello_centers(1).z,-(particles.X(index).y-jello_centers(1).y))*(-1);
            }
            break;}
        case 50:{
            T youngs_modulus=1e5*unit_p;
            T poissons_ratio=.45;
            T damping=0.01*s;
           // T g=0.8;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            for(int i=0;i<particles.X.m;i++)
                if(particles.X(i).x>=1.5)
                    externally_forced.Append(i);

            break;}
        case 51:{
            T youngs_modulus=1e5*unit_p;
            T poissons_ratio=.45;
            T damping=0.01*s;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            for(int i=0;i<particles.X.m;i++)
                if(particles.X(i).x>=1.5)
                    externally_forced.Append(i);
            break;}
        case 55:
        case 54:
        case 48:
        case 53:{
            T youngs_modulus=1e5*unit_p;
            T poissons_ratio=.45;
            T damping=0.01*s;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            if(test_number==55) particles.X(1).x=stretch;
            break;}
        case 61:
        case 62:
        case 63:
        case 64:
        case 65:{
            T youngs_modulus=1e5*unit_p;
            T poissons_ratio=.45;
            T damping=0.01*s;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            break;}
        case 68:{
            for(int i=0; i<5; i++){
              TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(i);
              T youngs_modulus = 1e5*unit_p*(i%2==0?10:1);
              Add_Constitutive_Model(tetrahedralized_volume,youngs_modulus,(T).45,(T).01*s);}
            Add_Gravity();
            break;}
        case 69:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e6,(T).45,(T).01*s);
            Add_Gravity();
            break;}
        case 100:{
            T youngs_modulus=1e5*unit_p;
            T poissons_ratio=.4;
            T damping=0.1*s;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            break;}
        case 200:
        case 202:
        case 201:{
            T youngs_modulus=1e5*unit_p;
            T poissons_ratio=.3;
            T damping=0;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            if(test_number==200) particles.X*=1.1;
            if(test_number==201) rand.Fill_Uniform(particles.X,-1,1);
            if(test_number==202) Add_Gravity();
            break;}
        case 701:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e5*unit_p,(T).45,(T).01*s);
            Add_Gravity();
            break;}
            // ./be_evolution 720 -last_frame 100 -no_self -rd_stiffness 1e5 -kry_tol 1e-4 -no_collisions_in_solve -newton_it 20 -dt .005 -resolution 5
            // ./be_evolution 721 -last_frame 100 -no_self -rd_stiffness 1e5 -kry_tol 1e-4 -no_collisions_in_solve -newton_it 20 -dt .005 -resolution 3
        case 720:case 721:case 722:{
            for(int i=0;;i++){
                auto st=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>(i);
                if(!st) break;
                Add_Constitutive_Model(*st,(T)1e6*unit_p,(T).45,(T).01*s);}
            Add_Gravity();
            break;}
        case 730:
        case 740:
        case 750:case 751:{
            for(int i=0;;i++){
                auto st=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>(i);
                if(!st) break;
                Add_Constitutive_Model(*st,(T)1e6*unit_p,(T).45,(T).01*s);}
            TV gravity=m/(s*s)*TV(0,-1.8,0);
            DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
            GRAVITY<TV>* g=new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true,gravity);
            solid_body_collection.Add_Force(g);
            break;}
        case 130:{
            for(int i=0;;i++){
                auto st=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>(i);
                if(!st) break;
                Add_Constitutive_Model(*st,(T)1e4*unit_p,(T).45,(T).01*s);}
            TV gravity=m/(s*s)*TV(0,-1.8,0);
            DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
            GRAVITY<TV>* g=new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true,gravity);
            solid_body_collection.Add_Force(g);
            break;}
        case 140:{
            TRIANGULATED_SURFACE<T>& cloth=deformable_body_collection.template Find_Structure<TRIANGULATED_SURFACE<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,
                solid_body_collection.rigid_body_collection,true,true));
            T linear_stiffness=100*stiffness_multiplier*10/(1+sqrt((T)2)),linear_damping=damping_multiplier*15;
            solid_body_collection.Add_Force(Create_Edge_Springs(cloth,linear_stiffness,linear_damping));
            T bending_stiffness_multiplier=100,bending_damping_multiplier=1;
            T bending_stiffness=bending_stiffness_multiplier*2/(1+sqrt((T)2)),bending_damping=bending_damping_multiplier*8;
            solid_body_collection.Add_Force(Create_Bending_Springs(cloth,bending_stiffness,bending_damping));
            for(int i=0;;i++){
                auto st=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>(i);
                if(!st) break;
                Add_Constitutive_Model(*st,(T)1e4*unit_p,(T).45,(T)0.01*s);}
            break;}
        default:
            LOG::cerr<<"Missing bodies implementation for test number "<<test_number<<std::endl;exit(1);}

    if(ether_drag) solid_body_collection.Add_Force(new ELASTIC_ETHER_DRAG<TV>(deformable_body_collection.particles,true,ether_drag,1,save_dt));

    if(solid_body_collection.deformable_body_collection.mpi_solids){
        VECTOR<int,3> processes_per_dimension(2,1,1);
        deformable_body_collection.mpi_solids->Simple_Partition(deformable_body_collection,solid_body_collection.rigid_body_collection,particles.X,processes_per_dimension);}

    if(use_penalty_collisions)
        for(int b=0;b<rigid_body_collection.rigid_body_particles.number;b++){
            IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> > *iot=rigid_body_collection.Rigid_Body(b).implicit_object;
            if(LEVELSET_IMPLICIT_OBJECT<TV>* lio=dynamic_cast<LEVELSET_IMPLICIT_OBJECT<TV>*>(iot->object_space_implicit_object))
                iot=new IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >(new SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV>(lio->levelset.grid,lio->levelset.phi),true,iot->transform);
            solid_body_collection.Add_Force(new IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>(particles,
                    iot,penalty_collisions_stiffness,penalty_collisions_separation,penalty_collisions_length));}
    if(use_constraint_collisions && backward_euler_evolution)
        for(int b=0;b<rigid_body_collection.rigid_body_particles.number;b++){
            IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> > *iot=rigid_body_collection.Rigid_Body(b).implicit_object;
            if(LEVELSET_IMPLICIT_OBJECT<TV>* lio=dynamic_cast<LEVELSET_IMPLICIT_OBJECT<TV>*>(iot->object_space_implicit_object)){
                lio->levelset.interpolation=new CUBIC_SPLINE_INTERPOLATION_UNIFORM<TV,T,FACE_LOOKUP_UNIFORM<TV> >();
                iot=new IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >(new SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV>(lio->levelset.grid,lio->levelset.phi),true,iot->transform);}
            backward_euler_evolution->minimization_objective.collision_objects.Append(iot);
            backward_euler_evolution->minimization_objective.coefficient_of_friction.Append(rigid_body_collection.Rigid_Body(b).coefficient_of_friction);}
    else if(!use_penalty_collisions)
        for(int i=0;i<deformable_body_collection.structures.m;i++){
            deformable_body_collection.collisions.collision_structures.Append(deformable_body_collection.structures(i));
            if(solids_parameters.triangle_collision_parameters.perform_self_collision)
                solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(deformable_body_collection.structures(i));}

    if(use_penalty_self_collisions){
        if(use_distance_based_self_collisions){
            for(int b=0;b<deformable_body_collection.structures.m;b++){
                DEFORMABLE_PARTICLES<TV>& undeformed_particles=*particles.Clone();
                TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(b);
                if(!tetrahedralized_volume.triangulated_surface) tetrahedralized_volume.Initialize_Triangulated_Surface();
                TRIANGULATED_SURFACE<T>* triangulated_surface=tetrahedralized_volume.triangulated_surface;
                TRIANGULATED_SURFACE<T>& undeformed_triangulated_surface=*(new TRIANGULATED_SURFACE<T>(triangulated_surface->mesh,undeformed_particles));
                LEVELSET_IMPLICIT_OBJECT<TV>& undeformed_levelset=*Initialize_Implicit_Surface(undeformed_triangulated_surface,10);
                DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>* coll=new DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>(particles,undeformed_particles,tetrahedralized_volume,
                    undeformed_triangulated_surface,undeformed_levelset,penalty_collisions_stiffness,penalty_collisions_separation);
                if(self_collide_surface_only){
                    coll->colliding_particles=tetrahedralized_volume.triangulated_surface->mesh.elements.Flattened();
                    Prune_Duplicates(coll->colliding_particles);}
                int force_id=solid_body_collection.Add_Force(coll);
                if(backward_euler_evolution) backward_euler_evolution->minimization_objective.deformables_forces_lazy.Set(force_id);}}
        else{
                LEVELSET_VOLUME_COLLISIONS<TV>* lvc=new LEVELSET_VOLUME_COLLISIONS<TV>(particles,penalty_collisions_stiffness);
                for(int b=0;b<deformable_body_collection.structures.m;b++){
                    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(b);
                    LEVELSET_IMPLICIT_OBJECT<TV>& undeformed_levelset=*Initialize_Implicit_Surface(tetrahedralized_volume.Get_Boundary_Object(),20);
                    lvc->Add_Mesh(tetrahedralized_volume,undeformed_levelset);
                    delete &undeformed_levelset;}
                int force_id=solid_body_collection.Add_Force(lvc);
                if(backward_euler_evolution)
                    backward_euler_evolution->minimization_objective.deformables_forces_lazy.Set(force_id);}}

    if(solids_parameters.triangle_collision_parameters.perform_self_collision){
        solid_body_collection.Add_Force(new TRIANGLE_REPULSIONS_PENALTY<TV>(particles,deformable_body_collection.triangle_repulsions.point_face_interaction_pairs));
        solid_body_collection.Add_Force(new TRIANGLE_REPULSIONS_PENALTY<TV>(particles,deformable_body_collection.triangle_repulsions.edge_edge_interaction_pairs));}

    if(enforce_definiteness) solid_body_collection.Enforce_Definiteness(true);
    for(int i=0;i<deformable_body_collection.deformables_forces.m;i++)
        if(COLLISION_FORCE<TV>* cf=dynamic_cast<COLLISION_FORCE<TV>*>(solid_body_collection.deformable_body_collection.deformables_forces(i)))
            cf->coefficient_of_friction=input_friction;

    this->After_Initialize_Bodies();
}

void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) override {}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    V.Subset(constrained_particles)=constrained_velocities;

    T final_time=50;
    if(test_number==24){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        TV velocity_x=velocity_time<final_time?TV(attachment_velocity,0,0):TV();
        TV velocity_y=velocity_time<final_time?TV(0,attachment_velocity,0):TV();
        TV velocity_z=velocity_time<final_time?TV(0,0,attachment_velocity):TV();
        for(int i=m/3+0;i<2*m/3;i++)for(int j=n/3;j<2*n/3;j++){V(i+m*j)=-velocity_z;V(i+m*j+(mn-1)*m*n)=velocity_z;}
        for(int i=m/3+0;i<2*m/3;i++)for(int ij=mn/3;ij<2*mn/3;ij++){V(i+m*n*ij)=-velocity_y;V(i+m*(n-1)+m*n*ij)=velocity_y;}
        for(int ij=mn/3;ij<2*mn/3;ij++)for(int j=n/3;j<2*n/3;j++){V(m*j+m*n*ij)=-velocity_x;V(m-1+m*j+m*n*ij)=velocity_x;}}
    if(test_number==25){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        V(0)=velocity_time<final_time?TV(-attachment_velocity,-attachment_velocity,-attachment_velocity):TV();
        V(m-1)=velocity_time<final_time?TV( attachment_velocity,-attachment_velocity,-attachment_velocity):TV();
        V(m*(n-1))=velocity_time<final_time?TV(-attachment_velocity, attachment_velocity,-attachment_velocity):TV();
        V(m*n-1)=velocity_time<final_time?TV( attachment_velocity, attachment_velocity,-attachment_velocity):TV();
        V((mn-1)*n*m)=velocity_time<final_time?TV(-attachment_velocity,-attachment_velocity,attachment_velocity):TV();
        V((mn-1)*n*m+m-1)=velocity_time<final_time?TV( attachment_velocity,-attachment_velocity,attachment_velocity):TV();
        V((mn-1)*n*m+m*(n-1))=velocity_time<final_time?TV(-attachment_velocity, attachment_velocity,attachment_velocity):TV();
        V(mn*n*m-1)=velocity_time<final_time?TV( attachment_velocity, attachment_velocity,attachment_velocity):TV();}
    if(test_number==26){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        TV velocity_x=velocity_time<final_time?TV(attachment_velocity,0,0):TV();
        TV velocity_y=velocity_time<final_time?TV(0,attachment_velocity,0):TV();
        for(int ij=0;ij<mn;ij++)for(int j=0;j<n;j++){V(m*j+m*n*ij)=-velocity_x;V(m-1+m*j+m*n*ij)=velocity_x;}
        for(int i=3*m/7+0;i<4*m/7;i++)for(int ij=0;ij<mn;ij++)for(int j=0;j<n;j++){V(i+m*j+m*n*ij)=-velocity_y;V(i+m*j+m*n*ij)=-velocity_y;}}
    if(test_number==27){
        final_time=70;
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        TV velocity_x=velocity_time<final_time?TV(attachment_velocity,0,0):TV();
        for(int ij=0;ij<mn;ij++)for(int j=0;j<n;j++){V(m*j+m*n*ij)=velocity_x;V(m-1+m*j+m*n*ij)=-velocity_x;}}
    if(test_number==23){
        final_time=12;
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        TV velocity_x=velocity_time<final_time?TV(attachment_velocity,0,0):TV();
        for(int ij=0;ij<mn;ij++)for(int j=0;j<n;j++){V(m*j+m*n*ij)=-(T)(velocity_time<=25)*velocity_x;V(m-1+m*j+m*n*ij)=(T)(velocity_time<=25)*velocity_x;}}
    if(test_number==28){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        TV velocity_x=velocity_time<final_time?TV(attachment_velocity,0,0):TV();
        TV velocity_y=velocity_time<final_time?TV(0,attachment_velocity,0):TV();
        for(int ij=0;ij<mn;ij++)for(int j=0;j<n;j++){V(m*j+m*n*ij)=-velocity_x;V(m-1+m*j+m*n*ij)=velocity_x;}}
    if(test_number==53){V(0)=TV(1,0,0);V(1).x=0;V(2).x=0;V(3).x=0;}
    if(test_number==54){V(0)=TV(1,0,0);V(1)=TV();V(2).x=0;V(3).x=0;}
    if(test_number==55){V(1)=V(0)=TV();V(2).x=0;V(3).x=0;}
    if(test_number==100){V(0)=TV(0,(particles.X(0).y<5-1e-4)*.5,0);V(1).y=0;V(2).y=0;V(3).y=0;}
    if(test_number==48){
        for(int i=0;i<stuck_particles.m;i++){
            int p=stuck_particles(i);
            V(p)=V(p).Projected_Orthogonal_To_Unit_Direction(particles.X(p).Normalized());}}
    if(test_number==17) V.Subset(stuck_particles).Fill(TV());
    if(test_number==61){
        TV axis(0,0,1),d_angle_axis=scalar_curve.Derivative(velocity_time)*axis;
        ROTATION<TV> rot(scalar_curve.Value(velocity_time),axis);
        for(int i=0;i<externally_forced.m;i++) V(externally_forced(i))=rot.Rotate(d_angle_axis.Cross(initial_positions(externally_forced(i))));}
    if(test_number==62){
        TV axis(1,0,0),d_angle_axis=scalar_curve.Derivative(velocity_time)*axis;
        ROTATION<TV> rot(scalar_curve.Value(velocity_time),axis);
        for(int i=0;i<externally_forced.m;i++) V(externally_forced(i))=rot.Rotate(d_angle_axis.Cross(initial_positions(externally_forced(i))));}
    if(test_number==63){
        TV axis(1,0,0),d_angle_axis=scalar_curve.Derivative(velocity_time)*axis;
        ROTATION<TV> rot(scalar_curve.Value(velocity_time),axis);
        for(int i=0;i<externally_forced.m;i++){
          V(externally_forced(i))=rot.Rotate(d_angle_axis.Cross(initial_positions(externally_forced(i))));
          if(velocity_time<9)
            V(externally_forced(i))-=(T)5*axis;}}
    if(test_number==64||test_number==65){
        TV axis(-1,0,0);
        for(int i=0;i<externally_forced.m;i++)
          V(externally_forced(i))=scalar_curve.Derivative(velocity_time)*axis;}
}
void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > rotation,const T time) override {}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) override
{
    if(test_number==61){
        ROTATION<TV> rot(scalar_curve.Value(time),TV(0,0,1));
        for(int i=0;i<externally_forced.m;i++) X(externally_forced(i))=rot.Rotate(initial_positions(externally_forced(i)));}
    if(test_number==62){
        ROTATION<TV> rot(scalar_curve.Value(time),TV(1,0,0));
        for(int i=0;i<externally_forced.m;i++) X(externally_forced(i))=rot.Rotate(initial_positions(externally_forced(i)));}
    if(test_number==63){
        ROTATION<TV> rot(scalar_curve.Value(time),TV(1,0,0));
        for(int i=0;i<externally_forced.m;i++){
          X(externally_forced(i))=rot.Rotate(initial_positions(externally_forced(i)));
          if(time<9)
            X(externally_forced(i))-=TV(5*time,0,0);
          else
            X(externally_forced(i))-=TV(45,0,0);}}
    if(test_number==64||test_number==65)
        for(int i=0;i<externally_forced.m;i++) X(externally_forced(i))=scalar_curve.Value(time)*TV(-1,0,0)+initial_positions(externally_forced(i));
}
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) override {}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    V.Subset(constrained_particles).Fill(TV());
    V.Subset(externally_forced).Fill(TV());

    if(test_number==24){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        for(int i=m/3+0;i<2*m/3;i++)for(int j=n/3;j<2*n/3;j++){V(i+m*j)=TV();V(i+m*j+(mn-1)*m*n)=TV();}
        for(int i=m/3+0;i<2*m/3;i++)for(int ij=mn/3;ij<2*mn/3;ij++){V(i+m*n*ij)=TV();V(i+m*(n-1)+m*n*ij)=TV();}
        for(int ij=mn/3;ij<2*mn/3;ij++)for(int j=n/3;j<2*n/3;j++){V(m*j+m*n*ij)=TV();V(m-1+m*j+m*n*ij)=TV();}}
    if(test_number==25){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        V(0)=TV();
        V(m-1)=TV();
        V(m*(n-1))=TV();
        V(m*n-1)=TV();
        V((mn-1)*n*m)=TV();
        V((mn-1)*n*m+m-1)=TV();
        V((mn-1)*n*m+m*(n-1))=TV();
        V(mn*n*m-1)=TV();}
    if(test_number==26){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        for(int ij=0;ij<mn;ij++)for(int j=0;j<n;j++){V(m*j+m*n*ij)=TV();V(m-1+m*j+m*n*ij)=TV();}
        for(int i=3*m/7+0;i<4*m/7;i++)for(int ij=0;ij<mn;ij++)for(int j=0;j<n;j++){V(i+m*j+m*n*ij)=TV();V(i+m*j+m*n*ij)=TV();}}
    if(test_number==23){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        for(int ij=0;ij<mn;ij++)for(int j=0;j<n;j++){V(m*j+m*n*ij)=TV();V(m-1+m*j+m*n*ij)=TV();}}
    if(test_number==27){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        for(int ij=0;ij<mn;ij++)for(int j=0;j<n;j++){V(m*j+m*n*ij)=TV();V(m-1+m*j+m*n*ij)=TV();}}
    if(test_number==28){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        for(int ij=0;ij<mn;ij++)for(int j=0;j<n;j++){V(m*j+m*n*ij)=TV();V(m-1+m*j+m*n*ij)=TV();}}
    if(test_number==53){V(0)=TV();V(1).x=0;V(2).x=0;V(3).x=0;}
    if(test_number==54 || test_number==55){V(0)=V(1)=TV();V(2).x=0;V(3).x=0;}
    if(test_number==24){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        for(int i=m/3+0;i<2*m/3;i++)for(int j=n/3;j<2*n/3;j++){V(i+m*j)=TV();V(i+m*j+(mn-1)*m*n)=TV();}
        for(int i=m/3+0;i<2*m/3;i++)for(int ij=mn/3;ij<2*mn/3;ij++){V(i+m*n*ij)=TV();V(i+m*(n-1)+m*n*ij)=TV();}}
    if(test_number==100){V(0)=TV();V(1).y=0;V(2).y=0;V(3).y=0;}
    if(test_number==48){
        for(int i=0;i<stuck_particles.m;i++){
            int p=stuck_particles(i);
            V(p)=V(p).Projected_Orthogonal_To_Unit_Direction(particles.X(p).Normalized());}}
    else V.Subset(stuck_particles).Fill(TV());
}
//#####################################################################
// Function Read_Output_Files_Solids
//#####################################################################
void Read_Output_Files_Solids() override
{
    BASE::Read_Output_Files_Solids();
    solid_body_collection.Update_Simulated_Particles();
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) override
{
    if(test_number==41)
    {
        if(id==kinematic_ids(0))
        {if(time>=.2) frame=FRAME<TV>(TV(1,plateau,0));else frame=curves(0).Value(time);}
        if(id==kinematic_ids(1))
        {if(time>=.27) frame=FRAME<TV>(TV(2,plateau,0));else frame=curves(1).Value(time);}
        if(id==kinematic_ids(2))
        {if(time>=.33) frame=FRAME<TV>(TV(3,plateau,0));else frame=curves(2).Value(time);}
        if(id==kinematic_ids(3))
        {if(time>=.4) frame=FRAME<TV>(TV(4,plateau,0));else frame=curves(3).Value(time);}
        if(id==kinematic_ids(4))
        {if(time>=.47) frame=FRAME<TV>(TV(5,plateau,0));else frame=curves(4).Value(time);}
        if(id==kinematic_ids(5))
        {if(time>=.53) frame=FRAME<TV>(TV(6,plateau,0));else frame=curves(5).Value(time);}
        if(id==kinematic_ids(6))
        {if(time>=.6) frame=FRAME<TV>(TV(7,plateau,0));else frame=curves(6).Value(time);}
        if(id==kinematic_ids(7))
        {if(time>=.67) frame=FRAME<TV>(TV(8,plateau,0));else frame=curves(7).Value(time);}
        return;}
    for(int i=0;i<kinematic_ids.m;i++)
        if(id==kinematic_ids(i)){
            frame=curves(i).Value(time);
            break;}
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) override
{
    for(int i=0;i<kinematic_ids.m;i++)
        if(id==kinematic_ids(i)){
            twist=curves(i).Derivative(time);
            return true;}
    return false;
}
//#####################################################################
// Function Preprocess_Substep
//#####################################################################
void Preprocess_Substep(const T dt,const T time) override
{
    BASE::Preprocess_Substep(dt,time);
    save_dt=dt;
    if(test_number==11) solid_body_collection.template Find_Force<GRAVITY<TV>&>().gravity=TV(0,-10*time*m/(s*s*s),0);
    if(test_number==33)
    {
        TETRAHEDRALIZED_VOLUME<T>& tet_volume=solid_body_collection.deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
        FINITE_VOLUME<TV,3>& fvm=solid_body_collection.deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,3>&>();

        int number_of_vertices=solid_body_collection.deformable_body_collection.collisions.check_collision.m;
        for(int i=0;i<number_of_vertices;i++) solid_body_collection.deformable_body_collection.collisions.check_collision(i)=false;

        for(int t=0;t<fvm.Fe_hat.m;t++)
            if(fvm.Fe_hat(t).x.x<3)
                for(int i=0;i<4;i++) solid_body_collection.deformable_body_collection.collisions.check_collision(tet_volume.mesh.elements(t)(i))=true;

        if (time >=1 ) solid_body_collection.template Find_Force<GRAVITY<TV>&>().gravity.y=-9.8*m/(s*s);}
    if(test_number==58)
    {
        solid_body_collection.deformable_body_collection.collisions.check_collision.Fill(true);
        for(int f=0;FINITE_VOLUME<TV,3>* fvm=solid_body_collection.deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,3>*>(f);f++)
            for(int t=0;t<fvm->Fe_hat.m;t++){
               // LOG::cout << "booya " << stretch_cutoff << std::endl;
                if(fvm->Fe_hat(t).x.x>=stretch_cutoff)
                    solid_body_collection.deformable_body_collection.collisions.check_collision.Subset(fvm->strain_measure.mesh_object.mesh.elements(t)).Fill(false);}}
    if(test_number==29) std::cout << "rame!" <<      solids_parameters.triangle_collision_parameters.perform_self_collision << std::endl;
    if(test_number==31) solid_body_collection.deformable_body_collection.collisions.check_collision.Subset(constrained_particles).Fill(false);
    if(test_number==51) fish_V=solid_body_collection.deformable_body_collection.particles.V;
}
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
void Update_Time_Varying_Material_Properties(const T time) override
{
    if(test_number==29 && time > .01){
        if(solids_parameters.triangle_collision_parameters.perform_self_collision && time<=1.3){
            solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=override_collisions;
            solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=override_collisions;
            solids_parameters.triangle_collision_parameters.perform_self_collision=override_collisions;}//This gets turned off later then back on

        T critical=(T)1.0;
        T critical2=(T)1.0+rebound_time;
        T critical3=(T)1.0+rebound_time+.3;

        T start_young=(T)0;T end_young=(T)rebound_stiffness;
        T pois=(T).45;if(input_poissons_ratio!=-1) pois=input_poissons_ratio;
        DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
        if(time<critical) forces_are_removed=true;
        if(time>critical && forces_are_removed) forces_are_removed=false;

        if(time>critical && time<critical2){
            FINITE_VOLUME<TV,3>& fv=deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,3>&>();
            CONSTITUTIVE_MODEL<T,3>& icm=fv.constitutive_model;
            T young=pow(10.0,start_young+(time-critical)/(critical2-critical)*(end_young-start_young));
            icm.Update_Lame_Constants(young,pois,(T).01);
            forces_are_removed=false;}
        if(time>critical && rebound_time < 1e-6){
            FINITE_VOLUME<TV,3>& fv=deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,3>&>();
            CONSTITUTIVE_MODEL<T,3>& icm=fv.constitutive_model;
            icm.Update_Lame_Constants(pow(10.0,end_young),pois,(T).01);
            forces_are_removed=false;}
        if(time>critical3 && self_collision_flipped==false){
            self_collision_flipped=true;
            FINITE_VOLUME<TV,3>& fv=deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,3>&>();
            CONSTITUTIVE_MODEL<T,3>& icm=fv.constitutive_model;
            T young=pow(10.0,end_young-rebound_drop);
            icm.Update_Lame_Constants(young,pois,(T).01);
            solids_parameters.triangle_collision_parameters.perform_self_collision=override_collisions;}}
}
//#####################################################################
// Function Postprocess_Substep
//#####################################################################
void Postprocess_Substep(const T dt,const T time) override
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    if(test_number==51) for(int i=0;i<particles.X.m;i++) if(particles.V(i).x>6) particles.V(i).x=6;
    if(test_number==48){
        stuck_particles.Remove_All();
        T r=std::max((T)0,hole-time*stretch);
        if(time>1.2*hole/stretch) r=std::min((time-(T)1.2*hole/stretch)*stretch,hole);
        for(int i=0;i<particles.X.m;i++)
            if(particles.X(i).Magnitude()>r){
                stuck_particles.Append(i);
                particles.X(i)*=r/particles.X(i).Magnitude();}}
    if(test_number==10){
            for(int i=0;i<tori_stack_height*sqr(tori_stack_width);i++){
              TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=solid_body_collection.deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(i);
              tetrahedralized_volume.Update_Bounding_Box();
              if(tetrahedralized_volume.bounding_box->min_corner.y>10*m){
                ARRAY<int> mesh_particles;
                Get_Unique(mesh_particles,tetrahedralized_volume.mesh.elements.Flattened());
                for(int i=0;i<mesh_particles.m;i++)
                  if(particles.V(mesh_particles(i)).y<-attachment_velocity)
                    particles.V(mesh_particles(i)).y=-attachment_velocity;}}}

    if(test_number==130 && time<4){
        RIGID_BODY<TV>& torus=solid_body_collection.rigid_body_collection.Rigid_Body(0);
        torus.Frame().t=TV(4,0.45,3.7);
        torus.Twist()=TWIST<TV>();}
    if(test_number==140 && time>2){
        int n=32+1;
        if(stuck_particles.m>n)
            stuck_particles.Resize(n);}

    BASE::Postprocess_Substep(dt,time);
}
//#####################################################################
// Function Bind_Intersecting_Particles
//#####################################################################
void Bind_Intersecting_Particles()
{
    if(nobind) return;
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    bool added_binding=false;
    for(int i=0;i<rigid_body_collection.static_and_kinematic_rigid_bodies.m;i++){int b=rigid_body_collection.static_and_kinematic_rigid_bodies(i);
        RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(b);
        for(int p=0;p<particles.X.m;p++){
            if(!binding_list.Binding(p) && rigid_body.Implicit_Geometry_Lazy_Inside(particles.X(p),1e-3)){
                added_binding=true;
                LOG::cout<<"adding binding: "<<b<<"  "<<p<<std::endl;
                binding_list.Add_Binding(new RIGID_BODY_BINDING<TV>(particles,p,rigid_body_collection,b,rigid_body.Object_Space_Point(particles.X(p))));}}}

    if(added_binding) solid_body_collection.Update_Simulated_Particles();
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) override
{
    if(test_number==32 && frame==1100) Bind_Intersecting_Particles();

    if(test_number==52){
        DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;

        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;

        for(int k=0;k<jello_centers.m;k++) if (frame==(k-1)*60+1){
            TV center=TV();
            for(int i=0;i<m*n*mn;i++){
                int index=i+(k-1)*m*n*mn;
                center+=particles.X(index);}
            center/=m*n*mn;
            for(int i=0;i<m*n*mn;i++){
                int index=i+(k-1)*m*n*mn;
                particles.X(index)+=TV(0,0.2,0)-center;}
            for(int i=0;i<m;i++)
            for(int j=0;j<n;j++)
            for(int ij=0;ij<mn;ij++){
                particles.V(i+m*j+m*n*ij+(k-1)*m*n*mn)=(T)0.1*TV(-sin(sin(187*k)*5*i/(T)m),-cos(cos(217*k)*6*j/(T)n),sin(5*ij*sin(471*k)/(T)mn));}

            for(int k_other=k+1;k_other<jello_centers.m;k_other++){
                TV center=TV();
                for(int i=0;i<m*n*mn;i++){
                    int index=i+(k_other-1)*m*n*mn;
                    center+=particles.X(index);}
                center/=m*n*mn;
                for(int i=0;i<m*n*mn;i++){
                    int index=i+(k_other-1)*m*n*mn;
                    particles.X(index).y+=50-center.y;}
                for(int i=0;i<m;i++)
                for(int j=0;j<n;j++)
                for(int ij=0;ij<mn;ij++) particles.V(i+m*j+m*n*ij+k_other*m*n*mn)=TV();}}}
}
void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) override {}
//#####################################################################
// Function Add_External_Forces
//#####################################################################
void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) override
{
    T v0=4,v1=6;
    if(test_number==50 || test_number==51){
        DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
        ARRAY<bool> use(particles.X.m);
        use.Subset(externally_forced).Fill(true);
        for(int p=0;p<particles.X.m;p++){
            T height=particles.X(p).x;
            if(!use(p) && height<10) continue;
            T force_multiplier=sqr(max((T)1,time-(T)1));
            if(height>8.8) force_multiplier*=2;
            if(height>14+time) continue;
            T vm=fish_V(p).Magnitude();
            if(vm>v1) continue;
            if(vm>v0) force_multiplier*=(vm-v0)/(v1-v0);
            F(p)+=TV(1.0*force_multiplier*(14+time-height),0,0);}}
}
//#####################################################################
// Function Add_Constitutive_Model
//#####################################################################
void Add_Constitutive_Model(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,T stiffness,T poissons_ratio,T damping)
{
    if(input_poissons_ratio!=-1) poissons_ratio=input_poissons_ratio;
    if(input_youngs_modulus!=0) stiffness=input_youngs_modulus;
    solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new COROTATED_FIXED<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier)));

    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    if(damping && damping_multiplier){
        DEFORMABLES_FORCES<TV>* force=Create_Finite_Volume(tetrahedralized_volume,new COROTATED_FIXED<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier));
        force->Update_Position_Based_State(0,true,true);
        solid_body_collection.Add_Force(new RALEIGH_DAMPING_FORCE<TV>(particles,force,damping*damping_multiplier,1,save_dt));}
}
//#####################################################################
// Function Add_Gravity
//#####################################################################
GRAVITY<TV>& Add_Gravity()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    GRAVITY<TV>* g=new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true);
    solid_body_collection.Add_Force(g);
    return *g;
}
};
}
#endif

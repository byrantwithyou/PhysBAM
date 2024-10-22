//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS_2D
//#####################################################################
//  100. What is it
//#####################################################################
#ifndef __STANDARD_TESTS_2D__
#define __STANDARD_TESTS_2D__

#include <Core/Log/LOG.h>
#include <Core/Math_Tools/pow.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Images/PNG_FILE.h>
#include <Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Images/TEX_FILE.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/SMOOTH_LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Deformables/Bindings/RIGID_BODY_BINDING.h>
#include <Deformables/Bindings/SOFT_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/ELASTIC_ETHER_DRAG.h>
#include <Deformables/Forces/FINITE_VOLUME.h>
#include <Deformables/Forces/RALEIGH_DAMPING_FORCE.h>
#include <Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
#include <Solids/Forces_And_Torques/ETHER_DRAG.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
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

extern int siggraph_hack_newton_iterations;

template<class TV> class STANDARD_TESTS;
template<class T_input>
class STANDARD_TESTS<VECTOR<T_input,2> >:public STANDARD_TESTS_BASE<VECTOR<T_input,2> >
{
    typedef T_input T;
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;typedef VECTOR<T,3> TV3;
public:
    typedef STANDARD_TESTS_BASE<TV> BASE;
    using BASE::solids_parameters;using BASE::viewer_dir;
    using BASE::last_frame;using BASE::frame_rate;
    using BASE::solid_body_collection;using BASE::stream_type;
    using BASE::solids_evolution;using BASE::test_number;
    using BASE::data_directory;using BASE::m;using BASE::s;using BASE::kg;
    using BASE::unit_J;using BASE::unit_rho;using BASE::unit_p;
    using BASE::backward_euler_evolution;using BASE::user_last_frame;
    
    SOLIDS_STANDARD_TESTS<TV> tests;

    T attachment_velocity;
    ARRAY<int> kinematic_ids;
    ARRAY<INTERPOLATION_CURVE<T,FRAME<TV> > > curves;
    INTERPOLATION_CURVE<T,T> scalar_curve;
    bool print_matrix;
    int resolution;
    T stiffness_multiplier;
    T damping_multiplier;
    ARRAY<int> externally_forced;
    ARRAY<int> constrained_particles;
    ARRAY<TV> constrained_velocities;
    TV_INT image_size;
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
    T save_dt;
    T final_x;
    T ether_drag;
    CELL_ITERATOR<TV>* cell_iterator;
    ARRAY<TV3,TV_INT> image;
    ARRAY<int,TV_INT> raw_image;
    GRID<TV> image_grid;
    ARRAY<INTERPOLATION_CURVE<T,TV> > kinematic_particle_positions;
    ARRAY<int> kinematic_particle_ids;
    ARRAY<TV> initial_positions;
    std::string image_file;
    std::string raw_image_file;
    int threads;

    STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :BASE(stream_type_input,parse_args),tests(stream_type_input,data_directory,solid_body_collection),
        print_matrix(false),resolution(0),stiffness_multiplier(1),damping_multiplier(1),image_size(500,500),rand_seed(1234),
        use_rand_seed(false),use_residuals(false),project_nullspace(false),
        use_penalty_collisions(false),use_constraint_collisions(true),penalty_collisions_stiffness((T)1e4),
        penalty_collisions_separation((T)1e-4),penalty_collisions_length(1),enforce_definiteness(false),
        density(pow<TV::m>(10)),final_x((T)-16.175),ether_drag(0),
        cell_iterator(0),image_file("image.png"),raw_image_file("image.txt")
    {
        parse_args.Add("-resolution",&resolution,"resolution","resolution used by multiple tests to change the parameters of the test");
        parse_args.Add("-stiffen",&stiffness_multiplier,"multiplier","stiffness multiplier for various tests");
        parse_args.Add("-dampen",&damping_multiplier,"multiplier","damping multiplier for various tests");
        parse_args.Add("-residuals",&use_residuals,"print residuals during timestepping");
        parse_args.Add("-print_matrix",&print_matrix,"print Krylov matrix");
        parse_args.Add("-project_nullspace",&project_nullspace,"project out nullspace");
        parse_args.Add("-seed",&rand_seed,&use_rand_seed,"seed","random seed to use");
        parse_args.Add("-ether_drag",&ether_drag,"drag","Ether drag");
        parse_args.Add("-image_size",&image_size,"size","image size for plots");
        parse_args.Add("-final_x",&final_x,"position","final x position");
        parse_args.Add("-use_penalty",&use_penalty_collisions,"use penalty collisions");
        parse_args.Add_Not("-no_constraints",&use_constraint_collisions,"disable constrained optimization for collisions");
        parse_args.Add("-penalty_stiffness",&penalty_collisions_stiffness,"tol","penalty collisions stiffness");
        parse_args.Add("-penalty_separation",&penalty_collisions_separation,"tol","penalty collisions separation");
        parse_args.Add("-penalty_length",&penalty_collisions_length,"tol","penalty collisions length scale");
        parse_args.Add("-enf_def",&enforce_definiteness,"enforce definiteness in system");
        parse_args.Add("-image_file",&image_file,"file","output image filename");
        parse_args.Add("-raw_image_file",&raw_image_file,"file","octave output image filename");
        parse_args.Parse();

        LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
        if(!this->user_output_directory)
            viewer_dir.output_directory=LOG::sprintf("Test_%d",test_number);
        if(use_rand_seed) rand.Set_Seed(rand_seed);
        solids_parameters.implicit_solve_parameters.project_nullspace_frequency=project_nullspace;
        penalty_collisions_length*=m;
        penalty_collisions_separation*=m;
        ether_drag/=s;
        penalty_collisions_stiffness*=unit_J;
        density*=unit_rho;

        if(use_constraint_collisions) use_penalty_collisions=false;
        if(use_penalty_collisions || use_constraint_collisions){
            solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=false;
            solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=false;
            solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;}

        solid_body_collection.Print_Residuals(use_residuals);
    }

    virtual ~STANDARD_TESTS()
    {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    bool automatically_add_to_collision_structures=true;
    // deformable bodies
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    switch(test_number){
        case 1:{
            if(!resolution) resolution=10;
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,4)*m));
            tests.Create_Mattress(GRID<TV>(TV_INT()+resolution+1,RANGE<TV>::Centered_Box()),true,&initial_state,density);
            tests.Add_Ground();
            break;}
        case 100:{
            image_grid.Initialize(image_size,RANGE<TV>::Centered_Box()*5,true);
            cell_iterator=new CELL_ITERATOR<TV>(image_grid);
            image.Resize(image_grid.Cell_Indices());
            raw_image.Resize(image_grid.Cell_Indices());
            this->fixed_dt=1./24;
            if(!this->user_frame_rate) frame_rate=1/(this->fixed_dt*image_size.Product());
            if(!user_last_frame) last_frame=1;
            for(int i=0;i<9;i++) if(i!=4) constrained_particles.Append(i);
            constrained_velocities.Resize(constrained_particles.m,use_init,TV());
            tests.Create_Mattress(GRID<TV>(TV_INT(3,2),RANGE<TV>(TV(),TV(2,1))),true,0,density);
            break;}
        case 101:{
            kinematic_particle_ids.Append(0);
            kinematic_particle_ids.Append(2);
            kinematic_particle_ids.Append(3);
            kinematic_particle_ids.Append(4);
            kinematic_particle_ids.Append(5);
            kinematic_particle_positions.Resize(5);
            kinematic_particle_positions(0).Add_Control_Point(0,TV(0,0));
            kinematic_particle_positions(1).Add_Control_Point(0,TV(2,0));
            kinematic_particle_positions(2).Add_Control_Point(0,TV(0,1));
            kinematic_particle_positions(3).Add_Control_Point(0,TV(1,1));
            kinematic_particle_positions(4).Add_Control_Point(0,TV(2,1));
            kinematic_particle_positions(0).Add_Control_Point(1,TV(0,0));
            kinematic_particle_positions(1).Add_Control_Point(1,TV(final_x,0));
            kinematic_particle_positions(2).Add_Control_Point(1,TV(0,1));
            kinematic_particle_positions(3).Add_Control_Point(1,TV(final_x/2,1));
            kinematic_particle_positions(4).Add_Control_Point(1,TV(final_x,1));
            tests.Create_Mattress(GRID<TV>(TV_INT(3,2),RANGE<TV>(TV(),TV(2,1))),true,0,density);
            break;}
        case 102:{
            image_grid.Initialize(image_size,RANGE<TV>::Centered_Box()*5,true);
            cell_iterator=new CELL_ITERATOR<TV>(image_grid);
            image.Resize(image_grid.Cell_Indices());
            raw_image.Resize(image_grid.Cell_Indices());
            if(!this->user_frame_rate) frame_rate=1/(this->fixed_dt*image_size.Product());
            if(!user_last_frame) last_frame=1;
            initial_positions=particles.X;
            break;}
        case 103:{
            INTERPOLATION_CURVE<T,VECTOR<T,2> > image_curve;
            image_curve.Add_Control_Point(0*24,VECTOR<T,2>(1,0));
            image_curve.Add_Control_Point(1*24,VECTOR<T,2>(2,0));
            image_curve.Add_Control_Point(2*24,VECTOR<T,2>(2,2));
            image_curve.Add_Control_Point(4*24,VECTOR<T,2>(-3,2));
            image_curve.Add_Control_Point(5*24,VECTOR<T,2>(-3,0));
            image_curve.Add_Control_Point(8*24,VECTOR<T,2>(1,-2));
            image_curve.Add_Control_Point(9*24,VECTOR<T,2>(1,0));
            TRIANGULATED_AREA<T>& ta=tests.Create_Mattress(GRID<TV>(TV_INT(3,6),RANGE<TV>(TV(-1,-1),TV(1,4))),true,0,density);
            SEGMENT_MESH& sm=ta.Get_Segment_Mesh();
            initial_positions=particles.X;
            Create_Directory("initial");
            ARRAY<VECTOR<T,3>,TV_INT> in_image;
            PNG_FILE<T>::Read("color_template.png",in_image);

            for(int i=0;i<=image_curve.control_points.Last().t;i++)
            {
                VECTOR<T,2> X=image_curve.Value(i);
//                VECTOR<T,2> X(-1+2*cos(i*pi*2/96),2*sin(i*pi*2/96));
                MATRIX<T,2> M(X,X.Orthogonal_Vector());
                for(int i=0;i<9;i++){
                    particles.X(i)=M*initial_positions(i);
                    particles.X(i+9)=particles.X(i)+TV(0,3);}
                char buff[100];
                sprintf(buff,"initial/initial-%03d.tex",i);
                TEX_FILE<T> tex(buff,RANGE<TV>::Unit_Box()*480);
                tex.cur_format.line_style=0;
                tex.cur_format.fill_style=1;
                tex.cur_format.fill_color=VECTOR<T,3>(.8,.8,.8);
//                tex.Draw_Object(RANGE<TV>::Centered_Box()*5.1);
                tex.cur_format.line_style=1;
                tex.cur_format.fill_style=0;

                tex.cur_format.line_width=.02;
                tex.cur_format.arrow_style="c-c";
                tex.cur_format.line_color=VECTOR<T,3>(0,0,0);
                for(int e=0;e<sm.elements.m;e++)
                    tex.Draw_Object(particles.X(sm.elements(e).x),particles.X(sm.elements(e).y));

                tex.cur_format.line_style=0;
                tex.cur_format.fill_style=1;
                tex.cur_format.fill_color=VECTOR<T,3>(0,0,0);
                for(int p=0;p<particles.X.m;p++)
                    tex.Draw_Object(particles.X(p),.05);
                tex.cur_format.line_style=1;
                tex.cur_format.fill_style=0;

                tex.cur_format.line_color=VECTOR<T,3>(0,0,0);
                VECTOR<T,3> col=in_image(TV_INT(X/5*240+240));
                tex.cur_format.fill_color=col*2;
                tex.cur_format.line_style=1;
                tex.cur_format.fill_style=1;
                tex.Draw_Object(X,(T).2);
                tex.cur_format.line_style=1;
                tex.cur_format.fill_style=0;
                tex.bounding_box=RANGE<TV>::Centered_Box()*5;
            }
            exit(0);
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
            TRIANGULATED_AREA<T>& ta=deformable_body_collection.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Gravity();
            Add_Constitutive_Model(ta,(T)1e5*unit_p,(T).45,(T)0*s);
            break;}
        case 100:{
            TRIANGULATED_AREA<T>& ta=deformable_body_collection.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(ta,(T)1e5*unit_p,(T).45,(T)0*s);
            break;}
        case 101:{
            TRIANGULATED_AREA<T>& ta=deformable_body_collection.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(ta,(T)1e5*unit_p,(T).45,(T)0*s);
            break;}
        case 102:{
            TRIANGULATED_AREA<T>& ta=deformable_body_collection.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(ta,(T)1e5*unit_p,(T).45,(T)0*s);
            break;}
        default:
            LOG::cerr<<"Missing bodies implementation for test number "<<test_number<<std::endl;exit(1);}

    if(ether_drag) solid_body_collection.Add_Force(new ELASTIC_ETHER_DRAG<TV>(deformable_body_collection.particles,true,ether_drag,1,save_dt));

    if(use_penalty_collisions)
        for(int b=0;b<rigid_body_collection.rigid_body_particles.number;b++){
            IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> > *iot=rigid_body_collection.Rigid_Body(b).implicit_object;
            if(LEVELSET_IMPLICIT_OBJECT<TV>* lio=dynamic_cast<LEVELSET_IMPLICIT_OBJECT<TV>*>(iot->object_space_implicit_object))
                iot=new IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >(new SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV>(lio->levelset.grid,lio->levelset.phi),true,iot->transform);
            solid_body_collection.Add_Force(new IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>(particles,
                    iot,penalty_collisions_stiffness,penalty_collisions_separation,penalty_collisions_length));}
    else if(use_constraint_collisions && backward_euler_evolution)
        for(int b=0;b<rigid_body_collection.rigid_body_particles.number;b++){
            IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> > *iot=rigid_body_collection.Rigid_Body(b).implicit_object;
            if(LEVELSET_IMPLICIT_OBJECT<TV>* lio=dynamic_cast<LEVELSET_IMPLICIT_OBJECT<TV>*>(iot->object_space_implicit_object)){
                lio->levelset.interpolation=new CUBIC_SPLINE_INTERPOLATION_UNIFORM<TV,T,FACE_LOOKUP_UNIFORM<TV> >();
                iot=new IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >(new SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV>(lio->levelset.grid,lio->levelset.phi),true,iot->transform);}
            backward_euler_evolution->minimization_objective.collision_objects.Append(iot);
            backward_euler_evolution->minimization_objective.coefficient_of_friction.Append(rigid_body_collection.Rigid_Body(b).coefficient_of_friction);}
    else
        for(int i=0;i<deformable_body_collection.structures.m;i++){
            deformable_body_collection.collisions.collision_structures.Append(deformable_body_collection.structures(i));
            if(solids_parameters.triangle_collision_parameters.perform_self_collision)
                solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(deformable_body_collection.structures(i));}

    if(enforce_definiteness) solid_body_collection.Enforce_Definiteness(true);

    this->After_Initialize_Bodies();
}

void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) override {}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override
{
    V.Subset(constrained_particles)=constrained_velocities;
    for(int i=0;i<kinematic_particle_ids.m;i++)
        V(kinematic_particle_ids(i))=kinematic_particle_positions(i).Derivative(velocity_time);
}
void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > rotation,const T time) override {}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) override
{
    for(int i=0;i<kinematic_particle_ids.m;i++)
        X(kinematic_particle_ids(i))=kinematic_particle_positions(i).Value(time);
}
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) override {}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override
{
    V.Subset(constrained_particles).Fill(TV());
    V.Subset(externally_forced).Fill(TV());
    V.Subset(kinematic_particle_ids).Fill(TV());
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

    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    save_dt=dt;

    if(test_number==100) particles.X(4)=cell_iterator->Location();
    if(test_number==102){
        Add_Debug_Particle(cell_iterator->Location(),TV3(1,0,0));
        MATRIX<T,2> M(cell_iterator->Location(),cell_iterator->Location().Orthogonal_Vector());
        for(int i=0;i<9;i++){
            particles.X(i)=M*initial_positions(i);
            particles.X(i+9)=particles.X(i)+TV(0,3);}
        particles.V.Fill(TV());}
}
//#####################################################################
// Function Postprocess_Substep
//#####################################################################
void Postprocess_Substep(const T dt,const T time) override
{
    if(test_number==100){
        INTERPOLATED_COLOR_MAP<T> cm;
        cm.Initialize_Colors(0,20,false,true,false);
        if(siggraph_hack_newton_iterations>=0) image(cell_iterator->index)=cm(siggraph_hack_newton_iterations);
        else image(cell_iterator->index)=TV3(.5,.5,.5);
        raw_image(cell_iterator->index)=siggraph_hack_newton_iterations;
        cell_iterator->Next();
        if(!cell_iterator->Valid()){
            PNG_FILE<T>::Write(image_file.c_str(),image);
            OCTAVE_OUTPUT<T>(raw_image_file.c_str()).Write("image",raw_image);}}
    if(test_number==102){
        INTERPOLATED_COLOR_MAP<T> cm;
        cm.Initialize_Colors(0,50,false,true,false);
        if(siggraph_hack_newton_iterations>=0) image(cell_iterator->index)=cm(siggraph_hack_newton_iterations);
        else image(cell_iterator->index)=TV3(.5,.5,.5);
        raw_image(cell_iterator->index)=siggraph_hack_newton_iterations;
        cell_iterator->Next();
        if(!cell_iterator->Valid()){
            PNG_FILE<T>::Write(image_file.c_str(),image);
            OCTAVE_OUTPUT<T>(raw_image_file.c_str()).Write("image",raw_image);}}

    BASE::Postprocess_Substep(dt,time);
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) override
{
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) override
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    for(int i=0;i<constrained_particles.m;i++) Add_Debug_Particle(particles.X(constrained_particles(i)),TV3(1,0,0));
    for(int i=0;i<externally_forced.m;i++) Add_Debug_Particle(particles.X(externally_forced(i)),TV3(0,1,0));
    for(int i=0;i<kinematic_particle_ids.m;i++) Add_Debug_Particle(particles.X(kinematic_particle_ids(i)),TV3(1,1,0));
}
//#####################################################################
// Function Add_Constitutive_Model
//#####################################################################
void Add_Constitutive_Model(TRIANGULATED_AREA<T>& ta,T stiffness,T poissons_ratio,T damping)
{
    solid_body_collection.Add_Force(Create_Finite_Volume(ta,new COROTATED_FIXED<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier)));

    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    if(damping && damping_multiplier){
        DEFORMABLES_FORCES<TV>* force=Create_Finite_Volume(ta,new COROTATED_FIXED<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier));
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

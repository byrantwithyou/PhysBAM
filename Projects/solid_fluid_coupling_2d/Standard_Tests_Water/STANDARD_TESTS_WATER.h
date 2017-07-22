//#####################################################################
// Copyright 2007-2008, Jon Gretarsson, Avi Robinson-Mosher, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS_WATER
//#####################################################################
//   1. Falling deformable sphere and smoke
//   3. Cork
//   4. Many rigid bodies
//   5. Heavy rigid block falling into water and forcing up a jet
//   6. thin shell rigid boat hit with one light sphere and one heavy sphere
//   7. thin shell cloth hit with a heavy sphere
//   8. balloon filling with water then being dropped
//   9. heavy ball flying into water
//  10. buoyant balls
//#####################################################################
#ifndef __STANDARD_TESTS_WATER__
#define __STANDARD_TESTS_WATER__

#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/FINITE_VOLUME.h>
#include <Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_2D.h>
#include <Deformables/Forces/LINEAR_SPRINGS.h>
#include <Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Dynamics/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include <Dynamics/Standard_Tests/WATER_STANDARD_TESTS_2D.h>
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS_WATER:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,2> >
{
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;using BASE::output_directory;using BASE::last_frame;using BASE::frame_rate;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions;using BASE::mpi_world; // silence -Woverloaded-virtual
    using BASE::Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization;using BASE::solid_body_collection;using BASE::test_number;using BASE::resolution;
    using BASE::Add_To_Fluid_Simulation;using BASE::Add_Volumetric_Body_To_Fluid_Simulation;using BASE::Add_Thin_Shell_To_Fluid_Simulation;using BASE::Adjust_Phi_With_Source;
    using BASE::data_directory;

    WATER_STANDARD_TESTS_2D<TV> water_tests;
    SOLIDS_STANDARD_TESTS<TV> solids_tests;
    GRID<TV> mattress_grid;
    int deformable_object_id;
    T solid_density;
    int light_sphere_index,lightish_sphere_index,neutral_sphere_index,heavy_sphere_index;
    T light_sphere_initial_height,heavy_sphere_initial_height;
    T light_sphere_drop_time,heavy_sphere_drop_time;
    T initial_water_level;
    ARRAY<RANGE<TV> > fountain_source;
    ARRAY<TV> fountain_source_velocity;
    MATRIX<T,3> world_to_source;
    int left_fixed_index,right_fixed_index;
    int bodies;
    int rigid_body_id;
    T solid_scale;
    ARRAY<int> rigid_bodies_to_simulate;
    
    T rotation_angle;
    T velocity_angle;
    T velocity_multiplier;
    T mass_multiplier;
    bool flow_particles;
    bool solid_node;
    bool mpi;

    STANDARD_TESTS_WATER(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :BASE(stream_type_input,parse_args,1,fluids_parameters.WATER),water_tests(*this,fluids_parameters,solid_body_collection.rigid_body_collection),
        solids_tests(stream_type_input,data_directory,solid_body_collection),deformable_object_id(0),solid_density(1),light_sphere_index(0),lightish_sphere_index(0),neutral_sphere_index(0),
        heavy_sphere_index(0),light_sphere_initial_height((T)2),heavy_sphere_initial_height((T)2.25),light_sphere_drop_time((T).7),heavy_sphere_drop_time((T)1.25),
        left_fixed_index(0),right_fixed_index(0),bodies(5),solid_scale((T).02),velocity_multiplier(1),mass_multiplier(1),flow_particles(false)
    {
        solid_node=mpi_world->initialized && !mpi_world->rank;
        mpi=mpi_world->initialized;

        parse_args.Add("-bodies",&bodies,"value","bodies");
        parse_args.Add("-mass",&solid_density,"value","solid_density");
        parse_args.Add("-scale",&solid_scale,"value","solid_scale");
        parse_args.Add("-slip",&fluids_parameters.use_slip,"use slip");
        parse_args.Parse();

        solids_tests.data_directory=data_directory;
        water_tests.Initialize(Water_Test_Number(test_number),resolution);
        LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
        last_frame=1000;

        //mattress_grid=GRID_2D<T>(4,2,(T).1,(T).49,(T).4,(T).6);
        solids_parameters.triangle_collision_parameters.perform_self_collision=false;
        solids_parameters.rigid_body_collision_parameters.use_push_out=false;

        fluids_parameters.solve_neumann_regions=false;
        fluids_parameters.incompressible_iterations=400;
        *fluids_parameters.grid=water_tests.grid;
        fluids_parameters.fluid_affects_solid=fluids_parameters.solid_affects_fluid=true;
        fluids_parameters.second_order_cut_cell_method=false;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.use_trapezoidal_rule_for_velocities=false;
        if(solid_node || !mpi) solids_parameters.use_rigid_deformable_contact=true;

        fluids_parameters.use_vorticity_confinement=false;
        fluids_parameters.use_preconditioner_for_slip_system=true;
        
        fluids_parameters.use_removed_positive_particles=true;fluids_parameters.use_removed_negative_particles=true;
        fluids_parameters.write_removed_positive_particles=true;fluids_parameters.write_removed_negative_particles=true;
        fluids_parameters.store_particle_ids=true;
        // T default_removed_positive_particle_buoyancy_constant=fluids_parameters.removed_positive_particle_buoyancy_constant;
        fluids_parameters.removed_positive_particle_buoyancy_constant=0;
        
        fluids_parameters.gravity=TV();
        
        solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_symmqmr;
        
        switch(test_number){
            case 1:
                fluids_parameters.density=(T)1000;
                solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
                mattress_grid=GRID<TV>(TV_INT(20,10),RANGE<TV>(TV((T).4,(T).8),TV((T).8,(T)1.04)));
                solid_density=(T)1200;
                break;
            case 2:
                last_frame=1000;
                fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;
                break;
            case 3:
                fluids_parameters.density=(T)1000;
                solids_parameters.implicit_solve_parameters.cg_iterations=400;
                fluids_parameters.domain_walls[0][0]=fluids_parameters.domain_walls[0][1]=fluids_parameters.domain_walls[1][0]=true;//fluids_parameters.domain_walls[1][1]=true;

                (*fluids_parameters.grid).Initialize(TV_INT(2*resolution+1,5*resolution+1),RANGE<TV>(TV((T)-2,(T)0),TV((T)2,(T)10)));
                break;
            case 4:
                solids_parameters.rigid_body_collision_parameters.use_push_out=true;
                fluids_parameters.density=(T)1000;
                solids_parameters.rigid_body_collision_parameters.enforce_rigid_rigid_contact_in_cg=true;
                fountain_source.Append(RANGE<TV>(TV((T).45,(T).3),TV((T).55,(T).35)));
                //fountain_source.Append(RANGE<TV>((T).45,(T).55,(T).25,(T).3));
                //fountain_source.Append(RANGE<TV>((T).975,(T)1,(T).1,(T).4)); // TODO: make total input/output sourcing add up to zero
                //fountain_source.Append(RANGE<TV>((T)0,(T).025,(T).1,(T).4));
                //fountain_source.Append(RANGE<TV>((T).1,(T).4,(T).0,(T).05));
                //fountain_source.Append(RANGE<TV>((T).6,(T).9,(T).0,(T).05));
                fountain_source_velocity.Append(TV((T)0,(T)4));
                //fountain_source_velocity.Append(TV((T)0,(T)-2));
                //fountain_source_velocity.Append(TV((T)0,(T)-2));
                //fountain_source_velocity.Append(TV((T)2,(T)0));
                //fountain_source_velocity.Append(TV((T)-2,(T)0));
                world_to_source=MATRIX<T,3>::Identity_Matrix();
                fluids_parameters.reseeding_frame_rate=10;
                //solid_density=(T)1100;
                (*fluids_parameters.grid).Initialize(TV_INT(10*resolution+1,20*resolution+1),RANGE<TV>(TV((T)0,(T)0),TV((T)1,(T)2)));
                break;
            case 5:
                fluids_parameters.density=(T)1000;
                solid_density=(T)5000;
                solids_parameters.implicit_solve_parameters.cg_iterations=400;
                (*fluids_parameters.grid).Initialize(TV_INT(10*resolution+1,20*resolution+1),RANGE<TV>(TV((T)0,(T)0),TV((T)1,(T)2)));
                fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;
                break;
            case 6:
                fluids_parameters.density=(T)1000;
                solids_parameters.implicit_solve_parameters.cg_iterations=400;
                solid_density=(T)10000;
                (*fluids_parameters.grid).Initialize(TV_INT(10*resolution+1,20*resolution+1),RANGE<TV>(TV((T)0,(T)0),TV((T)1,(T)2)));
                fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=false;
                break;
            case 7:
                solids_parameters.triangle_collision_parameters.perform_self_collision=true;
                fluids_parameters.density=(T)1000;
                solid_density=(T)5000;
                solids_parameters.implicit_solve_parameters.cg_iterations=400;
                (*fluids_parameters.grid).Initialize(TV_INT(10*resolution+1,25*resolution+1),RANGE<TV>(TV((T)0,(T)0),TV((T)1,(T)2.5)));
                fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=false;
                break;
            case 8:
                solids_parameters.triangle_collision_parameters.perform_self_collision=true;
                light_sphere_initial_height=(T)1.15;
                light_sphere_drop_time=(T).7;
                heavy_sphere_drop_time=(T)2;
                fluids_parameters.density=(T)10000;
                solids_parameters.implicit_solve_parameters.cg_iterations=400;
                (*fluids_parameters.grid).Initialize(TV_INT(10*resolution+1,25*resolution+1),RANGE<TV>(TV((T)0,(T)0),TV((T)1,(T)2.5)));
                fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=false;
                break;
            case 9:
                last_frame=200;
                (*fluids_parameters.grid).Initialize(TV_INT(15*resolution+1,10*resolution+1),RANGE<TV>(TV(0,0),TV((T)1.5,1)));
                heavy_sphere_drop_time=(T).05;
                light_sphere_drop_time=(T).2;
                break;
            case 10:
                last_frame=200;
                fluids_parameters.density=(T)1000;
                (*fluids_parameters.grid).Initialize(TV_INT(40*resolution+1,20*resolution+1),RANGE<TV>(TV((T)-2,(T)0),TV((T)2,(T)2)));
                fluids_parameters.domain_walls[1][1]=false;
                break;
            case 12:
                fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;
                (*fluids_parameters.grid).Initialize(TV_INT(4*resolution+1,resolution+1),RANGE<TV>(TV((T)0,(T)0),TV((T)4,(T)1)));
                // (*fluids_parameters.grid).Initialize(TV_INT(4*resolution+1,2*resolution+1),RANGE<TV>(TV((T)0,(T)-0.5),TV((T)4,(T)1.5)));
                rotation_angle=-pi/16;
                solids_parameters.implicit_solve_parameters.cg_iterations=200;
            
                flow_particles=true;
            
                break;
            case 13:
                fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;fluids_parameters.domain_walls[1][0]=false;fluids_parameters.domain_walls[1][1]=false;
                (*fluids_parameters.grid).Initialize(TV_INT(2*resolution+1,resolution+1),RANGE<TV>(TV((T)0,(T)0),TV((T)2,(T)1)));
                velocity_angle=pi/16;
                solids_parameters.implicit_solve_parameters.cg_iterations=10000;
                fluids_parameters.density=(T)1000;
                break;
            default:
                LOG::cerr<<"Unrecognized test number "<<test_number<<std::endl;exit(1);}

        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Rigid_Body_Walls(*this);
        // output_directory=LOG::sprintf("Standard_Tests_Water/Test_%d_Resolution_%d_density_%d",test_number,resolution,solid_density);
        if(fluids_parameters.use_slip)
            output_directory=LOG::sprintf("Standard_Tests_Water/Test_%d_Resolution_%d_slip",test_number,resolution);
        else
            output_directory=LOG::sprintf("Standard_Tests_Water/Test_%d_Resolution_%d",test_number,resolution);
    }

    // Unused callbacks
    void Preprocess_Frame(const int frame) override {}
    void Postprocess_Frame(const int frame) override {
        
    if(flow_particles)
    {
        fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).removed_positive_particles(TV_INT(1,1))=fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).template_removed_particles.Clone();
        int particle_count=40;
        for(int i=0;i<particle_count;i++){
            int particle_index=fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).removed_positive_particles(TV_INT(1,1))->Add_Element();
            fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).removed_positive_particles(TV_INT(1,1))->X(particle_index)=TV((T).05,.05+(i/(T)particle_count)*0.9);
        }
        fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).Update_Particle_Cells(fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).removed_positive_particles);
    }
        
        }
    void Postprocess_Solids_Substep(const T time,const int substep) override {}
    void Apply_Constraints(const T dt,const T time) override {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) override {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) override {}
    void Update_Time_Varying_Material_Properties(const T time) override {}
    void Limit_Solids_Dt(T& dt,const T time) override {}
    void Update_Solids_Parameters(const T time) override {}
    void Preprocess_Solids_Substep(const T time,const int substep) override {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) override {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) override {}

void After_Initialization() override {BASE::After_Initialization();}
//#####################################################################
// Function Water_Test_Number
//#####################################################################
static int Water_Test_Number(const int test_number)
{
    switch(test_number){
        case 1:
        case 2:
        case 3:
        case 4:
            return 1;
        default:
            return 1;}
}
//#####################################################################
// Function Find_Placement
//#####################################################################
FRAME<TV> Find_Placement(RANDOM_NUMBERS<T>& random,const RANGE<TV>& bounding_box,ARRAY<ORIENTED_BOX<TV> >& bounding_boxes,const RANGE<TV>& world,bool want_rotate)
{
    for(int i=0;i<10000;i++){
        FRAME<TV> frame;
        if(want_rotate) frame.r=random.template Get_Rotation<TV>();
        ORIENTED_BOX<TV> oriented_box(bounding_box,frame.r);
        RANGE<TV> new_box(oriented_box.Bounding_Box());
        frame.t=random.Get_Uniform_Vector(world.min_corner-new_box.min_corner,world.max_corner-new_box.max_corner);
        oriented_box.corner+=frame.t;
        bool okay=true;
        for(int j=0;j<bounding_boxes.m;j++) if(oriented_box.Intersection(bounding_boxes(j))){okay=false;break;}
        if(okay){
            bounding_boxes.Append(oriented_box);
            return frame;}}
    PHYSBAM_FATAL_ERROR("Could not find suitable placement");
}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() override
{
    fluids_parameters.Use_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initial_Phi
//#####################################################################
T Initial_Phi(const TV& X) const
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    T phi=1;
    switch(test_number){
        case 1:
        case 6:
            phi=X.y-(T).6*fluids_parameters.grid->domain.max_corner.y;
            break;
        case 3:
            phi=X.y-fluids_parameters.grid->domain.max_corner.y+.5;//-(T)1;
            break;
        case 4:
            phi=X.y-(T).6;
            break;
        case 5:
            phi=X.y-(T).4*fluids_parameters.grid->domain.max_corner.y;
            break;
        case 7:
            phi=X.y-(T).45*fluids_parameters.grid->domain.max_corner.y;
            break;
        case 8:
            phi=X.y;
            break;
        case 9:
            phi=min(X.y-(T).400235234,rigid_body_collection.Rigid_Body(heavy_sphere_index).Implicit_Geometry_Extended_Value(X));
            break;
        case 10:{
            T surface_phi=X.y-(T)1.35;
            T geometry_phi=-min(min(rigid_body_collection.Rigid_Body(light_sphere_index).Implicit_Geometry_Extended_Value(X),
                    rigid_body_collection.Rigid_Body(lightish_sphere_index).Implicit_Geometry_Extended_Value(X)),
                min(rigid_body_collection.Rigid_Body(neutral_sphere_index).Implicit_Geometry_Extended_Value(X),
                    rigid_body_collection.Rigid_Body(heavy_sphere_index).Implicit_Geometry_Extended_Value(X)));                
            phi=max(surface_phi,geometry_phi);
            break;}
        case 12:
            phi=-1;
            break;
        case 13:
            phi=-1;
            break;
        default:
            phi=water_tests.Initial_Phi(X);}
    return phi;
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
bool Adjust_Phi_With_Sources(const T time) override
{
    if(test_number==8)
        if(time<light_sphere_drop_time){
            RANGE<TV> source_box(TV((T).45,(T)1.15),TV((T).55,(T)1.25));
            Adjust_Phi_With_Source(source_box,MATRIX<T,3>::Identity_Matrix());return true;}
    
    ARRAY<T,VECTOR<int,2> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    if(test_number==12)
        phi.Fill(-1);
    if(test_number==13)
        phi.Fill(-1);
    return false;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() override
{
    GRID<TV>& grid=*fluids_parameters.grid;
    ARRAY<T,VECTOR<int,2> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) phi(iterator.Cell_Index())=Initial_Phi(iterator.Location());
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() override
{
    for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) 
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=water_tests.Initial_Velocity(iterator.Location())[iterator.Axis()];
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,const T time) override
{
    for(int i=0;i<fountain_source.m;i++) BASE::Get_Source_Velocities(fountain_source(i),world_to_source,fountain_source_velocity(i));
    //if(test_number==4) for(int i=0;i<fountain_source.m;i++){
    //    for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) if(fountain_source(i).Lazy_Inside(world_to_source.Homogeneous_Times(iterator.Location()))){
    //        int axis=iterator.Axis();//fluids_parameters.incompressible->projection.elliptic_solver->psi_N.Component(axis)(iterator.Face_Index())=true;
    //        fluid_collection.incompressible_fluid_collection.face_velocities.Component(axis)(iterator.Face_Index())=fountain_source_velocity(i)[axis];}}
    if(test_number==12)
    {
        for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,0);iterator.Valid();iterator.Next())
        {
            int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            face_velocities(axis,face_index)=velocity_multiplier;
            psi_N(axis,face_index)=true;
        }
        // for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,2);iterator.Valid();iterator.Next())
        // {
            // int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            // face_velocities(axis,face_index)=0;
            // psi_N(axis,face_index)=true;
        // }
        // for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,3);iterator.Valid();iterator.Next())
        // {
            // int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            // face_velocities(axis,face_index)=0;
            // psi_N(axis,face_index)=true;
        // }
    }
    else if(test_number==13)
    {
        for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,0);iterator.Valid();iterator.Next())
        {
            int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            if(iterator.Location().y>(0.5-sin(velocity_angle)))
                face_velocities(axis,face_index)=velocity_multiplier*cos(velocity_angle);
            else
                face_velocities(axis,face_index)=-velocity_multiplier*cos(velocity_angle);
            if(iterator.Location().y>(0.5-sin(velocity_angle)) || iterator.Location().y<(0.5-sin(velocity_angle)-fluids_parameters.grid->DX()(2)))
                psi_N(axis,face_index)=true;
        }
        for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,2);iterator.Valid();iterator.Next())
        {
            int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            face_velocities(axis,face_index)=-velocity_multiplier*sin(velocity_angle);
            psi_N(axis,face_index)=true;
        }
        for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,1);iterator.Valid();iterator.Next())
        {
            int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            if(iterator.Location().y<(0.5+sin(velocity_angle)))
                face_velocities(axis,face_index)=-velocity_multiplier*cos(velocity_angle);
            else
                face_velocities(axis,face_index)=velocity_multiplier*cos(velocity_angle);
            if(iterator.Location().y<(0.5+sin(velocity_angle)) || iterator.Location().y>(0.5+sin(velocity_angle)+fluids_parameters.grid->DX()(2)))
                psi_N(axis,face_index)=true;
        }
        for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,3);iterator.Valid();iterator.Next())
        {
            int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            face_velocities(axis,face_index)=velocity_multiplier*sin(velocity_angle);
            psi_N(axis,face_index)=true;
        }
    }
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    int first_coupled_rigid_body=INT_MAX;

    switch(test_number){
        case 1:{
            TRIANGULATED_AREA<T>& triangulated_area=*TRIANGULATED_AREA<T>::Create();
            triangulated_area.Initialize_Herring_Bone_Mesh_And_Particles(mattress_grid);
            triangulated_area.Check_Signed_Areas_And_Make_Consistent(true);
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(triangulated_area,solid_density,true);
            solids_tests.Copy_And_Add_Structure(triangulated_area);
            solids_tests.Add_Ground();
            break;}
        case 2:{
            RIGID_BODY<TV>& rigid_body_square=solids_tests.Add_Rigid_Body("square",(T).1,(T)0);
            rigid_body_square.Frame()=FRAME<TV>(TV((T).5,(T).6));
            rigid_body_square.Set_Coefficient_Of_Restitution((T)0);
            rigid_body_square.name="square";
            rigid_body_square.Set_Mass(100);
            break;}
        case 3:{
            T solid_density=fluids_parameters.density*(T)1.0;
            T sphere_mass=sqr((T).4)*(T).2*solid_density;//(T)pi*sqr(scale)*solid_density;
            RIGID_BODY<TV>& rigid_body_sphere=solids_tests.Add_Rigid_Body("circle",(T).4,(T)0);
            rigid_body_id=rigid_body_sphere.particle_index;
            rigid_body_sphere.Update_Bounding_Box();
            rigid_body_sphere.Frame().t=TV((T)0,(T).5);
            rigid_body_sphere.Set_Mass(sphere_mass);
            break;}
        case 4:{
            // place a bunch of rigid bodies - blocks and spheres
            int num_spheres=bodies;int num_blocks=bodies;
            T scale=(T).02;
            T sphere_mass=(T)pi*sqr(scale)*solid_density,block_mass=(T)sqr(2*scale)*solid_density;
            ARRAY<ORIENTED_BOX<TV> > bounding_boxes;
            const GRID<TV>& grid=*fluids_parameters.grid;
            
            //RANGE<TV> world(grid.domain.min_corner.x,grid.domain.max_corner.x,grid.domain.min_corner.y+(grid.domain.max_corner.y-grid.domain.min_corner.y)*(T).75,(T)1.5*grid.domain.max_corner.y);
            RANGE<TV> world(TV(grid.domain.min_corner.x,(T).75),TV(grid.domain.max_corner.x,(T)1.5));

            first_coupled_rigid_body=rigid_body_collection.rigid_body_particles.Size()+1;
            RANDOM_NUMBERS<T> random;
            random.Set_Seed(1234);

            for(int i=0;i<num_spheres;i++){
                RIGID_BODY<TV>& rigid_body_sphere=solids_tests.Add_Rigid_Body("circle",scale,(T)0);
                rigid_body_sphere.Update_Bounding_Box();
                rigid_body_sphere.Frame()=Find_Placement(random,rigid_body_sphere.axis_aligned_bounding_box,bounding_boxes,world,true);
                rigid_body_sphere.Set_Mass(sphere_mass);}

            for(int i=0;i<num_blocks;i++){
                RIGID_BODY<TV>& rigid_body_block=solids_tests.Add_Rigid_Body("square",scale,(T)0);
                rigid_body_block.Update_Bounding_Box();
                rigid_body_block.Frame()=Find_Placement(random,rigid_body_block.axis_aligned_bounding_box,bounding_boxes,world,true);
                rigid_body_block.Set_Mass(block_mass);}
            break;}
        case 5:{
            RIGID_BODY<TV>& rigid_body_square=solids_tests.Add_Rigid_Body("square",(T).4,(T)0);
            rigid_body_square.Frame()=FRAME<TV>(TV((T).4,(T)1.4));
            rigid_body_square.Set_Coefficient_Of_Restitution((T)0);
            rigid_body_square.name="square";
            rigid_body_square.Update_Bounding_Box();
            T volume=sqr((T).4);
            rigid_body_square.Set_Mass(volume*solid_density);
            break;}
        case 6:{
            // boat
            RIGID_BODY<TV>& boat=solids_tests.Add_Rigid_Body("boat",(T).5,(T)0,false); // don't read implicit because it's a volumetric level set
            boat.thin_shell=true;
            boat.Frame()=FRAME<TV>(TV((T).5,(T)1.3));
            T boat_density=(T)25;
            boat.Set_Mass(boat.simplicial_object->Total_Length()*boat_density);
            // light sphere
            RIGID_BODY<TV>& light_sphere=solids_tests.Add_Rigid_Body("circle",(T).125,(T)0);
            light_sphere_index=light_sphere.particle_index;
            light_sphere.Frame()=FRAME<TV>(TV((T).35,light_sphere_initial_height));
            T light_sphere_density=(T)100;
            light_sphere.Set_Mass(light_sphere.Volume()*light_sphere_density);
            // heavy sphere
            RIGID_BODY<TV>& heavy_sphere=solids_tests.Add_Rigid_Body("circle",(T).125,(T)0);
            heavy_sphere_index=heavy_sphere.particle_index;
            heavy_sphere.Frame()=FRAME<TV>(TV((T).55,heavy_sphere_initial_height));
            heavy_sphere.Set_Mass(heavy_sphere.Volume()*solid_density);
            break;}
        case 7:{
            // heavy rigid body
            RIGID_BODY<TV>& heavy_square=solids_tests.Add_Rigid_Body("square",(T).125,(T)0);
            heavy_sphere_index=heavy_square.particle_index;
            heavy_square.Frame()=FRAME<TV>(TV((T).5,(T)1.5));
            heavy_square.Update_Bounding_Box();
            heavy_square.Set_Mass(500);

            // light cloth
            deformable_object_id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Deformable_Object(deformable_body_collection,100,TV((T).2,(T)1.15),TV((T)1,(T)1.15));
            SEGMENTED_CURVE_2D<T>& segmented_curve=dynamic_cast<SEGMENTED_CURVE_2D<T>&>(*deformable_body_collection.structures(deformable_object_id));
            THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Set_Mass(segmented_curve,10);
            segmented_curve.Initialize_Hierarchy();
            break;}
        case 8:{
            // balloon
            deformable_object_id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Deformable_Object(deformable_body_collection,100,TV((T).4,light_sphere_initial_height),TV((T).6,(T)light_sphere_initial_height));
            SEGMENTED_CURVE_2D<T>& segmented_curve=dynamic_cast<SEGMENTED_CURVE_2D<T>&>(*deformable_body_collection.structures(deformable_object_id));
            THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Set_Mass(segmented_curve,10);
            left_fixed_index=1;right_fixed_index=segmented_curve.particles.Size();
            break;}
        case 9:{
            RIGID_BODY<TV>& sphere=solids_tests.Add_Rigid_Body("circle",(T).1,(T)0);
            sphere.Frame().t=TV((T)1.25,(T).55);
            sphere.Set_Coefficient_Of_Restitution((T)0);
            sphere.coefficient_of_friction=(T)1;
            sphere.Set_Mass((T)1e10);
            heavy_sphere_index=sphere.particle_index;
            break;}
        case 10:{
            T initial_height=.5;
            T fluid_mass=(T)pi*sqr((T).325)*fluids_parameters.density;
            RIGID_BODY<TV>& light_cork=solids_tests.Add_Rigid_Body("circle",(T).325,(T)0);
            light_sphere_index=light_cork.particle_index;
            light_cork.Frame()=FRAME<TV>(TV((T)-1.5,initial_height));
            light_cork.name="Light_Cork";
            light_cork.Set_Mass((T).1*fluid_mass);
            light_cork.Set_Coefficient_Of_Restitution((T)1);
            light_cork.coefficient_of_friction=(T)1;
            rigid_bodies_to_simulate.Append(light_cork.particle_index);

            RIGID_BODY<TV>& lightish_cork=solids_tests.Add_Rigid_Body("circle",(T).325,(T)0);
            lightish_sphere_index=lightish_cork.particle_index;
            lightish_cork.Frame()=FRAME<TV>(TV((T)-0.5,initial_height));
            lightish_cork.name="Lightish_Cork";
            lightish_cork.Set_Mass((T).5*fluid_mass);
            lightish_cork.Set_Coefficient_Of_Restitution((T)1);
            lightish_cork.coefficient_of_friction=(T)1;
            rigid_bodies_to_simulate.Append(lightish_cork.particle_index);

            RIGID_BODY<TV>& buoyant_cork=solids_tests.Add_Rigid_Body("circle",(T).325,(T)0);
            neutral_sphere_index=buoyant_cork.particle_index;
            buoyant_cork.Frame()=FRAME<TV>(TV((T)0.5,initial_height));
            buoyant_cork.name="Buoyant_Cork";
            buoyant_cork.Set_Mass((T).9*fluid_mass);
            buoyant_cork.Set_Coefficient_Of_Restitution((T)1);
            buoyant_cork.coefficient_of_friction=(T)1;
            rigid_bodies_to_simulate.Append(buoyant_cork.particle_index);

            RIGID_BODY<TV>& heavy_cork=solids_tests.Add_Rigid_Body("circle",(T).325,(T)0);
            heavy_sphere_index=heavy_cork.particle_index;
            heavy_cork.Frame()=FRAME<TV>(TV((T)1.5,initial_height));
            heavy_cork.name="Heavy_Cork";
            heavy_cork.Set_Mass((T)10*fluid_mass);
            heavy_cork.Set_Coefficient_Of_Restitution((T)1);
            heavy_cork.coefficient_of_friction=(T)1;
            rigid_bodies_to_simulate.Append(heavy_cork.particle_index);
            break;}
        case 12:
            {
            // RIGID_BODY<TV>& rigid_body=solids_tests.Add_Rigid_Body("ground",(T).004,(T)0,false,true);
            // rigid_body.Frame()=FRAME<TV>(TV((T)1.0,(T)0.5));
            // rigid_body.Frame().r=ROTATION<TV>::From_Angle(rotation_angle);
            // rigid_body.Set_Coefficient_Of_Restitution((T)0);
            // rigid_body.name="square");
            // T density=2000;
            // rigid_body.Set_Mass((T)pi*(T).01*(T)density*mass_multiplier);
            // rigid_body.thin_shell=true;
            // rigid_body.is_static=true;
            
            T scale=(T).15;
            RIGID_BODY<TV>& rigid_body=solids_tests.Add_Rigid_Body("circle",scale,(T)0,true,true);
            rigid_body.Frame()=FRAME<TV>(TV((T)1.0,(T)0.5));
            // rigid_body.Frame()=FRAME<TV>(TV((T)1.0,(T)0.505));
            // rigid_body.Frame()=FRAME<TV>(TV((T)5.0,(T)0.5));
            rigid_body.Set_Coefficient_Of_Restitution((T)0);
            rigid_body.name="circle_fixed";
            T density=1000;
            rigid_body.Set_Mass(scale*scale*(T)pi*(T).01*(T)density*mass_multiplier);
            // rigid_body.thin_shell=false;
            rigid_body.is_static=true;
            
            RIGID_BODY<TV>& rigid_body_dummy=solids_tests.Add_Rigid_Body("circle",(T).05,(T)0);
            rigid_body_dummy.Frame()=FRAME<TV>(TV((T)10,(T)0));
            rigid_body_dummy.Set_Coefficient_Of_Restitution((T)0);
            rigid_body_dummy.name="circle";
            rigid_body_dummy.Set_Mass(10);
            }
            break;
        case 13:
            {
            RIGID_BODY<TV>& rigid_body=solids_tests.Add_Rigid_Body("ground",(T).012,(T)0,false,true);
            rigid_body.Frame()=FRAME<TV>(TV((T)1.0,(T)0.5));
            rigid_body.Frame().r=ROTATION<TV>::From_Angle(velocity_angle);
            rigid_body.Set_Coefficient_Of_Restitution((T)0);
            rigid_body.name="square";
            T density=2000;
            rigid_body.Set_Mass((T)pi*(T).01*(T)density*mass_multiplier);
            rigid_body.thin_shell=true;
            rigid_body.is_static=true;
            
            RIGID_BODY<TV>& rigid_body_dummy=solids_tests.Add_Rigid_Body("circle",(T).05,(T)0);
            rigid_body_dummy.Frame()=FRAME<TV>(TV((T)10,(T)0));
            rigid_body_dummy.Set_Coefficient_Of_Restitution((T)0);
            rigid_body_dummy.name="circle";
            rigid_body_dummy.Set_Mass(10);
            }
            break;
        default: LOG::cerr<<"Missing implementation for test number "<<test_number<<std::endl;exit(1);}

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.structures);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.structures.m;i++) deformable_body_collection.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();
    
#if 1
    solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true));
#endif
    switch(test_number){
        case 1:{
            TRIANGULATED_AREA<T>& triangulated_area=deformable_body_collection.template Find_Structure<TRIANGULATED_AREA<T>&>();
            std::cout<<"Setting density "<<solid_density<<std::endl;
            solid_body_collection.Add_Force(Create_Edge_Springs(triangulated_area,(T)5e2/(1+sqrt((T)2)),(T)10));
            solid_body_collection.Add_Force(Create_Altitude_Springs(triangulated_area,(T)5e2));
            DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(triangulated_area.Get_Boundary_Object());
            deformable_collisions.object.Initialize_Hierarchy();
            Add_To_Fluid_Simulation(deformable_collisions);
            break;}
        case 3:{
            Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_id));
            break;}
        case 4:{
            for(int i=first_coupled_rigid_body;i<=rigid_body_collection.rigid_body_particles.Size();i++) Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(i));
            solids_tests.Add_Ground();
            break;}
        case 5:
        case 2:{
            Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particles.Size()));
            solids_tests.Add_Ground();
            break;}
        case 6:{
            for(int p=0;p<rigid_body_collection.rigid_body_particles.Size();p++) Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(p));
            solids_tests.Add_Ground();
            break;}
        case 7:{
            // Set up cloth forces
            SEGMENTED_CURVE_2D<T>& segmented_curve=dynamic_cast<SEGMENTED_CURVE_2D<T>&>(*deformable_body_collection.structures(deformable_object_id));
            segmented_curve.Initialize_Hierarchy();
            solid_body_collection.Add_Force(Create_Edge_Springs(segmented_curve,(T)3e4,(T)1,false,(T).1,true,(T)0,true));
            solid_body_collection.Add_Force(Create_Segment_Bending_Springs(segmented_curve,(T)2/(1+sqrt((T)2)),(T)2,false,(T).1,true,(T)0,true));
            DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(segmented_curve);

            Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(heavy_sphere_index));
            Add_To_Fluid_Simulation(deformable_collisions);
            solids_tests.Add_Ground();
            this->solids_evolution->fully_implicit=true;
            break;}
        case 8:{
            SEGMENTED_CURVE_2D<T>& segmented_curve=dynamic_cast<SEGMENTED_CURVE_2D<T>&>(*deformable_body_collection.structures(deformable_object_id));
            segmented_curve.Initialize_Hierarchy();
            solid_body_collection.Add_Force(Create_Edge_Springs(segmented_curve,(T)2e3,(T)4,false,(T).1,true,(T)0,true));
            solid_body_collection.Add_Force(Create_Segment_Bending_Springs(segmented_curve,(T)2/(1+sqrt((T)2)),(T)2,false,(T).1,true,(T)0,true));
            DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(segmented_curve);
            Add_To_Fluid_Simulation(deformable_collisions);
            solids_tests.Add_Ground();
            this->solids_evolution->fully_implicit=true;
            break;}
        case 9:{
            Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(heavy_sphere_index));
            solids_tests.Add_Ground();
            break;}
        case 10:{
            for(int i=0;i<rigid_bodies_to_simulate.m;i++){
                RIGID_BODY<TV>& rigid_body_to_add=rigid_body_collection.Rigid_Body(rigid_bodies_to_simulate(i));
                if(rigid_body_to_add.thin_shell) Add_Thin_Shell_To_Fluid_Simulation(rigid_body_to_add); else Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_to_add);}
            break;}
        case 12:{
            // solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true,solid_gravity));
            // Add_Thin_Shell_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particles.Size()-1));
            Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particles.Size()-1));
            Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particles.Size()));
            break;}
        case 13:{
            // solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true,solid_gravity));
            // Add_Thin_Shell_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particles.Size()-1));
            Add_Thin_Shell_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particles.Size()-1));
            Add_Volumetric_Body_To_Fluid_Simulation(rigid_body_collection.Rigid_Body(rigid_body_collection.rigid_body_particles.Size()));
            
            // TRIANGULATED_AREA<T>& triangulated_area=deformable_body_collection.template Find_Structure<TRIANGULATED_AREA<T>&>();
            // triangulated_area.Initialize_Hierarchy();
            // solid_body_collection.Add_Force(Create_Finite_Volume(triangulated_area,new NEO_HOOKEAN<T,2>((T)1e3,(T).45,(T).01)));
            // DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(triangulated_area.Get_Boundary_Object());
            // deformable_collisions.object.Initialize_Hierarchy();
            // Add_To_Fluid_Simulation(deformable_collisions,true,false);
            break;}
        default:break;}
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) override
{
    if(test_number==6){
        if(time<light_sphere_drop_time) frame(light_sphere_index).t.y=light_sphere_initial_height; // drop ball 1 at time .35 takes ~.5 seconds to fall
        if(time<heavy_sphere_drop_time) frame(heavy_sphere_index).t.y=heavy_sphere_initial_height;} // drop ball 2 at time 1
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) override
{
    if(test_number==8)
        if(time<heavy_sphere_drop_time)
        {X(left_fixed_index)=TV((T).4,light_sphere_initial_height);X(right_fixed_index)=TV((T)0.6,light_sphere_initial_height);}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) override
{
    if(test_number==9)
        if(velocity_time<heavy_sphere_drop_time) twist(heavy_sphere_index).linear=-TV((T)2.25,(T)2.25);
    if(test_number==6){
        if(velocity_time<light_sphere_drop_time) twist(light_sphere_index).linear.y=(T)0;
        if(velocity_time<heavy_sphere_drop_time) twist(heavy_sphere_index).linear.y=(T)0;}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override
{
    if(test_number==8)
        if(velocity_time<heavy_sphere_drop_time)
        {V(left_fixed_index)=TV((T)0,(T)0);V(right_fixed_index)=TV((T)0,(T)0);}
}
//#####################################################################
};
}
#endif

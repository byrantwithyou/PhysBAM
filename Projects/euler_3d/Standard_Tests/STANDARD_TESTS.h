//#####################################################################
// Copyright 2009 Jon Gretarsson, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <fstream>
#include <iostream>
#include "math.h"

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_REFLECTION_ATTENUATION.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Tessellation/RANGE_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER_HASH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Fracture/FRACTURE_PATTERN.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Fracture/FRACTURE_REGION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGIDS_NEWMARK_COLLISION_CALLBACKS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_RF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS_SMOOTH_TRANSITION_INCOMPRESSIBLE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_SLIP.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>

namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
public:
    typedef T_input T;typedef VECTOR<T,3> TV;typedef GRID<TV> T_GRID;typedef VECTOR<int,3> TV_INT;typedef VECTOR<T,T_GRID::dimension+2> TV_DIMENSION;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef VECTOR<T,2*T_GRID::dimension> T_FACE_VECTOR;typedef VECTOR<TV,2*T_GRID::dimension> TV_FACE_VECTOR;
    typedef VECTOR<bool,2*T_GRID::dimension> T_FACE_VECTOR_BOOL;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_SCALAR T_LINEAR_INTERPOLATION_SCALAR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;

    using BASE::initial_time;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;using BASE::restart;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::solids_fluids_parameters;
    using BASE::solid_body_collection; using BASE::data_directory;using BASE::parse_args;using BASE::test_number;using BASE::resolution;
    using BASE::solids_evolution; using BASE::stream_type;

    enum SHOCK_TYPE {SPHERICAL,VERTICAL};

    SOLIDS_STANDARD_TESTS<TV> solid_tests;
    FRACTURE_PATTERN<T> fracture_pattern;
    int fp;
    
    int shock_type;
    TV_DIMENSION state_inside,state_outside; // (density,velocity_x,velocity_y,velocity_z,pressure)
    RANGE<TV> inside_shock_box;
    TV shock_center;
    T shock_radius;
    T solid_mass;
    TV sphere_initial_position;
    T sphere_scale;
    int sphere;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    T solid_gravity;
    bool use_solids_gravity;
    bool no_solids;
    bool fracture_walls;
    RIGID_BODY_COLLISION_MANAGER_HASH* collision_manager;

    EOS_SMOOTH_TRANSITION_INCOMPRESSIBLE<EOS_GAMMA<T> >* eos_smooth_transition;
    bool transition_to_incompressible;
    bool incompressible;
    T vorticity_confinement;

    bool use_soot;
    bool use_soot_sourcing;
    bool use_soot_sourcing_from_shock;
    bool use_soot_fuel_combustion;
    bool use_smoke_sourcing;
    CYLINDER<T> smoke_source,soot_source;
    T source_density_value,source_temperature_value,source_soot_value;
    TV source_velocity_value;

    bool read_soot_from_file;
    T_LINEAR_INTERPOLATION_SCALAR soot_interpolation;
    std::string soot_data_dir;

    bool simulate_rigids,simulate_deformable;
    bool use_fixed_farfield_boundary;
    bool run_self_tests;
    bool print_poisson_matrix;
    bool print_index_map;
    bool print_matrix;
    bool print_rhs;
    bool print_each_matrix;
    bool output_iterators;

    ARRAY<DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>*> deformable_objects_to_simulate;

    /***************
    example explanation:
    1. Circular higher density region.
    2. Same as 1 + a two-way coupled circular object.
    3. O'brien's shock hitting a wall
    4. Sphere resting on ground and hit by spherical shock.
    5. Sphere resting on ground and hit by vertical shock.
    6. Sphere resting on ground with wall on right, hit by vertical shock.
    7. Stack of bodies resting on ground with all walls except left, hit by vertical shock.
    8. Stack of squares resting on ground with all walls except left, hit by vertical shock.
    9. Vertical shock hitting a light wall
    10. Vertical shock hitting a heavy wall
    11. Spherical shock hitting a fracture wall.
    12. Shock fractures a wall, interacting with the ground.
    13. Spherical shock in a room with 4 walls.
    14. Spherical shock hitting a light wall
    15. Spherical shock hitting a heavy wall
    16. Spherical shock in a box.
    17. Deformable sphere hit by a spherical shock
    18. Trinity Test
    20. Smoke
    21. Smoke hit by shock
    ***************/

    STANDARD_TESTS(const STREAM_TYPE stream_type,const bool incompressible_input)
        :BASE(stream_type,0,incompressible_input?fluids_parameters.SMOKE:fluids_parameters.COMPRESSIBLE),
        solid_tests(*this,solid_body_collection),rigid_body_collection(solid_body_collection.rigid_body_collection),collision_manager(0),
        incompressible(incompressible_input)
    {
    }
    
    virtual ~STANDARD_TESTS() {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Integer_Argument("-eno_scheme",2,"eno_scheme","eno scheme");
    parse_args->Add_Integer_Argument("-eno_order",2,"eno_order","eno order");
    parse_args->Add_Integer_Argument("-rk_order",3,"rk_order","runge kutta order");
    parse_args->Add_Double_Argument("-cfl",(T).6,"CFL","cfl number");
    parse_args->Add_Option_Argument("-weak_shock","Use stronger shock with temperature ratio 2900:290, p ratio 1000:1");
    parse_args->Add_Option_Argument("-no_solids","ensure no solid is added");
    parse_args->Add_Double_Argument("-mass",(T).0625,"solid_mass","the mass of the solid in the simulation");
    parse_args->Add_Option_Argument("-timesplit","split time stepping into an explicit advection part, and an implicit non-advection part");
    parse_args->Add_Option_Argument("-slip","use slip/spd for coupling");
    parse_args->Add_Option_Argument("-exact","output a fully-explicit sim to (output_dir)_exact");

    parse_args->Add_Option_Argument("-test_system");
    parse_args->Add_Option_Argument("-print_poisson_matrix");
    parse_args->Add_Option_Argument("-print_index_map");
    parse_args->Add_Option_Argument("-print_matrix");
    parse_args->Add_Option_Argument("-print_rhs");
    parse_args->Add_Option_Argument("-print_each_matrix");
    parse_args->Add_Option_Argument("-output_iterators");
    parse_args->Add_Option_Argument("-no_preconditioner");
    parse_args->Add_Option_Argument("-preconditioner");

    parse_args->Add_Option_Argument("-transition_to_incompressible","transition to incompressible in a time window");
    parse_args->Add_Double_Argument("-time_start_transition",(T).5,"time to start transitioning to incompressible flow");
    parse_args->Add_Double_Argument("-time_end_transition",(T).7,"time to end transitioning to incompressible flow");
    parse_args->Add_Double_Argument("-one_over_c_incompressible",(T)0,"one over incompressible sound speed");
    parse_args->Add_Double_Argument("-vorticity_confinement",(T)0,"vorticity_confinement","Vorticity Confinement");
    parse_args->Add_Option_Argument("-no_soot","advect around soot");
    parse_args->Add_Option_Argument("-use_soot_sourcing","source soot");
    parse_args->Add_Option_Argument("-use_soot_sourcing_from_shock","source soot from initial shock place");
    parse_args->Add_Option_Argument("-combustion","source soot");

    parse_args->Add_Double_Argument("-burn_temperature_threshold",(T)600,"temperature above which fuel will burn");
    parse_args->Add_Double_Argument("-burn_rate",(T)100,"rate of fuel burning");
    parse_args->Add_Double_Argument("-calorific_value",(T)54000,"heat generated per fuel burnt");

    parse_args->Add_Option_Argument("-use_smoke_sourcing","source smoke (denisty,temperature,velocity)");

    parse_args->Add_Option_Argument("-read_soot_from_file","read soot from input file");
    parse_args->Add_String_Argument("-soot_data_dir","directory containing grid, soot, velocity files");
    
    parse_args->Add_Option_Argument("-use_fixed_farfield_boundary","use fixed farfield values for outflow boundaries");
    parse_args->Add_Option_Argument("-use_incompressible_gravity","add gravity on incompressible fluid");
    parse_args->Add_Option_Argument("-no_solids_gravity","add gravity on solids");

    parse_args->Add_Integer_Argument("-fp",0,"fracture pattern","fracture pattern");
    parse_args->Add_Option_Argument("-createpattern");
    parse_args->Add_Option_Argument("-fracture_walls");
    
    parse_args->Add_Option_Argument("-faster_frame_rate");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();

    if(parse_args->Get_Option_Value("-createpattern")){
       if(test_number==11) Create_Wall_Pattern();
       exit(0);}

    int eno_scheme=parse_args->Get_Integer_Value("-eno_scheme");
    int eno_order=parse_args->Get_Integer_Value("-eno_order");
    int rk_order=parse_args->Get_Integer_Value("-rk_order");
    T cfl_number=(T)parse_args->Get_Double_Value("-cfl");
    bool strong_shock=!parse_args->Is_Value_Set("-weak_shock");
    no_solids=parse_args->Is_Value_Set("-no_solids");
    solid_mass=(T)parse_args->Get_Double_Value("-mass");
    bool timesplit=parse_args->Is_Value_Set("-timesplit") && !parse_args->Is_Value_Set("-exact");
    bool use_slip=parse_args->Is_Value_Set("-slip");
    run_self_tests=parse_args->Is_Value_Set("-test_system");
    print_poisson_matrix=parse_args->Is_Value_Set("-print_poisson_matrix");
    print_index_map=parse_args->Is_Value_Set("-print_index_map");
    print_matrix=parse_args->Is_Value_Set("-print_matrix");
    print_rhs=parse_args->Is_Value_Set("-print_rhs");
    print_each_matrix=parse_args->Is_Value_Set("-print_each_matrix");
    output_iterators=parse_args->Is_Value_Set("-output_iterators");
    transition_to_incompressible=parse_args->Is_Value_Set("-transition_to_incompressible");
    T time_start_transition=(T)parse_args->Get_Double_Value("-time_start_transition");
    T time_end_transition=(T)parse_args->Get_Double_Value("-time_end_transition");
    T one_over_c_incompressible=(T)parse_args->Get_Double_Value("-one_over_c_incompressible");
    T vorticity_confinement=(T)parse_args->Get_Double_Value("-vorticity_confinement");
    use_soot=!parse_args->Is_Value_Set("-no_soot");
    use_soot_sourcing=parse_args->Is_Value_Set("-use_soot_sourcing");
    use_soot_fuel_combustion=parse_args->Is_Value_Set("-combustion");
    T burn_temperature_threshold=parse_args->Get_Double_Value("-burn_temperature_threshold");
    T burn_rate=parse_args->Get_Double_Value("-burn_rate");
    T soot_fuel_calorific_value=parse_args->Get_Double_Value("-calorific_value");
    use_soot_sourcing_from_shock=parse_args->Is_Value_Set("-use_soot_sourcing_from_shock");
    use_smoke_sourcing=parse_args->Is_Value_Set("-use_smoke_sourcing");
    read_soot_from_file=parse_args->Is_Value_Set("-read_soot_from_file");
    soot_data_dir=parse_args->Get_String_Value("-soot_data_dir");
    use_fixed_farfield_boundary=parse_args->Is_Value_Set("-use_fixed_farfield_boundary");
    bool use_incompressible_gravity=parse_args->Is_Value_Set("-use_incompressible_gravity");
    use_solids_gravity=!parse_args->Is_Value_Set("-no_solids_gravity");
    fracture_walls=parse_args->Is_Value_Set("-fracture_walls");
    fp=parse_args->Get_Integer_Value("-fp");
    bool faster_frame_rate=parse_args->Is_Value_Set("-faster_frame_rate");

    //grid
    int cells=resolution;
    T grid_size=(T)1.;
    if(test_number==7||test_number==17){
        fluids_parameters.grid->Initialize(3*cells+1,2*cells+1,2*cells+1,RANGE<TV>(TV((T)-4.5,(T)-3,(T)-3),TV((T)4.5,(T)3,(T)3))*grid_size);}
    else if(test_number==3){
        fluids_parameters.grid->Initialize(cells+1,cells+1,cells+1,RANGE<TV>(TV((T)-10,(T)-10,(T)-10),TV((T)10,(T)10,(T)10))*grid_size);}
    else if(test_number==9||test_number==10||test_number==11){
        fluids_parameters.grid->Initialize(3*cells+1,2*cells+1,2*cells+1,RANGE<TV>(TV((T)-1.5,(T)-1,(T)-1),TV((T)1.5,(T)1,(T)1))*grid_size);}
    else if(test_number==12){
        fluids_parameters.grid->Initialize(cells+1,cells+1,cells+1,RANGE<TV>(TV((T)-1,(T)0,(T)-1),TV((T)1,(T)2,(T)1))*grid_size);}
    else if(test_number==13){
        fluids_parameters.grid->Initialize(2*cells+1,3*cells+1,2*cells+1,RANGE<TV>(TV((T)-3,(T)0,(T)-3),TV((T)5,(T)12,(T)5))*grid_size);}
    else if(test_number==16){
        fluids_parameters.grid->Initialize(cells+1,cells+1,cells+1,RANGE<TV>(TV((T)-3,(T)0,(T)-3),TV((T)5,(T)8,(T)5))*grid_size);}
    else if(test_number==14||test_number==15){
        fluids_parameters.grid->Initialize(3*cells+1,2*cells+1,2*cells+1,RANGE<TV>(TV((T)-120,(T)-80,(T)-80),TV((T)120,(T)80,(T)80))*grid_size);}
    else if(test_number==18){
        fluids_parameters.grid->Initialize(2*cells+1,cells+1,2*cells+1,RANGE<TV>(TV((T)-100,(T)0,(T)-100),TV((T)100,(T)100,(T)100))*grid_size);}
    else if(test_number==19){
        fluids_parameters.grid->Initialize(5*cells+1,2*cells+1,2*cells+1,RANGE<TV>(TV((T)-8,(T)-10,(T)-10),TV((T)42,(T)10,(T)10))*grid_size);}
    else if(test_number==20){
        fluids_parameters.grid->Initialize(cells+1,2*cells+1,cells+1,RANGE<TV>(TV((T)-5,(T)-7.5,(T)-5),TV((T)5,(T)7.5,(T)5))*grid_size);}
    else if(test_number==21){
        fluids_parameters.grid->Initialize(2*cells+1,cells+1,cells+1,RANGE<TV>(TV((T)-1,(T)0,(T)-.25),TV((T)1,(T)1,(T).75))*grid_size);}
    else{
        fluids_parameters.grid->Initialize(cells+1,cells+1,cells+1,RANGE<TV>(TV((T)-1,(T)-1,(T)-1),TV((T)1,(T)1,(T)1))*grid_size);}

    *fluids_parameters.grid=fluids_parameters.grid->Get_MAC_Grid();
    fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[1][2]=false;fluids_parameters.domain_walls[2][1]=false;
    fluids_parameters.domain_walls[2][2]=false;fluids_parameters.domain_walls[3][1]=false;fluids_parameters.domain_walls[3][2]=false;
    if(test_number==3 || test_number==4 || test_number==5) fluids_parameters.domain_walls[2][1]=true;
    if(test_number==6){
        fluids_parameters.domain_walls[1][2]=true;
        fluids_parameters.domain_walls[2][1]=true;}
    if(test_number==7||test_number==8||test_number==9||test_number==10||test_number==11||test_number==14||test_number==15||test_number==17){
        fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[1][2]=true;
        fluids_parameters.domain_walls[2][1]=true;fluids_parameters.domain_walls[2][2]=true;
        fluids_parameters.domain_walls[3][1]=true;fluids_parameters.domain_walls[3][2]=true;}
    if(test_number==13||test_number==16||test_number==18||test_number==19){
        fluids_parameters.domain_walls[2][1]=true;}
    if(test_number==19){
        fluids_parameters.domain_walls[1][2]=true;fluids_parameters.domain_walls[2][1]=true;}
    if(test_number==20){
        fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[1][2]=true;
        fluids_parameters.domain_walls[2][1]=true;fluids_parameters.domain_walls[2][2]=false;
        fluids_parameters.domain_walls[3][1]=true;fluids_parameters.domain_walls[3][2]=true;}
    if(test_number==21){
        fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[1][2]=false;
        fluids_parameters.domain_walls[2][1]=true;fluids_parameters.domain_walls[2][2]=true;
        fluids_parameters.domain_walls[3][1]=true;fluids_parameters.domain_walls[3][2]=true;}


    //time
    initial_time=(T)0.;last_frame=1000;frame_rate=(T)1e4;

    switch(test_number){
        case 3: last_frame=1;
                frame_rate=(T)75;
                break;
        case 7: frame_rate=(T)5e4;
                break;
        case 14:
        case 15: frame_rate=(T)1e4;
                 break;
        case 17: frame_rate=(T)1e5;
                 break;
        case 18:last_frame=200;frame_rate=(T)2.5e4;
                break;
        case 19:frame_rate=(T)1e4;
                break;
        case 20:frame_rate=32;
                break;
        case 21:frame_rate=(T)6e5;
                break;
        default:frame_rate=(T)1e4;last_frame=1000;}

    if(faster_frame_rate) frame_rate=96;
    fluids_parameters.cfl=cfl_number;
    //custom stuff . . . 
    if(transition_to_incompressible){
        eos_smooth_transition=new EOS_SMOOTH_TRANSITION_INCOMPRESSIBLE<EOS_GAMMA<T> >(time_start_transition,time_end_transition,one_over_c_incompressible,true,3);
        fluids_parameters.compressible_eos=eos_smooth_transition;}
    else fluids_parameters.compressible_eos = new EOS_GAMMA<T>;
    if(eno_scheme==1) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,false,false);
    else if(eno_scheme==2) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,true,false);
    else fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+2>(true,true,true);
    fluids_parameters.compressible_spatial_order=eno_order;
    fluids_parameters.compressible_conservation_method->Save_Fluxes();
    //fluids_parameters.compressible_conservation_method->Scale_Outgoing_Fluxes_To_Clamp_Variable(true,1,(T)1e-5);
    fluids_parameters.compressible_rungekutta_order=rk_order;
    fluids_parameters.compressible_timesplit=timesplit;
    fluids_parameters.use_slip=use_slip;

    fluids_parameters.use_preconditioner_for_slip_system=true;
    if(parse_args->Is_Value_Set("-preconditioner")) fluids_parameters.use_preconditioner_for_slip_system=true;
    if(parse_args->Is_Value_Set("-no_preconditioner")) fluids_parameters.use_preconditioner_for_slip_system=false;

    if(use_soot){
        fluids_parameters.use_soot=true;
        if(use_soot_fuel_combustion){
            fluids_parameters.use_soot_fuel_combustion=true;
            fluids_parameters.burn_temperature_threshold=burn_temperature_threshold;
            fluids_parameters.burn_rate=burn_rate;
            fluids_parameters.soot_fuel_calorific_value=soot_fuel_calorific_value;}
        else fluids_parameters.use_soot_fuel_combustion=false;
        fluids_parameters.use_fixed_soot_boundary=true;
        fluids_parameters.ambient_soot=(T)0;
        fluids_parameters.soot_boundary=new BOUNDARY_REFLECTION_ATTENUATION<T_GRID,T>(VECTOR_UTILITIES::Complement(fluids_parameters.domain_walls),fluids_parameters.ambient_soot,(T)1);}

    if(incompressible){
        if(vorticity_confinement>0){
            fluids_parameters.use_vorticity_confinement=true;
            fluids_parameters.confinement_parameter=vorticity_confinement*fluids_parameters.grid->Minimum_Edge_Length();}
        else fluids_parameters.use_vorticity_confinement=false;
        solids_fluids_parameters.use_leakproof_solve=false;
        fluids_parameters.use_body_force=false;
        fluids_parameters.density=(T)1; // not used
        fluids_parameters.use_density=true;
        fluids_parameters.use_temperature=true;
        fluids_parameters.use_fixed_density_boundary=true;
        fluids_parameters.use_fixed_temperature_boundary=true;
        if(use_incompressible_gravity) fluids_parameters.gravity=(T)9.8;
        else fluids_parameters.gravity=(T)0;}

    // setup solids
    bool simulate_rigids=(test_number==2||test_number==4||test_number==5||test_number==6||test_number==7||test_number==8||test_number==9||test_number==10||test_number==11||test_number==12||test_number==13||test_number==14||test_number==15||test_number==16||test_number==21);
    if(use_slip) simulate_rigids=true;
    bool simulate_deformable=(test_number==17);

    if(simulate_rigids||simulate_deformable){
        solids_fluids_parameters.use_leakproof_solve=false;
        fluids_parameters.solid_affects_fluid=true;
        fluids_parameters.fluid_affects_solid=true;
        solids_parameters.use_post_cg_constraints=true;
        solids_parameters.triangle_collision_parameters.perform_self_collision=false;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.rigid_body_collision_parameters.use_push_out=false;
        solids_parameters.rigid_body_collision_parameters.use_legacy_push_out=true;
        solids_parameters.rigid_body_collision_parameters.use_fracture_particle_optimization=false;
        solids_parameters.rigid_body_evolution_parameters.maximum_rigid_body_time_step_fraction=(T)1;
        solids_parameters.use_trapezoidal_rule_for_velocities=false;
        solids_parameters.verbose_dt=true;
        solid_body_collection.print_residuals=false;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-16;
        solids_parameters.implicit_solve_parameters.cg_restart_iterations=200;
        solids_parameters.implicit_solve_parameters.cg_iterations=1000;
    
        solids_parameters.use_rigid_deformable_contact=true;
        solids_parameters.rigid_body_collision_parameters.enforce_rigid_rigid_contact_in_cg=true;

        if(solids_parameters.rigid_body_collision_parameters.enforce_rigid_rigid_contact_in_cg){
            solids_parameters.implicit_solve_parameters.cg_projection_iterations=20;}
        else solids_parameters.implicit_solve_parameters.cg_projection_iterations=0;
        solids_parameters.rigid_body_collision_parameters.contact_iterations=20;
        solids_parameters.rigid_body_collision_parameters.contact_project_iterations=20;
        solids_parameters.rigid_body_collision_parameters.perform_contact=true;
        solids_parameters.rigid_body_collision_parameters.use_shock_propagation=true;

        if(simulate_deformable) solid_body_collection.deformable_body_collection.simulate=true;
        
        sphere_initial_position=TV((T).65,(T).11,(T)0);
        sphere_scale=(T).25;
        solid_gravity=(T)9.8;}

    if(fracture_walls) solids_fluids_parameters.use_fluid_rigid_fracture=true;

    if(test_number==4){sphere_initial_position=TV((T)-1+sphere_scale*(T)2,-(T)1+sphere_scale,(T)0);}
    if(test_number==5){sphere_initial_position=TV((T).2+sphere_scale*(T)2,-(T)1+sphere_scale,(T)0);}
    if(test_number==6){sphere_initial_position=TV(-sphere_scale*(T)2,-(T)1+sphere_scale,(T)0);}
 
    // shock type
    if(test_number==1||test_number==2||test_number==3||test_number==4||test_number==11||test_number==12||test_number==13||test_number==14||test_number==15||test_number==16||test_number==17) shock_type=SPHERICAL;
    else if(test_number==5||test_number==6||test_number==7||test_number==8||test_number==9||test_number==10||test_number==19||test_number==21) shock_type=VERTICAL;

    // spherical shock positioning
    T_GRID& grid=(fluids_parameters.mpi_grid)?fluids_parameters.mpi_grid->global_grid:*fluids_parameters.grid;
    shock_radius=(T).6;
    shock_center=TV();
    if(test_number==3){
        shock_radius=(T).5;
        shock_center.y-=6.5;
        strong_shock=true;}
    if(test_number==4) shock_radius=(T).1;
    if(test_number==11)shock_radius=.4;
    if(test_number==12){
        shock_radius=.4;
        shock_center.y=(T)1;}
    if(test_number==13||test_number==16){
        shock_center=TV(1,sphere_scale,(T)1);
        shock_radius=.4;}
    if(test_number==14||test_number==15){
        shock_radius=.3*grid.Domain().max_corner.y;
        shock_center=TV(grid.Domain().max_corner.x*.5-shock_radius*2,grid.Domain().min_corner.y+shock_radius*1.2,0);}
    if(test_number==17){
        shock_radius=(T)grid.Domain().max_corner.x*.15;
        shock_center=TV(-shock_radius*(1.2),-shock_radius,0);}
    if(test_number==18){
        shock_radius=(T)6.5;
        shock_center=TV(0,10,0);}

    // vertical shock positioning
    if(test_number==5){inside_shock_box=grid.Domain()-TV(1,0,0);}
    if(test_number==6){inside_shock_box=grid.Domain()-TV(1.7,0,0);}
    if(test_number==7||test_number==21){
        T x_delta=grid.Domain().max_corner.x*(T)1.2;
        inside_shock_box=grid.Domain()-TV(x_delta,0,0);}
    if(test_number==8){inside_shock_box=grid.Domain()-TV(1.2,0,0);}
    if(test_number==9||test_number==10){inside_shock_box=grid.Domain()-TV(1.5,0,0);}
    if(test_number==19) inside_shock_box=RANGE<TV>(TV(36,-2,-2),TV(42,2,2));
    
    // shock strength
    state_inside=TV_DIMENSION((T)1,(T)0,(T)0,(T)0,(T)1);
    state_outside=TV_DIMENSION((T).125,(T)0,(T)0,(T)0,(T).1);
    if(strong_shock){
        T p_atm=(T)101325.;
        T temperature_inside=(T)2900,temperature_outside=(T)290;
        T p_inside=(T)1000*p_atm,p_outside=p_atm;

        T rho_inside=fluids_parameters.compressible_eos->rho_From_p_And_T(p_inside,temperature_inside);
        T rho_outside=fluids_parameters.compressible_eos->rho_From_p_And_T(p_outside,temperature_outside);

        state_inside=TV_DIMENSION(rho_inside,(T)0,(T)0,(T)0,p_inside);
        state_outside=TV_DIMENSION(rho_outside,(T)0,(T)0,(T)0,p_outside);}

    if(test_number==18){
        T p_atm=(T)101325.;
        T temperature_inside=(T)2.62497e8,temperature_outside=(T)290;
        T p_inside=(T)9.41831e10,p_outside=p_atm;

        T rho_inside=fluids_parameters.compressible_eos->rho_From_p_And_T(p_inside,temperature_inside);
        T rho_outside=fluids_parameters.compressible_eos->rho_From_p_And_T(p_outside,temperature_outside);

        state_inside=TV_DIMENSION(rho_inside,(T)0,(T)0,(T)0,p_inside);
        state_outside=TV_DIMENSION(rho_outside,(T)0,(T)0,(T)0,p_outside);}


    if(incompressible){
         // Set ambient density and temperature
         EOS<T>& eos=*fluids_parameters.compressible_eos;
         T rho_outside=state_outside(1);
         T p_outside=state_outside(5);
         T e_outside=eos.e_From_p_And_rho(p_outside,rho_outside);
         T temperature_outside=eos.T(rho_outside,e_outside);
         fluids_parameters.ambient_density=rho_outside;
         fluids_parameters.ambient_temperature=temperature_outside;
         T linear_attenuation;
         if(use_fixed_farfield_boundary) linear_attenuation=(T)1;
         else linear_attenuation=(T).1;
         fluids_parameters.density_boundary=new BOUNDARY_REFLECTION_ATTENUATION<T_GRID,T>(VECTOR_UTILITIES::Complement(fluids_parameters.domain_walls),fluids_parameters.ambient_density,linear_attenuation);
         fluids_parameters.temperature_boundary=new BOUNDARY_REFLECTION_ATTENUATION<T_GRID,T>(VECTOR_UTILITIES::Complement(fluids_parameters.domain_walls),fluids_parameters.ambient_temperature,linear_attenuation);}

     if(use_smoke_sourcing){
         if(test_number==20){
            smoke_source.radius=(T)1;
            smoke_source.Set_Endpoints(TV((T)0,grid.Domain().min_corner.y,0),TV((T)0,grid.Domain().min_corner.y+smoke_source.radius,0));

            source_density_value=(T)state_outside(1)*.01;
            source_temperature_value=(T)1000;
            source_velocity_value=TV(0,0,0);}}
    if(use_soot_sourcing){
        if(test_number==20){
            soot_source=smoke_source;
            source_soot_value=(T)1;}
        if(test_number==21){
            source_soot_value=(T)1;
            soot_source.radius=T(.05);
            VECTOR<T,3> endpoint1=VECTOR<T,3>::Constant_Vector(0.25),endpoint2=VECTOR<T,3>::Constant_Vector(0.25);
            endpoint1(2)=T(0);endpoint2(2)=T(0.05);
            soot_source.Set_Endpoints(endpoint1,endpoint2);}}

    // output directory
    if(timesplit) output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d__Resolution_%d_%d_%d_semiimplicit",
        test_number,(fluids_parameters.grid->counts.x),(fluids_parameters.grid->counts.y),(fluids_parameters.grid->counts.z));
    else output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d__Resolution_%d_%d_%d_explicit",test_number,
        (fluids_parameters.grid->counts.x),(fluids_parameters.grid->counts.y),(fluids_parameters.grid->counts.z));
    if(eno_scheme==2) output_directory+="_density_weighted";
    else if(eno_scheme==3) output_directory+="_velocity_weighted";
    if(use_slip) output_directory+="_slip";
    if(transition_to_incompressible) output_directory+="_transition_incompressible";
    if(use_soot) output_directory+="_soot";
    if(use_fixed_farfield_boundary) output_directory+="_fixedFF";
    if(strong_shock) output_directory+="_strong";
    if(fracture_walls) output_directory+="_fracture";
    output_directory+=STRING_UTILITIES::string_sprintf("_mass_%f",solid_mass);
    output_directory+=STRING_UTILITIES::string_sprintf("_fp_%d",fp);
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();
    if(!(test_number==13 || test_number==19))
        fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection.rigid_geometry_collection);
}
//#####################################################################
// Function Intialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    if(incompressible){
        fluids_parameters.Use_Fluid_Coupling_Defaults();
        return;}

    //set custom boundary
    VECTOR<VECTOR<bool,2>,T_GRID::dimension> valid_wall;
    for(int axis=1;axis<=T_GRID::dimension;axis++) for(int axis_side=1;axis_side<=2;axis_side++)
        valid_wall[axis][axis_side]=(fluids_parameters.mpi_grid?!fluids_parameters.mpi_grid->Neighbor(axis,axis_side):true) && !fluids_parameters.domain_walls[axis][axis_side];

    TV far_field_velocity=TV(state_outside(2),state_outside(3),state_outside(4));
    if(use_fixed_farfield_boundary){
        fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<T_GRID>(fluids_parameters.euler,
            T_FACE_VECTOR(state_outside(1),state_outside(1),state_outside(1),state_outside(1),state_outside(1),state_outside(1)),
            T_FACE_VECTOR(state_outside(5),state_outside(5),state_outside(5),state_outside(5),state_outside(5),state_outside(5)),
            TV_FACE_VECTOR(far_field_velocity,far_field_velocity,far_field_velocity,far_field_velocity,far_field_velocity,far_field_velocity),
            (T).5,valid_wall,true,T_FACE_VECTOR(1,1,1,1,1,1),T_FACE_VECTOR_BOOL(true,true,true,true,true,true));}
    else{
        fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<T_GRID>(fluids_parameters.euler,
            T_FACE_VECTOR(state_outside(1),state_outside(1),state_outside(1),state_outside(1),state_outside(1),state_outside(1)),
            T_FACE_VECTOR(state_outside(5),state_outside(5),state_outside(5),state_outside(5),state_outside(5),state_outside(5)),
            TV_FACE_VECTOR(far_field_velocity,far_field_velocity,far_field_velocity,far_field_velocity,far_field_velocity,far_field_velocity),
            (T).5,valid_wall);}
}
//#####################################################################
// Function Inside_Shock
//#####################################################################
bool Inside_Shock(const TV location)
{
    if(shock_type==SPHERICAL) return (location-shock_center).Magnitude()<shock_radius;
    else if(shock_type==VERTICAL) return inside_shock_box.Lazy_Inside(location);
    else PHYSBAM_FATAL_ERROR("Undefined shock type");
}
//#####################################################################
// Function Initialize_Euler_State
//#####################################################################
void Initialize_Euler_State()
{
    if(incompressible) return;

    T_GRID& grid=fluids_parameters.euler->grid;
    ARRAY<VECTOR<T,5> ,VECTOR<int,3> >& U=fluids_parameters.euler->U;
    EOS<T> *eos = fluids_parameters.euler->eos;

    if(transition_to_incompressible) fluids_parameters.euler->euler_projection.use_neumann_condition_for_outflow_boundaries=false;

    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(fluids_parameters.euler->grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        T rho,u_vel,v_vel,w_vel,p;
        if(Inside_Shock(grid.X(cell_index))){
            rho=state_inside(1);u_vel=state_inside(2);v_vel=state_inside(3);w_vel=state_inside(4);p=state_inside(5);}
        else{rho=state_outside(1);u_vel=state_outside(2);v_vel=state_outside(3);w_vel=state_outside(4);p=state_outside(5);}

        U(cell_index)(1)=rho;U(cell_index)(2)=rho*u_vel;U(cell_index)(3)=rho*v_vel;U(cell_index)(4)=rho*w_vel;
        U(cell_index)(5)=rho*(eos->e_From_p_And_rho(p,rho)+(sqr(u_vel)+sqr(v_vel)+sqr(w_vel))/(T)2.);}

    if(read_soot_from_file) Read_Soot_Velocities(); // read in velocities in soot domain
}
//#####################################################################
// Function Read_Soot_Velocities
//#####################################################################
void Read_Soot_Velocities()
{
    T_GRID& grid=*fluids_parameters.grid;
    ARRAY<VECTOR<T,5> ,VECTOR<int,3> >& U=fluids_parameters.euler->U;
    std::string soot_grid_file=soot_data_dir+"/grid";
    std::string soot_velocity_file=soot_data_dir+"/mac_velocities";
    T_GRID soot_grid;
    T_ARRAYS_SCALAR soot_values;
    T_FACE_ARRAYS_SCALAR soot_mac_velocities;
    FILE_UTILITIES::Read_From_File(stream_type,soot_grid_file,soot_grid);
    FILE_UTILITIES::Read_From_File(stream_type,soot_velocity_file,soot_mac_velocities);
    T soot_dx_over_2=soot_grid.dX.Max()*(T).5;
    RANGE<TV> soot_domain=soot_grid.Domain();
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        TV location=iterator.Location();
        if(soot_domain.Inside(location,soot_dx_over_2)){
            T velocity;
            for(int axis=1;axis<=T_GRID::dimension;axis++){
                velocity=soot_interpolation.Clamped_To_Array_Face_Component(axis,soot_grid,FACE_LOOKUP_UNIFORM<T_GRID>(soot_mac_velocities),location);
                U(cell_index)(axis+1)=U(cell_index)(1)*velocity;}}}
}
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Adjust_Density_And_Temperature_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    T_GRID& grid=*fluids_parameters.grid;
    if(!incompressible) PHYSBAM_FATAL_ERROR("this shouldn't be called in compressible case");
    if(!use_smoke_sourcing && time>0) return;

    if(!use_smoke_sourcing){ // use shock data to initialize density and temperature
        EOS<T>& eos=*fluids_parameters.compressible_eos;
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT cell_index=iterator.Cell_Index();
            T rho,temperature,p,e;
            if(Inside_Shock(grid.X(cell_index))){
                rho=state_inside(1);
                p=state_inside(5);
                e=eos.e_From_p_And_rho(p,rho);
                temperature=eos.T(rho,e);}
            else{
                rho=state_outside(1);
                p=state_outside(5);
                e=eos.e_From_p_And_rho(p,rho);
                temperature=eos.T(rho,e);}
            fluids_parameters.density_container.density(cell_index)=rho;
            fluids_parameters.temperature_container.temperature(cell_index)=temperature;}}
    else{
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT cell_index=iterator.Cell_Index();
            if(smoke_source.Lazy_Inside(grid.X(cell_index))){
                fluids_parameters.density_container.density(cell_index)=source_density_value;
                fluids_parameters.temperature_container.temperature(cell_index)=source_temperature_value;}}}
}
//#####################################################################
// Function Clear_Inside_Solid_Soot
//#####################################################################
void Clear_Inside_Solid_Soot()
{
    T_GRID& grid=*fluids_parameters.grid;
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>& collision_bodies_affecting_fluid=*fluids_parameters.collision_bodies_affecting_fluid;

    for(COLLISION_GEOMETRY_ID id(1);id<=collision_bodies_affecting_fluid.collision_geometry_collection.bodies.m;id++)
        if(collision_bodies_affecting_fluid.collision_geometry_collection.Is_Active(id)){
            COLLISION_GEOMETRY<TV>& collision_body=*(collision_bodies_affecting_fluid.collision_geometry_collection.bodies(id));
            collision_body.Update_Bounding_Box();
            RANGE<TV_INT> grid_clamped_bounding_box=grid.Clamp_To_Cell(collision_body.Axis_Aligned_Bounding_Box(),0);
            for(CELL_ITERATOR iterator(grid,grid_clamped_bounding_box);iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                T phi_value=collision_body.Implicit_Geometry_Extended_Value(iterator.Location());
                if(phi_value<(T)0){
                    fluids_parameters.soot_container.density(cell_index)=0;
                    if(use_soot_fuel_combustion) fluids_parameters.soot_fuel_container.density(cell_index)=0;}}}
}
//#####################################################################
// Function Adjust_Soot_With_Sources
//#####################################################################
void Adjust_Soot_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    T_GRID& grid=*fluids_parameters.grid;
    if(!use_soot) PHYSBAM_FATAL_ERROR("this shouldn't be called in non use_soot case");

    if(time>0) Clear_Inside_Solid_Soot();

    if(!use_soot_sourcing && !use_soot_sourcing_from_shock && time>0) return;

    if(read_soot_from_file && time==(T)0){
        std::string soot_grid_file=soot_data_dir+"/grid";
        std::string soot_file=soot_data_dir+"/density";
        T_GRID soot_grid;
        T_ARRAYS_SCALAR soot_values;
        FILE_UTILITIES::Read_From_File(stream_type,soot_grid_file,soot_grid);
        FILE_UTILITIES::Read_From_File(stream_type,soot_file,soot_values);
        T soot_dx_over_2=soot_grid.dX.Max()*(T).5;
        RANGE<TV> soot_domain=soot_grid.Domain();
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            TV location=iterator.Location();
            if(soot_domain.Inside(location,soot_dx_over_2)){
                fluids_parameters.soot_container.density(cell_index)=
                    soot_interpolation.Clamped_To_Array(soot_grid,soot_values,location);}
            else fluids_parameters.soot_container.density(cell_index)=(T)0;}}
    else if(use_soot_sourcing){
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT cell_index=iterator.Cell_Index();
            if(soot_source.Lazy_Inside(grid.X(cell_index))){
                fluids_parameters.soot_container.density(cell_index)=source_soot_value;}}}
    else if(use_soot_sourcing_from_shock){ // use shock region to source soot
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT cell_index=iterator.Cell_Index();
            if(Inside_Shock(grid.X(cell_index))) fluids_parameters.soot_container.density(cell_index)=(T)1;}}
    else{ // use shock data to initialize soot
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT cell_index=iterator.Cell_Index();
            T soot,fuel;
            if(Inside_Shock(grid.X(cell_index))){
                if(use_soot_fuel_combustion){soot=(T).5;fuel=(T).3;}
                else {soot=(T)1;fuel=(T)0;}}
            else{soot=(T)0;fuel=(T)0;}
            fluids_parameters.soot_container.density(cell_index)=soot;
            if(use_soot_fuel_combustion) fluids_parameters.soot_fuel_container.density(cell_index)=fuel;}}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(T_FACE_ARRAYS_SCALAR& face_velocities,T_FACE_ARRAYS_BOOL& psi_N,const T time) PHYSBAM_OVERRIDE
{
    if(!use_smoke_sourcing) return;
    if(!incompressible) PHYSBAM_FATAL_ERROR("this shouldn't be called in compressible case");
    for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
        if(smoke_source.Lazy_Inside(iterator.Location())){
            int axis=iterator.Axis();
            if(source_velocity_value(axis)){
                psi_N.Component(axis)(iterator.Face_Index())=true;
                face_velocities.Component(axis)(iterator.Face_Index())=source_velocity_value(axis);}}}
}
//#####################################################################
// Function Set_Rigid_Body_Parameters
//#####################################################################
void Set_Rigid_Body_Parameters(int rigid_body_index,const std::string& name,const TV position,const T mass,const T mu)
{
    rigid_body_collection.Rigid_Body(rigid_body_index).name=name;
    rigid_body_collection.Rigid_Body(rigid_body_index).Set_Coefficient_Of_Restitution((T)1);
    rigid_body_collection.Rigid_Body(rigid_body_index).Set_Coefficient_Of_Friction((T)mu);
    rigid_body_collection.rigid_body_particle.X(rigid_body_index)=position;
    rigid_body_collection.Rigid_Body(rigid_body_index).Set_Mass(mass);
    rigid_body_collection.rigid_body_particle.kinematic(rigid_body_index)=false;
    rigid_body_collection.Rigid_Body(rigid_body_index).simplicial_object->mesh.Initialize_Adjacent_Elements();
}
//#####################################################################
// Function Add_Sphere
//#####################################################################
void Add_Sphere()
{
    sphere=rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/sphere_refined",sphere_scale,true,true,false);
    Set_Rigid_Body_Parameters(sphere,"sphere",sphere_initial_position,solid_mass,.5);
}
//#####################################################################
// Function Add_Ground
//#####################################################################
void Add_Ground(bool is_static)
{
    T_GRID& grid=(fluids_parameters.mpi_grid)?fluids_parameters.mpi_grid->global_grid:*fluids_parameters.grid;
    T scale=(T)grid.Domain().max_corner.x*10;
    RIGID_BODY<TV>& ground=solid_tests.Add_Analytic_Box(TV::All_Ones_Vector()*scale);
    if(test_number==12) Set_Rigid_Body_Parameters(ground.particle_index,"ground",TV(0,-scale*(T).5,0),(T)1e10,.5);
    else if(test_number==13||test_number==16) Set_Rigid_Body_Parameters(ground.particle_index,"ground",TV(1,-scale*(T).5,1),(T)1e10,.5);
    else Set_Rigid_Body_Parameters(ground.particle_index,"ground",TV(0,grid.Domain().min_corner.y-scale*(T).5,0),(T)1e10,.5);
    if(is_static) ground.is_static=true;
}
//#####################################################################
// Function Add_Wall
//#####################################################################
void Add_Wall()
{
    T mu=1;
    T x_scale=(T)2,y_scale=(T)6,z_scale=(T)30;
    RIGID_BODY<TV>& wall=solid_tests.Add_Analytic_Box(TV(x_scale,y_scale,z_scale));
    Set_Rigid_Body_Parameters(wall.particle_index,"wall",TV(4+x_scale*(T).5,y_scale*(T).5-10,0),(T)1e10,mu);
}
//#####################################################################
// Function Add_Destructive_Wall
//#####################################################################
void Add_Destructive_Wall()
{
    if(fracture_walls) FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/fracture_pattern-%d",data_directory.c_str(),fp),fracture_pattern);

    TV edge_lengths((T)1.5,(T).1,(T)1);
    TV_INT dimensions(31,3,21);BOX<TV> box((T)-.5*edge_lengths,(T).5*edge_lengths);TV dx=edge_lengths/TV(dimensions-1);
    TV_INT levelset_resolution(151,11,101);TV levelset_dx=edge_lengths/TV(levelset_resolution-1);int ghost_cells=2;
    RIGID_BODY<TV>* rigid_body=new RIGID_BODY<TV>(rigid_body_collection,true);

    rigid_body->Add_Structure(*TESSELLATION::Generate_Triangles(BOX<TV>((T)-.5*edge_lengths,(T).5*edge_lengths)));
    LEVELSET_IMPLICIT_OBJECT<TV>& implicit_object=*LEVELSET_IMPLICIT_OBJECT<TV>::Create();
    implicit_object.levelset.grid.Initialize(levelset_resolution+2*ghost_cells,box.Thickened(levelset_dx.x*ghost_cells),false);
    implicit_object.levelset.phi.Resize(implicit_object.levelset.grid.Domain_Indices());implicit_object.levelset.phi.Fill(FLT_MAX);
    implicit_object.Update_Box();implicit_object.Update_Minimum_Cell_Size();
    for(NODE_ITERATOR iterator(implicit_object.levelset.grid);iterator.Valid();iterator.Next())
        implicit_object.levelset.phi(iterator.index)=box.Signed_Distance(iterator.Location());
    rigid_body->Add_Structure(implicit_object);
    if(test_number==11){
        rigid_body->X()=TV(.8,.2,0);
        rigid_body->Rotation()=ROTATION<TV>((T)-pi/2,TV(0,0,(T)1))*ROTATION<TV>((T)pi/2,TV(0,(T)1,0));}
    else if(test_number==12)
        rigid_body->X()=TV(0,.1,0);
    rigid_body_collection.Add_Rigid_Body_And_Geometry(rigid_body);
    rigid_body->fracture_threshold=(T).005;
}
//#####################################################################
// Function Add_Room
//#####################################################################
void Add_Room()
{
    if(fracture_walls) FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/fracture_pattern-%d",data_directory.c_str(),fp),fracture_pattern);

    for(int i=1;i<=4;++i){
        TV edge_lengths((T)1,(T).2,(T)1);
        TV_INT dimensions(31,3,21);BOX<TV> box((T)-.5*edge_lengths,(T).5*edge_lengths);TV dx=edge_lengths/TV(dimensions-1);
        TV_INT levelset_resolution(151,11,101);TV levelset_dx=edge_lengths/TV(levelset_resolution-1);int ghost_cells=2;
        RIGID_BODY<TV>* rigid_body=new RIGID_BODY<TV>(rigid_body_collection,true);

        rigid_body->Add_Structure(*TESSELLATION::Generate_Triangles(BOX<TV>((T)-.5*edge_lengths,(T).5*edge_lengths)));
        rigid_body->simplicial_object->mesh.Initialize_Adjacent_Elements();
        LEVELSET_IMPLICIT_OBJECT<TV>& implicit_object=*LEVELSET_IMPLICIT_OBJECT<TV>::Create();
        implicit_object.levelset.grid.Initialize(levelset_resolution+2*ghost_cells,box.Thickened(levelset_dx.x*ghost_cells),false);
        implicit_object.levelset.phi.Resize(implicit_object.levelset.grid.Domain_Indices());implicit_object.levelset.phi.Fill(FLT_MAX);
        implicit_object.Update_Box();implicit_object.Update_Minimum_Cell_Size();
        for(NODE_ITERATOR iterator(implicit_object.levelset.grid);iterator.Valid();iterator.Next())
            implicit_object.levelset.phi(iterator.index)=box.Signed_Distance(iterator.Location());
        rigid_body->Add_Structure(implicit_object);
        ROTATION<TV> rotation((T)i*(T)pi/2,TV(0,1,0));
        rigid_body->Rotation()=ROTATION<TV>((T)i*(T)pi/2,TV(0,(T)1,0))*ROTATION<TV>((T)-pi/2,TV(0,0,(T)1))*ROTATION<TV>((T)pi/2,TV(0,(T)1,0));
        rigid_body_collection.Add_Rigid_Body_And_Geometry(rigid_body);
        if((i==3) && fracture_walls){
            rigid_body->fracture_threshold=(T)18;
            Set_Rigid_Body_Parameters(rigid_body->particle_index,STRING_UTILITIES::string_sprintf("wall %d",i).c_str(),TV(1,0,1)+rotation.Rotate(TV((T).8,(T).5,0)),
                                      (T)1.5e1,(T).5);}
        else{
            Set_Rigid_Body_Parameters(rigid_body->particle_index,STRING_UTILITIES::string_sprintf("wall %d",i).c_str(),TV(1,0,1)+rotation.Rotate(TV((T).8,(T).5,0)),
                                      (T)1e10,(T)1);}
    }
}
//#####################################################################
// Function Add_Enclosed_Room
//#####################################################################
void Add_Enclosed_Room()
{
    if(!fp) fp = 7; // This one looks the best.
    if(fracture_walls) FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/fracture_pattern-%d",data_directory.c_str(),fp),fracture_pattern);

    for(int i=1;i<=4;++i){
        TV edge_lengths((T)1,(T).2,(T)1);
        TV_INT dimensions(31,3,31);BOX<TV> box((T)-.5*edge_lengths,(T).5*edge_lengths);TV dx=edge_lengths/TV(dimensions-1);
        TV_INT levelset_resolution(101,21,101);TV levelset_dx=edge_lengths/TV(levelset_resolution-1);int ghost_cells=2;
        RIGID_BODY<TV>* rigid_body=new RIGID_BODY<TV>(rigid_body_collection,true);

        rigid_body->Add_Structure(*TESSELLATION::Generate_Triangles(BOX<TV>((T)-.5*edge_lengths,(T).5*edge_lengths)));
        rigid_body->simplicial_object->mesh.Initialize_Adjacent_Elements();
        LEVELSET_IMPLICIT_OBJECT<TV>& implicit_object=*LEVELSET_IMPLICIT_OBJECT<TV>::Create();
        implicit_object.levelset.grid.Initialize(levelset_resolution+2*ghost_cells,box.Thickened(levelset_dx.x*ghost_cells),false);
        implicit_object.levelset.phi.Resize(implicit_object.levelset.grid.Domain_Indices());implicit_object.levelset.phi.Fill(FLT_MAX);
        implicit_object.Update_Box();implicit_object.Update_Minimum_Cell_Size();
        for(NODE_ITERATOR iterator(implicit_object.levelset.grid);iterator.Valid();iterator.Next())
            implicit_object.levelset.phi(iterator.index)=box.Signed_Distance(iterator.Location());
        rigid_body->Add_Structure(implicit_object);
        ROTATION<TV> rotation((T)i*(T)pi/2,TV(0,1,0));
        rigid_body->Rotation()=rotation*ROTATION<TV>((T)-pi/2,TV(0,0,(T)1))*ROTATION<TV>((T)pi/2,TV(0,(T)1,0));
        rigid_body_collection.Add_Rigid_Body_And_Geometry(rigid_body);
        if(fracture_walls){
            rigid_body->fracture_threshold=(T)0;
            Set_Rigid_Body_Parameters(rigid_body->particle_index,STRING_UTILITIES::string_sprintf("wall %d",i).c_str(),TV(1,0,1)+rotation.Rotate(TV((T).5,(T).5,(T).1)),
                                      (T)1e3,(T).5);}
        else{
            Set_Rigid_Body_Parameters(rigid_body->particle_index,STRING_UTILITIES::string_sprintf("wall %d",i).c_str(),TV(1,0,1)+rotation.Rotate(TV((T).5,(T).5,(T).1)),
                                      (T)1e10,(T)1);}}
    if(0){
        TV edge_lengths((T)1,(T).2,(T)1);
        TV_INT dimensions(31,3,31);BOX<TV> box((T)-.5*edge_lengths,(T).5*edge_lengths);TV dx=edge_lengths/TV(dimensions-1);
        TV_INT levelset_resolution(101,21,101);TV levelset_dx=edge_lengths/TV(levelset_resolution-1);int ghost_cells=2;
        RIGID_BODY<TV>* rigid_body=new RIGID_BODY<TV>(rigid_body_collection,true);

        rigid_body->Add_Structure(*TESSELLATION::Generate_Triangles(BOX<TV>((T)-.5*edge_lengths,(T).5*edge_lengths)));
        rigid_body->simplicial_object->mesh.Initialize_Adjacent_Elements();
        LEVELSET_IMPLICIT_OBJECT<TV>& implicit_object=*LEVELSET_IMPLICIT_OBJECT<TV>::Create();
        implicit_object.levelset.grid.Initialize(levelset_resolution+2*ghost_cells,box.Thickened(levelset_dx.x*ghost_cells),false);
        implicit_object.levelset.phi.Resize(implicit_object.levelset.grid.Domain_Indices());implicit_object.levelset.phi.Fill(FLT_MAX);
        implicit_object.Update_Box();implicit_object.Update_Minimum_Cell_Size();
        for(NODE_ITERATOR iterator(implicit_object.levelset.grid);iterator.Valid();iterator.Next())
            implicit_object.levelset.phi(iterator.index)=box.Signed_Distance(iterator.Location());
        rigid_body->Add_Structure(implicit_object);
        ROTATION<TV> rotation((T)pi/2,TV(1,0,0));
        rigid_body->Rotation()=ROTATION<TV>((T)pi/2,TV(0,(T)1,0));
        rigid_body_collection.Add_Rigid_Body_And_Geometry(rigid_body);
        if(fracture_walls){
            rigid_body->fracture_threshold=(T)0;
            Set_Rigid_Body_Parameters(rigid_body->particle_index,"ceiling",TV(1,0,1)+rotation.Rotate(TV((T)0,(T)0,(T)-1)),
                                      (T)1e3,(T).5);}
        else{
            Set_Rigid_Body_Parameters(rigid_body->particle_index,"ceiling",TV(1,0,1)+rotation.Rotate(TV((T)0,(T)0,(T)-1)),
                                      (T)1e10,(T)1);}
    }
}
//#####################################################################
// Function Post_Initialization
//#####################################################################
void Post_Initialization() PHYSBAM_OVERRIDE
{
    RIGID_BODY_COLLISIONS<TV>& collisions=*solids_evolution->rigid_body_collisions;
    collisions.Set_Push_Out_Level_Iterations(1);
    if(test_number==19) collisions.collision_manager=collision_manager;

    if(test_number==13 || test_number==19){ // HACK -- adds fluid_collision_bodies after they've been read in from file.
        LOG::cout<<"Adding bodies to the simulation"<<std::endl;
        for(int i=1;i<=deformable_objects_to_simulate.m;++i){
            LOG::cout<<"\tadded DEFORMABLE body "<<i<<std::endl;
            DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& collision_structure=*deformable_objects_to_simulate(i);
            collision_structure.object.Initialize_Hierarchy();
            Add_To_Fluid_Simulation(collision_structure);}
        for(int i=1;i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++)
            if(rigid_body_collection.Is_Active(i) && rigid_body_collection.Rigid_Body(i).simplicial_object)
                rigid_body_collection.Rigid_Body(i).simplicial_object->Update_Triangle_List();
        fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection.rigid_geometry_collection);}
}
//#####################################################################
// Function Bunny
//#####################################################################
void Bunny()
{
    int cannon_id,ammo_id;
    if(!fp) fp = 6; // This one looks the best.
    if(fracture_walls) FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/fracture_pattern-%d",data_directory.c_str(),fp),fracture_pattern);

    { // The cannon.
        RIGID_BODY<TV>* rigid_body=new RIGID_BODY<TV>(rigid_body_collection,true);
        TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
        FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/Rigid_Bodies/cannon.tri.gz",*surface);
        surface->mesh.Initialize_Adjacent_Elements();
        LEVELSET_IMPLICIT_OBJECT<TV>* levelset=LEVELSET_IMPLICIT_OBJECT<TV>::Create();
        FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/Rigid_Bodies/cannon.phi.gz",*levelset);
        rigid_body->Add_Structure(*surface);
        rigid_body->Add_Structure(*levelset);
        rigid_body->Rotation()=ROTATION<TV>((T)pi,TV(0,(T)1,0))*ROTATION<TV>((T)-pi/2,TV(0,0,(T)1))*ROTATION<TV>((T)pi/5,TV(0,(T)1,0));
        rigid_body_collection.Add_Rigid_Body_And_Geometry(rigid_body);
        Set_Rigid_Body_Parameters(rigid_body->particle_index,"Tube",TV(38,0,0),(T)1e10,(T)1);
        cannon_id=rigid_body->particle_index;
    }
    { // The pedistal.
        RIGID_BODY<TV>& ground=solid_tests.Add_Analytic_Box(TV((T)3,(T)7.5,(T)3));
        ground.simplicial_object->mesh.Initialize_Adjacent_Elements();
        Set_Rigid_Body_Parameters(ground.particle_index,"ground",TV(38,(T)-6.25,0),(T)1e10,1);
    }
    { // The Ammo.
        RIGID_BODY<TV>& rigid_body=solid_tests.Add_Rigid_Body("sphere",(T)1.45,(T).5,true,false);
        ammo_id=rigid_body.particle_index;
        rigid_body.simplicial_object->mesh.Initialize_Adjacent_Elements();
        rigid_body.X()=TV((T)35,0,0);
        rigid_body.Set_Mass((T)1e3);
        collision_manager=new RIGID_BODY_COLLISION_MANAGER_HASH;
        collision_manager->hash.Insert(PAIR<int,int>(cannon_id,ammo_id));
        collision_manager->hash.Insert(PAIR<int,int>(ammo_id,cannon_id));
    }
    { // The bunny.
        if(fracture_walls){fp=6;FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/fracture_pattern-%d",data_directory.c_str(),fp),fracture_pattern);}
        RIGID_BODY<TV>& rigid_body=solid_tests.Add_Rigid_Body("bunny",6,(T).5,true,false);(void)rigid_body;
        rigid_body.simplicial_object->mesh.Initialize_Adjacent_Elements();
        rigid_body.X()=TV(0,(T)-1,0);
        rigid_body.Rotation()=ROTATION<TV>((T)pi/2,TV(0,1,0))*rigid_body.Rotation();
        rigid_body.fracture_threshold=(T)1e0;
    }
    { // The pedistal.
        RIGID_BODY<TV>& ground=solid_tests.Add_Analytic_Box(TV((T)3,(T)7,(T)3));
        ground.simplicial_object->mesh.Initialize_Adjacent_Elements();
        Set_Rigid_Body_Parameters(ground.particle_index,"ground",TV(0,(T)-6.5,0),(T)1e10,1);
    }
}
//#####################################################################
// Function Add_Wall_High_Resolution
//#####################################################################
void Add_Wall_High_Resolution(const T x_delta,const T unit_mass)
{
    T_GRID& grid=(fluids_parameters.mpi_grid)?fluids_parameters.mpi_grid->global_grid:*fluids_parameters.grid;
    const char* wallfile="short_plank_subdivided";
    T scale=((T)1/(T)1.5)*grid.Domain().max_corner.z;
    T mass=unit_mass*scale*scale*scale;
    T x_center=grid.Domain().max_corner.x*.5+x_delta,y_top=grid.Domain().min_corner.y,z_center=0;
    T x_pos=x_center,y_pos=y_top+scale,z_pos=z_center;
    T mu=(T).5;

    RIGID_BODY<TV>& wall=solid_tests.Add_Rigid_Body(wallfile,scale,mu);
    Set_Rigid_Body_Parameters(wall.particle_index,"wall",TV(x_pos,y_pos,z_pos),mass,mu);
    wall.Rotation()=ROTATION<TV>((T)pi/2,TV(0,0,1));
}
//#####################################################################
// Function Add_Stack_Of_Squares
//#####################################################################
void Add_Stack_Of_Squares()
{
    int parameter=0;
    const char* boxfile=parameter?"box":"subdivided_box";
    T stack_mu=(T).5;

    T box_scale_x=(T)1,box_scale_y=(T)1,box_scale_z=(T)1;

    T x_center=0,y_top=-1,z_center=0;
    T scale=.1;
    T epsilony=0;

    T x_scale=box_scale_x*scale,y_scale=box_scale_y*scale,z_scale=box_scale_z*scale;
    T x_pos=x_center,y_pos=y_top+y_scale,z_pos=z_center;
    RIGID_BODY<TV>& box1=solid_tests.Add_Rigid_Body(boxfile,scale,stack_mu);
    Set_Rigid_Body_Parameters(box1.particle_index,"box1",TV(x_pos,y_pos,z_pos),solid_mass*20,stack_mu);

    y_top=y_pos+y_scale+epsilony;
    
    x_scale=box_scale_x*scale;y_scale=box_scale_y*scale;z_scale=box_scale_z*scale;
    x_pos=x_center,y_pos=y_top+y_scale,z_pos=z_center;
    RIGID_BODY<TV>& box4=solid_tests.Add_Rigid_Body(boxfile,scale,stack_mu);
    Set_Rigid_Body_Parameters(box4.particle_index,"box4",TV(x_pos,y_pos,z_pos),solid_mass*10,stack_mu);
    
    y_top=y_pos+y_scale+epsilony;

    x_scale=box_scale_x*scale;y_scale=box_scale_y*scale;z_scale=box_scale_z*scale;
    x_pos=x_center,y_pos=y_top+y_scale,z_pos=z_center;
    RIGID_BODY<TV>& box5=solid_tests.Add_Rigid_Body(boxfile,scale,stack_mu);
    Set_Rigid_Body_Parameters(box5.particle_index,"box5",TV(x_pos,y_pos,z_pos),solid_mass*10,stack_mu);
    y_top=y_pos+y_scale+epsilony;
}
//#####################################################################
// Function Add_Stack_Of_Bodies
//#####################################################################
void Add_Stack_Of_Bodies()
{
    T_GRID& grid=(fluids_parameters.mpi_grid)?fluids_parameters.mpi_grid->global_grid:*fluids_parameters.grid;
    int parameter=0;
    const char* boxfile=parameter?"box":"subdivided_box";
    const char* plankfile=parameter?"unsubdivided_plank":"plank_refined";
    T stack_mu=(T)1;
    T mass=(T).08*grid.Domain().Size();

    T box_scale_x=(T)1,box_scale_y=(T)1,box_scale_z=(T)1;
    T plank_scale_x=(T)1,plank_scale_y=(T).25,plank_scale_z=(T)5;

    T x_center=grid.Domain().max_corner.x*.133,y_top=grid.Domain().min_corner.y,z_center=0;
    T scale=grid.Domain().max_corner.x*.067;
    T epsilonx=.05,epsilony=0,epsilonz=.05;

    T x_scale=box_scale_x*scale,y_scale=box_scale_y*scale,z_scale=box_scale_z*scale;;
    T x_pos=x_center-x_scale-epsilonx,y_pos=y_top+y_scale,z_pos=z_center;
    RIGID_BODY<TV>& box1=solid_tests.Add_Rigid_Body(boxfile,scale,stack_mu);
    Set_Rigid_Body_Parameters(box1.particle_index,"box1",TV(x_pos,y_pos,z_pos),mass*20,stack_mu);

    x_scale=box_scale_x*scale;y_scale=box_scale_y*scale;z_scale=box_scale_z*scale;
    x_pos=x_center+x_scale+epsilonx,y_pos=y_top+y_scale,z_pos=z_center;
    RIGID_BODY<TV>& box2=solid_tests.Add_Rigid_Body(boxfile,scale,stack_mu);
    Set_Rigid_Body_Parameters(box2.particle_index,"box2",TV(x_pos,y_pos,z_pos),mass*20,stack_mu);
    y_top=y_pos+y_scale+epsilony;
    
    x_scale=plank_scale_x*scale;y_scale=plank_scale_y*scale;z_scale=plank_scale_z*scale;
    x_pos=x_center,y_pos=y_top+y_scale,z_pos=z_center;
    RIGID_BODY<TV>& box3=solid_tests.Add_Rigid_Body(plankfile,scale,stack_mu);
    Set_Rigid_Body_Parameters(box3.particle_index,"plank1",TV(x_pos,y_pos,z_pos),mass*20,stack_mu);
    y_top=y_pos+y_scale+epsilony;

    x_scale=box_scale_x*scale;y_scale=box_scale_y*scale;z_scale=box_scale_z*scale;
    x_pos=x_center,y_pos=y_top+y_scale,z_pos=z_center-z_scale;
    RIGID_BODY<TV>& box4=solid_tests.Add_Rigid_Body(boxfile,scale,stack_mu);
    Set_Rigid_Body_Parameters(box4.particle_index,"box3",TV(x_pos,y_pos,z_pos),mass*10,stack_mu);

    x_scale=box_scale_x*scale;y_scale=box_scale_y*scale;z_scale=box_scale_z*scale;
    x_pos=x_center,y_pos=y_top+y_scale,z_pos=z_center+z_scale+epsilonz;
    RIGID_BODY<TV>& box5=solid_tests.Add_Rigid_Body(boxfile,scale,stack_mu);
    Set_Rigid_Body_Parameters(box5.particle_index,"box4",TV(x_pos,y_pos,z_pos),mass*10,stack_mu);
    y_top=y_pos+y_scale+epsilony;

    x_scale=plank_scale_x*scale;y_scale=plank_scale_y*scale;z_scale=plank_scale_z*scale;
    x_pos=x_center,y_pos=y_top+y_scale,z_pos=z_center;
    RIGID_BODY<TV>& box6=solid_tests.Add_Rigid_Body(plankfile,scale,stack_mu);
    Set_Rigid_Body_Parameters(box6.particle_index,"plank2",TV(x_pos,y_pos,z_pos),mass*20,stack_mu);
    y_top=y_pos+y_scale+epsilony;

    x_pos=x_center,y_pos=y_top+scale,z_pos=z_center;
    int sphere=rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/sphere_refined",scale,true,true,false);
    Set_Rigid_Body_Parameters(sphere,"sphere",TV(x_pos,y_pos,z_pos),mass*100,stack_mu);
}
//#####################################################################
// Function Add_Deformable_Ball
//#####################################################################
void Add_Deformable_Ball()
{
    T_GRID& grid=(fluids_parameters.mpi_grid)?fluids_parameters.mpi_grid->global_grid:*fluids_parameters.grid;
    T solid_density=10;
    T scale=.75;
    T y_top=grid.Domain().min_corner.y;
    T x_center=shock_radius*(T).8,y_center=y_top+(T)scale;

    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=solid_tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(x_center,y_center,0))),true,true,solid_density,scale);
    solid_tests.Set_Mass_Of_Particles(tetrahedralized_volume,solid_density,false);
}
//#####################################################################
// Function Finalize_Deformable_Bodies
//#####################################################################
void Finalize_Deformable_Bodies()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    //solid_body_collection.collision_body_list.Add_Bodies(*rigid_body_collection.rigid_geometry_collection.collision_body_list);

    // correct number nodes
    for(int i=1;i<=solid_body_collection.deformable_body_collection.deformable_geometry.structures.m;i++) solid_body_collection.deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    // Add forces
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
    solid_body_collection.Add_Force(Create_Edge_Springs(tetrahedralized_volume,(T)1e4,(T)5));
    solid_body_collection.Add_Force(Create_Altitude_Springs(tetrahedralized_volume,(T)1e4));

    // Add to fluid simulation
    TRIANGULATED_SURFACE<T>& triangulated_surface=tetrahedralized_volume.Get_Boundary_Object();
    triangulated_surface.mesh.Initialize_Incident_Elements();
    triangulated_surface.mesh.Initialize_Adjacent_Elements();

    solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&triangulated_surface);
    DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>* deformable_collisions=new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(triangulated_surface);
    if(test_number == 13 || test_number == 19)
        deformable_objects_to_simulate.Append(deformable_collisions);
    else
        Add_To_Fluid_Simulation(*deformable_collisions);
}
//#####################################################################
// Function Intialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{   
    if(no_solids || test_number==1 || test_number==20 || test_number==21) return;

    switch (test_number){
        case 2: Add_Sphere();
                break;
        case 3: Add_Wall();
                break;
        case 4: Add_Sphere();
                break;
        case 5: Add_Sphere();
                break;
        case 6: Add_Sphere();
                break;
        case 7: Add_Stack_Of_Bodies();
                break;
        case 8: Add_Stack_Of_Squares();
                break;
        case 9: Add_Wall_High_Resolution(0,(T)1);
                break;
        case 10:Add_Wall_High_Resolution(0,(T)850);
                break;
        case 11:Add_Destructive_Wall();
                break;
        case 12:Add_Destructive_Wall();
                break;
        case 13:Add_Room();
                break;
        case 14:Add_Wall_High_Resolution((T)-1,(T)10);
                break;
        case 15:Add_Wall_High_Resolution((T)-1,(T)500);
                break;
        case 16:Add_Enclosed_Room();
                break;
        case 17:Add_Deformable_Ball();
                break;
        case 18: break;
        case 19: Bunny();
                 break;
        default: PHYSBAM_FATAL_ERROR("wrong test number");}

    if(test_number==17) Finalize_Deformable_Bodies();

    if(test_number!=3 && test_number!=11 && test_number!=12 && test_number!=13 && test_number!=14 && test_number!=15 && !test_number==16 && test_number!=19)
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Rigid_Body_Walls(*this);
    else Add_Ground(true);

    if(use_solids_gravity) solid_body_collection.Add_Force(new GRAVITY<TV>(solid_body_collection.deformable_body_collection.particles,rigid_body_collection,true,true,solid_gravity));

    // Half forces for SPD
    for(int i=1;i<=solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->compute_half_forces=true;
    for(int k=1;k<=solid_body_collection.deformable_body_collection.deformables_forces.m;k++) solid_body_collection.deformable_body_collection.deformables_forces(k)->compute_half_forces=true;
    for(int i=1;i<=solid_body_collection.rigid_body_collection.rigids_forces.m;i++) solid_body_collection.rigid_body_collection.rigids_forces(i)->compute_half_forces=true;
}
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    if(fluids_parameters.use_slip){
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).run_self_tests=run_self_tests;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).output_iterators=output_iterators;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_matrix=print_matrix;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_rhs=print_rhs;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_each_matrix=print_each_matrix;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_poisson_matrix=print_poisson_matrix;
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution).print_index_map=print_index_map;}
    else if(fluids_parameters.fluid_affects_solid){
        dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION<TV>&>(*solids_evolution).print_matrix_rhs_and_solution=print_matrix;}
}
void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
{
    if(transition_to_incompressible){
        eos_smooth_transition->Set_Current_Time(time);
        //if(time>eos_smooth_transition->t_start_transition) fluids_parameters.euler->euler_projection.Set_Transition_To_Using_Implicit_Pressure(true);
    }
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
//    if(fracture_walls) FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/fracture_pattern-%d",data_directory.c_str(),fp),fracture_pattern);
    if(fluids_parameters.use_slip && solids_fluids_parameters.use_fluid_rigid_fracture){
        FRACTURE_PATTERN<T>* fp=(dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution)).fracture_pattern;
        if(solids_fluids_parameters.use_fluid_rigid_fracture && fp && !fp->regions.m){
            for(int i=1;i<=fracture_pattern.regions.m;i++) fp->regions.Append(fracture_pattern.regions(i));}}
}
//#####################################################################
// Function Shrink_Levelset
//#####################################################################
void Shrink_Levelset(GRID<TV>& grid,ARRAY<T,VECTOR<int,3> >& phi,int boundary,TV_INT& center)
{
    TV_INT min_adjust;
    RANGE<TV_INT> inside=RANGE<TV_INT>::Zero_Box().Thickened(-INT_MAX),domain=grid.Domain_Indices();
    for(NODE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) if(phi(iterator.index)<0) inside.Enlarge_To_Include_Point(iterator.index);
    inside=RANGE<TV_INT>::Intersect(inside.Thickened(boundary),domain);
    min_adjust=inside.min_corner-domain.min_corner;
    center-=min_adjust;
    GRID<TV> grid_new(inside.Edge_Lengths()+1,BOX<TV>(grid.Node(inside.min_corner),grid.Node(inside.max_corner)),false);
    ARRAY<T,VECTOR<int,3> > phi_new(grid_new.Domain_Indices());
    for(NODE_ITERATOR iterator(grid_new);iterator.Valid();iterator.Next()) phi_new(iterator.index)=phi(iterator.index+min_adjust);
    grid=grid_new;
    phi=phi_new;
}
//#####################################################################
// Function Create_Wall_Pattern
//#####################################################################
void Create_Wall_Pattern()
{
    FRACTURE_PATTERN<T> fp;
    TV_INT resolution(301,21,201);TV edge_lengths((T)1,(T)3,(T)4);TV dx=edge_lengths/TV(resolution-1);int ghost_cells=2;
    TV original_half_edge_length(edge_lengths/2);
    TV_INT actual_resolution=resolution+ghost_cells*2;TV actual_edge_lengths=TV(actual_resolution-1)*dx;
    TV half_edge_length(actual_edge_lengths/2);

    for(int i=2;i<=26;i++){
        GRID<TV>& local_grid=*new GRID<TV>(actual_resolution,BOX<TV>(-half_edge_length,half_edge_length),false);
        ARRAY<T,VECTOR<int,3> >& local_phi=*new ARRAY<T,VECTOR<int,3> >(local_grid.Domain_Indices());local_phi.Fill(FLT_MAX);
        LEVELSET_IMPLICIT_OBJECT<TV>* refined_lio=LEVELSET_IMPLICIT_OBJECT<TV>::Create();
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/wall/fragment.%d.phi",data_directory.c_str(),i),refined_lio->levelset);
        LEVELSET_IMPLICIT_OBJECT<TV>* lio=new LEVELSET_IMPLICIT_OBJECT<TV>(local_grid,local_phi);
        for(NODE_ITERATOR iterator(local_grid);iterator.Valid();iterator.Next())
            local_phi(iterator.index)=refined_lio->levelset.Extended_Phi(iterator.Location());
        TV_INT center_index=local_grid.Domain_Indices().Center();
        Shrink_Levelset(local_grid,local_phi,2,center_index);
        TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/wall/fragment.%d.tri",data_directory.c_str(),i),*surface);
        surface->Update_Number_Nodes();
        FRACTURE_REGION<T>* fr=new FRACTURE_REGION<T>(surface,lio,false);
        fr->fracture_offset=center_index;
        fp.regions.Append(fr);}
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/Fracture_Patterns/fracture_pattern-%d",data_directory.c_str(),test_number),fp);
    return;
}
//#####################################################################
};
}
#endif

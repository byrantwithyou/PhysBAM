//#####################################################################
// Copyright 2004-2005, Eran Guendelman, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_FLUID_COUPLING_TEST 
//##################################################################### 
#ifndef __SOLID_FLUID_COUPLING_TEST__
#define __SOLID_FLUID_COUPLING_TEST__

#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include "../WATER_MELTING_EXAMPLE_3D.h"
#include <Level_Sets/EXTRAPOLATION_3D.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>
namespace PhysBAM{

template<class T,class RW=T>
class MELTING_TWO_WAY:public WATER_MELTING_EXAMPLE_3D<T,RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::solids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::verbose_dt;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_time;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_frame_title;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::data_directory;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::abort_when_dt_below;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Use_Thin_Shells_Fluid_Coupling_Defaults;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::set_fluid_boundary_conditions_for_pressure_jump_solve;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::need_to_compute_pressure_coupling;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::collision_body_affected_by_fluid;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::use_object_pseudo_velocity_for_boundary_conditions;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::rigid_body_fluid_forces;
    using MELTING_EXAMPLE_3D<T,RW,SOLIDS_FLUIDS_EXAMPLE_3D<RW> >::melting_parameters;using WATER_MELTING_EXAMPLE_3D<T,RW>::Construct_Levelsets_For_Objects;using WATER_MELTING_EXAMPLE_3D<T,RW>::maximum_velocity;
    using WATER_MELTING_EXAMPLE_3D<T,RW>::phi_objects;
    using MELTING_EXAMPLE_3D<T,RW,SOLIDS_FLUIDS_EXAMPLE_3D<RW> >::initialized;using MELTING_EXAMPLE_3D<T,RW,SOLIDS_FLUIDS_EXAMPLE_3D<RW> >::write_frame_title;
    using WATER_MELTING_EXAMPLE_3D<T,RW>::V_objects;using WATER_MELTING_EXAMPLE_3D<T,RW>::phi_objects;using WATER_MELTING_EXAMPLE_3D<T,RW>::use_melting;

    int example_number;

    T initial_water_level;

    PARAMETER_LIST parameter_list;
    T domain_scale;

    int number_of_objects;
    ARRAY<VECTOR<T,3> > initial_position;
    ARRAY<QUATERNION<T> > initial_orientation;
    ARRAY<VECTOR<T,3> > initial_velocity;
    ARRAY<VECTOR<T,3> > initial_angular_velocity;

    ARRAY<int,VECTOR<int,3> > id_objects;
    //ARRAY<T,VECTOR<int,3> > phi_objects;
    //ARRAY<VECTOR<T,3> ,VECTOR<int,3> > V_objects;
    bool add_gravity_to_object;
    bool use_next_dt_for_velocity_adjustment;
    bool subtract_last_pressure_force_from_object;
    bool extrapolate_combined_solid_fluid_velocity;
    bool set_faces_to_invalid;
    T velocity_extrapolation_band_width;
    bool use_extrapolated_velocity_for_laplace;
    bool only_use_object_velocity_in_mixed_face;
    bool average_velocities_in_mixed_face;
    bool mixed_face_is_solid;
    bool use_effective_velocity_for_poisson;

    ARRAY<RIGID_BODY<TV>*> bodies_in_fluid;
    ARRAY<T> body_densities;

    T melting_speed;
    T start_melting_time;
    T mu;
    T epsilon;

    T scaling;
    int m,n,mn;

    LEVELSET_3D<GRID<TV> >* icecube_levelset;
    VECTOR<T,3> levelset_center;

    MELTING_TWO_WAY(const int example_number_input)
        :example_number(example_number_input),domain_scale(1),m(11),n(11),mn(11)
    {
        // use_melting will be turned on later
        use_melting=true;

        data_directory="/n/bcc/data/losasso/PhysBAM-Full/Public_Data";


        printf("Running *** MELTING_TWO_WAY *** number %d\n",example_number);
        {std::string filename=STRING_UTILITIES::string_sprintf("Melting_Two_Way/example_%d.param",example_number);
        if(FILE_UTILITIES::File_Exists(filename)){std::cout << "Reading parameter file '" << filename << "'" << std::endl;parameter_list.Read(filename);}}

        // TODO: is lockstep necessary?
        add_gravity_to_object=false;                      Set_Parameter(add_gravity_to_object,"add_gravity_to_object");
        use_next_dt_for_velocity_adjustment=false;        Set_Parameter(use_next_dt_for_velocity_adjustment,"use_next_dt_for_velocity_adjustment");
        subtract_last_pressure_force_from_object=true;    Set_Parameter(subtract_last_pressure_force_from_object,"subtract_last_pressure_force_from_object");
        extrapolate_combined_solid_fluid_velocity=true;   Set_Parameter(extrapolate_combined_solid_fluid_velocity,"extrapolate_combined_solid_fluid_velocity");
        set_faces_to_invalid=false;                       Set_Parameter(set_faces_to_invalid,"set_faces_to_invalid");
        velocity_extrapolation_band_width=0;
        use_extrapolated_velocity_for_laplace=false;      Set_Parameter(use_extrapolated_velocity_for_laplace,"use_extrapolated_velocity_for_laplace");
        only_use_object_velocity_in_mixed_face=false;     Set_Parameter(only_use_object_velocity_in_mixed_face,"only_use_object_velocity_in_mixed_face");
        average_velocities_in_mixed_face=false;           Set_Parameter(average_velocities_in_mixed_face,"average_velocities_in_mixed_face");
        mixed_face_is_solid=false;                        Set_Parameter(mixed_face_is_solid,"mixed_face_is_solid");
        use_effective_velocity_for_poisson=false;         Set_Parameter(use_effective_velocity_for_poisson,"use_effective_velocity_for_poisson");

        Set_Parameter(domain_scale,"domain_scale");

        mu=0;epsilon=.2;
        Set_Parameter(mu,"mu");
        Set_Parameter(epsilon,"epsilon");


        scaling=.17;

        T object_height_offset=0;T object_horizontal_offset=0;
        if(example_number==1)object_height_offset=0;
        if(example_number==2)object_height_offset=0;
        if(example_number==3)object_height_offset=-.4;
        if(example_number==4)object_height_offset=-.4;
        if(example_number==5)object_height_offset=-.4;
        if(example_number==6)object_height_offset=-.4;
        if(example_number==7)object_height_offset=-.4;

        // Common parameters
        output_directory="Melting_Two_Way/output";
        first_frame=0;last_frame=1000;frame_rate=60;
        restart=false;restart_frame=0;

        // Fluids parameters
        fluids_parameters.use_solid_fluid_coupling=false;
        fluids_parameters.use_volumetric_solid_fluid_coupling=true; Set_Parameter(fluids_parameters.use_volumetric_solid_fluid_coupling,"fluids_parameters.use_volumetric_solid_fluid_coupling");
        if(example_number==1)fluids_parameters.grid.Initialize(31,51,31,0,   domain_scale*0.6, 0,domain_scale, 0,domain_scale*0.6);
        if(example_number==2)fluids_parameters.grid.Initialize(85,141,85,0,  domain_scale*0.6, 0,domain_scale, 0,domain_scale*0.6);
        if(example_number==3)fluids_parameters.grid.Initialize(101,101,101,0,domain_scale*0.6, 0,domain_scale*0.6, 0,domain_scale*0.6);
        if(example_number==4)fluids_parameters.grid.Initialize(121,121,121,0,domain_scale*0.6, 0,domain_scale*0.6, 0,domain_scale*0.6);
        if(example_number==5)fluids_parameters.grid.Initialize(101,101,101,0,domain_scale*0.6, 0,domain_scale*0.6, 0,domain_scale*0.6);
        if(example_number==6)fluids_parameters.grid.Initialize(101,101,101,0,domain_scale*0.6, 0,domain_scale*0.6, 0,domain_scale*0.6);
        if(example_number==7)fluids_parameters.grid.Initialize(141,141,141,0,domain_scale*0.6, 0,domain_scale*0.6, 0,domain_scale*0.6);
//        if(example_number==5)fluids_parameters.grid.Initialize(151,151,151,0,domain_scale*0.6, 0,domain_scale*0.6, 0,domain_scale*0.6);
//        if(example_number==6)fluids_parameters.grid.Initialize(151,151,151,0,domain_scale*0.8, 0,domain_scale*0.8, 0,domain_scale*0.8);
        fluids_parameters.number_particles_per_cell=8;
        fluids_parameters.reseeding_frame_rate=10;
        fluids_parameters.incompressible_iterations=10000;
        fluids_parameters.cfl=0.9;
        fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[1][2]=true;fluids_parameters.domain_walls[2][1]=true;fluids_parameters.domain_walls[2][2]=false;
        fluids_parameters.use_removed_positive_particles=false;fluids_parameters.use_removed_negative_particles=true;
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;fluids_parameters.variable_viscosity=false;
        fluids_parameters.second_order_pressure=true;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;
        fluids_parameters.write_particles=true;fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=true;
        fluids_parameters.store_particle_ids=true;
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.enforce_divergence_free_extrapolation=false;
        fluids_parameters.levelset_substeps=(T).95;
        fluids_parameters.restart_data_write_rate=1;
        fluids_parameters.simulate=true;

        // Solids parameters
        solids_parameters.gravity=9.8;
        solids_parameters.cfl=(T).9;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies=true;
        solids_parameters.deformable_body_parameters.write=true;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        verbose_dt=true;

        solids_parameters.rigid_body_evolution_parameters.minimum_rigid_body_time_step_fraction=.001;

        // Melting parameters
        melting_parameters.maximum_depth=2;
//        melting_parameters.maximum_depth=1;
        melting_parameters.refine_near_interface=true;
        melting_parameters.refine_for_high_deformation=false;
        melting_parameters.rigid_body_coupling_density=800;

        initial_water_level=domain_scale*0.6123+object_height_offset;         Set_Parameter(initial_water_level,"initial_water_level");

        // Debugging
        fluids_parameters.write_debug_data=true;
        abort_when_dt_below=1e-7;
        write_time=true;
        write_frame_title=true;

        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Parameters_From_Parameter_List(*this,parameter_list);

        icecube_levelset=new LEVELSET_3D<GRID<TV> >(*(new GRID<TV>()),*(new ARRAY<T,VECTOR<int,3> >()));

        FILE_UTILITIES::Read_From_File<RW>("/n/bcc/data/losasso/PhysBAM-Full/Public_Data/Rigid_Bodies/icecube_smooth.phi",*icecube_levelset);

        VECTOR<T,3> positions[]={VECTOR<T,3>(.15333,.78123,0.2),VECTOR<T,3>(.20333,.80123,0.45),VECTOR<T,3>(.35333,.73123,0.2),
                                  VECTOR<T,3>(.45333,.93123,0.2),VECTOR<T,3>(.40333,.78123,0.4),VECTOR<T,3>(.25333,.98123,0.37)};

        number_of_objects=6;

        initial_position.Resize(number_of_objects);
        initial_orientation.Resize(number_of_objects);
        initial_velocity.Resize(number_of_objects);
        initial_angular_velocity.Resize(number_of_objects);

        RANDOM_NUMBERS random;random.Set_Seed(1234);// set the seed so that the example is reproducible
        for(int i=1;i<=number_of_objects;i++){
            initial_position(i)=positions[i-1]+VECTOR<T,3>(object_horizontal_offset,object_height_offset,object_horizontal_offset);
            do{initial_orientation(i).s=random.Get_Uniform_Number(0,1);initial_orientation(i).v=random.Get_Uniform_Vector(VECTOR<T,3>(-1,-1,-1),VECTOR<T,3>(1,1,1));}while(initial_orientation(i).Magnitude()>1);
            initial_orientation(i).Normalize();
            initial_angular_velocity(i)=random.Get_Uniform_Vector(VECTOR<T,3>(-2,-2,-2),VECTOR<T,3>(2,2,2));}


        melting_speed=.03;
        if(example_number==5)
            melting_speed=.015;
        if(example_number==2||example_number==6)
            melting_speed=.005;
        if(example_number==3||example_number==4||example_number==7)
            melting_speed=.003;
    }

    ~MELTING_TWO_WAY() 
    {}


//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi(const int object,ARRAY<T>& phi)
{
    RED_GREEN_GRID_3D<T>& grid=melting_parameters.levelsets(object)->grid;
    ARRAY<VECTOR<T,3> >& node_locations=grid.Node_Locations();
    for(int p=1;p<=phi.m;p++){
        phi(p)=icecube_levelset->Extended_Phi(((T)1/scaling)*(node_locations(p)-levelset_center)+levelset_center);}
}
//#####################################################################
// Function Initialize_Levelset_Velocity
//#####################################################################
void Initialize_Levelset_Velocity(const int object,ARRAY<VECTOR<T,3> >& V)
{
    for(int p=1;p<=V.m;p++) V(p)=VECTOR<T,3>(0,0,0);    
}
//#####################################################################
// Function Initialize_Levelset_Velocity
//#####################################################################
void Melting_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
{
    for(int object=1;object<=melting_parameters.levelsets.m;object++){
        LEVELSET_TETRAHEDRALIZED_VOLUME<T>& levelset=*melting_parameters.levelsets(object);
        RED_GREEN_GRID_3D<T>& grid=levelset.grid;
        int index=melting_parameters.body_index(object);
        if(!index)continue;
        RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list(index);
        FRAME<T> frame=rigid_body.Frame()*melting_parameters.rigid_body_grid_frames(object).Inverse();
        // map data to the cell based indices
        ARRAY<VECTOR<T,3> >& node_locations=grid.Node_Locations();
        
        LINEAR_INTERPOLATION<T,T> interpolation;
        for(int i=1;i<=node_locations.m;i++){
            T phi=interpolation.Clamped_To_Array(fluids_parameters.grid,fluids_parameters.particle_levelset_evolution.phi,frame*node_locations(i));

            T V=melting_speed*LEVELSET<T,VECTOR<T,3> >::Heaviside(-phi,2*fluids_parameters.grid.min_dx_dy_dz);
            melting_parameters.levelsets(object)->phi(i)+=V*dt;
        }
        }
    WATER_MELTING_EXAMPLE_3D<T,RW>::Melting_Substep(dt,time);
}
//#####################################################################
// Function Initialize_Particle_Positions_And_Velocities
//#####################################################################
void Initialize_Particle_Positions_And_Velocities(const int object)
{
    RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list(melting_parameters.body_index(object));
    //rigid_body.is_static=true;
    rigid_body.position=initial_position(object);
    rigid_body.velocity=initial_velocity(object);
    rigid_body.orientation=initial_orientation(object);
    rigid_body.angular_velocity=initial_angular_velocity(object);
    rigid_body.Update_Angular_Momentum();
}
//#####################################################################
// Function Initialize_Particle_Positions_And_Velocities
//#####################################################################
void Initialize_Grid(const int object,RED_GREEN_GRID_3D<T>& grid)
{
    printf("-------- %d ---------\n",object);

    GRID<TV>& levelset_grid=icecube_levelset->grid;
    levelset_center=levelset_grid.Domain().Center();
    T radius=scaling*.55*levelset_grid.Domain().Size().Max();
    grid.Initialize(GRID_3D<T>(4*m+1,4*n+1,4*mn+1,BOX_3D<T>(levelset_center).Thickened(radius)),melting_parameters.maximum_depth);
}




//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time)
{
    WATER_MELTING_EXAMPLE_3D<T,RW>::Update_Fluid_Parameters(dt,time);
/*    if(time>.8)
    {
        use_melting=true;melting_parameters.use_constant_melting_speed=true;
        melting_parameters.constant_melting_speed=min((time-2)*(T).1,melting_speed);}*/
}
//#####################################################################
// Function Set_Parameter
//#####################################################################
void Update_Melting_Substep_Parameters(const T dt,const T time)
{
    WATER_MELTING_EXAMPLE_3D<T,RW>::Update_Melting_Substep_Parameters(dt,time);
}

//#####################################################################
// Function Set_Parameter
//#####################################################################
template<class T2> void Set_Parameter(T2& parameter,const std::string& name)
{
    if(parameter_list.Is_Defined(name)){
        parameter=parameter_list.template Get_Parameter<T2>(name);
        std::cout << "[param] set " << name << " = " << parameter << std::endl;}
}
//#####################################################################
// Function Initialize_Forces
//#####################################################################
void Initialize_Forces()
{
    collision_body_affected_by_fluid.Clean_Memory();
    rigid_body_fluid_forces.Clean_Memory();
    fluids_parameters.collision_bodies_affecting_fluid.collision_bodies.Clean_Memory();
    bodies_in_fluid.Clean_Memory();
    body_densities.Clean_Memory();

    printf("\n---------- number of objects in collision lists %d %d\n\n",collision_body_affected_by_fluid.m,fluids_parameters.collision_bodies_affecting_fluid.collision_bodies.m);

    for(int object=1;object<=melting_parameters.body_index.m;object++){
        int index=melting_parameters.body_index(object);if(!index)continue;
        RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list(index);
        std::cout<<rigid_body.position<<std::endl;
        rigid_body.Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

        Add_Volumetric_Body_To_Fluid_Simulation(rigid_body,true,true,1);
        bodies_in_fluid.Append(&rigid_body);body_densities.Append(rigid_body.Volumetric_Density());
        rigid_body.Set_Coefficient_Of_Friction(mu);
        rigid_body.Set_Coefficient_Of_Restitution(epsilon);
    }
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Deformable_And_Rigid_Bodies()
{
    printf("XXX Initialize_Deformable_And_Rigid_Bodies()\n");
    for(int i=1;i<=number_of_objects;i++){
        Add_Melting_Object(melting_parameters.RIGID,0);}

    int id;

    VECTOR<T,3> center=fluids_parameters.grid.Domain().Center();

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies/ground",(T)domain_scale*.005);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction(mu);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Restitution(epsilon);
    solids_parameters.rigid_body_parameters.list(id)->position=VECTOR<T,3>(center.x,fluids_parameters.grid.ymin,center.z);
    solids_parameters.rigid_body_parameters.list(id)->is_static=true;

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies/ground",(T)domain_scale*.005);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction(mu);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Restitution(epsilon);
    solids_parameters.rigid_body_parameters.list(id)->position=VECTOR<T,3>(fluids_parameters.grid.xmax,center.y,center.z);
    solids_parameters.rigid_body_parameters.list(id)->orientation=QUATERNION<T>(.5*pi,VECTOR<T,3>(0,0,1));
    solids_parameters.rigid_body_parameters.list(id)->is_static=true;

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies/ground",(T)domain_scale*.005);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction(mu);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Restitution(epsilon);
    solids_parameters.rigid_body_parameters.list(id)->position=VECTOR<T,3>(fluids_parameters.grid.xmin,center.y,center.z);
    solids_parameters.rigid_body_parameters.list(id)->orientation=QUATERNION<T>(-.5*pi,VECTOR<T,3>(0,0,1));
    solids_parameters.rigid_body_parameters.list(id)->is_static=true;

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies/ground",(T)domain_scale*.005);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction(mu);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Restitution(epsilon);
    solids_parameters.rigid_body_parameters.list(id)->position=VECTOR<T,3>(center.x,center.y,fluids_parameters.grid.zmin);
    solids_parameters.rigid_body_parameters.list(id)->orientation=QUATERNION<T>(.5*pi,VECTOR<T,3>(1,0,0));
    solids_parameters.rigid_body_parameters.list(id)->is_static=true;

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies/ground",(T)domain_scale*.005);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction(mu);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Restitution(epsilon);
    solids_parameters.rigid_body_parameters.list(id)->position=VECTOR<T,3>(center.x,center.y,fluids_parameters.grid.zmax);
    solids_parameters.rigid_body_parameters.list(id)->orientation=QUATERNION<T>(-.5*pi,VECTOR<T,3>(1,0,0));
    solids_parameters.rigid_body_parameters.list(id)->is_static=true;
}
//#####################################################################
// Function Construct_Levelsets_For_Objects
//#####################################################################
void Construct_Levelsets_For_Objects(const T time)
{
    GRID<TV>& grid=fluids_parameters.grid;
    phi_objects.Resize(grid,3);ARRAY<T,VECTOR<int,3> >::copy(5*grid.max_dx_dy_dz,phi_objects);
    id_objects.Resize(grid,3);ARRAY<int,VECTOR<int,3> >::copy(0,id_objects);
    for(int object=1;object<=bodies_in_fluid.m;object++){
        RIGID_BODY<TV>& rigid_body=*bodies_in_fluid(object);
        for(int i=-2;i<=grid.m+3;i++)for(int j=-2;j<=grid.n+3;j++)for(int ij=-2;ij<=grid.mn+3;ij++){
            T phi_object=rigid_body.Implicit_Surface_Extended_Value(grid.X(i,j,ij));
            if(phi_object<phi_objects(i,j,ij)){
                id_objects(i,j,ij)=object;phi_objects(i,j,ij)=phi_object;}}}

    // Make signed distance
    FAST_LEVELSET_3D<GRID<TV> > levelset_objects(grid,phi_objects);levelset_objects.Set_Band_Width(3);levelset_objects.Reinitialize();

    // use effective velocity whenever possible
    bool use_effective_velocity=false; if(bodies_in_fluid.m&&bodies_in_fluid(1)->saved_states.m>=2) use_effective_velocity=true;
    if(!use_effective_velocity) std::cout << "Construct_Levelsets_For_Objects: not using effectively velocity (possibly start of simulation)" << std::endl;
    Construct_Velocities_For_Objects(use_effective_velocity,false,false,0,time);
}
//#####################################################################
// Function Construct_Velocities_For_Objects
//#####################################################################
void Construct_Velocities_For_Objects(const bool use_effective_velocity,const bool subtract_last_pressure_force,const bool add_gravity,const T dt,const T time)
{
    T dt_for_velocity_adjustment=use_next_dt_for_velocity_adjustment?fluids_parameters.next_dt:dt;

    VECTOR<T,3> extra_velocity;
    if(add_gravity){
        if(use_next_dt_for_velocity_adjustment) std::cout << "Using next dt " << fluids_parameters.next_dt << " for gravity" << std::endl;
        extra_velocity=dt_for_velocity_adjustment*solids_parameters.gravity*solids_parameters.gravity_direction;
    }

    ARRAY<RIGID_BODY_STATE_3D<T> > state_for_object_velocity(bodies_in_fluid.m);
    if(subtract_last_pressure_force){
        if(use_next_dt_for_velocity_adjustment) std::cout << "Using next dt " << fluids_parameters.next_dt << " to subtract last pressure force" << std::endl;
        for(int i=1;i<=bodies_in_fluid.m;i++){
            assert(&rigid_body_fluid_forces(i)->rigid_body==bodies_in_fluid(i));
            state_for_object_velocity(i)=*bodies_in_fluid(i)->saved_states(COLLISION_BODY_3D<T>::THIN_SHELLS_NEW_STATE);
            std::cout << i << ":  VELOCITY BEFORE " << state_for_object_velocity(i).velocity << std::endl;
            state_for_object_velocity(i).velocity-=dt_for_velocity_adjustment*rigid_body_fluid_forces(i)->last_computed_force/bodies_in_fluid(i)->mass;
            state_for_object_velocity(i).angular_momentum-=dt_for_velocity_adjustment*rigid_body_fluid_forces(i)->last_computed_torque;
            bodies_in_fluid(i)->Update_Angular_Velocity(state_for_object_velocity(i));
            std::cout << i << ":  VELOCITY AFTER " << state_for_object_velocity(i).velocity << std::endl;}}

    GRID<TV>& grid=fluids_parameters.grid;
    V_objects.Resize(grid,3);ARRAY<VECTOR<T,3> ,VECTOR<int,3> >::copy(VECTOR<T,3>(0,0,0),V_objects);
    for(int i=-2;i<=grid.m+3;i++)for(int j=-2;j<=grid.n+3;j++)for(int ij=-2;ij<=grid.mn+3;ij++)if(phi_objects(i,j,ij)<0){
        int body_id=id_objects(i,j,ij);
        if(use_effective_velocity) V_objects(i,j,ij)=bodies_in_fluid(body_id)->Pointwise_Object_Pseudo_Velocity(0,grid.X(i,j,ij),
                                                    COLLISION_BODY_3D<T>::THIN_SHELLS_OLD_STATE,COLLISION_BODY_3D<T>::THIN_SHELLS_NEW_STATE);
        else if(subtract_last_pressure_force) V_objects(i,j,ij)=bodies_in_fluid(body_id)->Pointwise_Object_Velocity(grid.X(i,j,ij),state_for_object_velocity(body_id));
        else V_objects(i,j,ij)=extra_velocity+bodies_in_fluid(body_id)->Pointwise_Object_Velocity(grid.X(i,j,ij));}
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    // Not so good to set up a heaviside function here because then the interface will
    // be exactly between the two nodes which can lead to roundoff issues when setting dirichlet cells, etc.
    GRID<TV>& grid=fluids_parameters.grid;
    for(int i=1;i<=grid.m;i++) for(int j=1;j<=grid.n;j++) for(int ij=1;ij<=grid.mn;ij++){
        fluids_parameters.particle_levelset_evolution.phi(i,j,ij)=grid.y(j)-grid.ymin-initial_water_level;
        if(phi_objects(i,j,ij)<0) fluids_parameters.particle_levelset_evolution.phi(i,j,ij)=-phi_objects(i,j,ij);
    }
}
//#####################################################################
// Utilities
//#####################################################################
bool Mixed_X_Face(const int i,const int j,const int ij){return phi_objects(i,j,ij)<0 || phi_objects(i,j+1,ij)<0 || phi_objects(i,j,ij+1)<0 || phi_objects(i,j+1,ij+1)<0;}
bool Mixed_Y_Face(const int i,const int j,const int ij){return phi_objects(i,j,ij)<0 || phi_objects(i+1,j,ij)<0 || phi_objects(i,j,ij+1)<0 || phi_objects(i+1,j,ij+1)<0;}
bool Mixed_Z_Face(const int i,const int j,const int ij){return phi_objects(i,j,ij)<0 || phi_objects(i+1,j,ij)<0 || phi_objects(i,j+1,ij)<0 || phi_objects(i+1,j+1,ij)<0;}
//#####################################################################
// Function Average_Node_Velocities_To_Faces
//#####################################################################
void Average_Node_Velocities_To_Faces()
{
    if(set_faces_to_invalid){
        std::cout << "Setting psi_D_old to all true" << std::endl;
        ARRAY<bool,VECTOR<int,3> >::copy(true,fluids_parameters.incompressible.psi_D_old);
    }

    if(!set_fluid_boundary_conditions_for_pressure_jump_solve){
        // Average from nodes to faces but ignore nodes inside objects
        // Note: if all nodes on face are in object the face must be psi_N so we don't process it

        PROJECTION_3D<T>& projection=fluids_parameters.incompressible.projection;
        ARRAY<VECTOR<T,3> ,VECTOR<int,3> > &V=fluids_parameters.incompressible.V,&V_ghost=fluids_parameters.incompressible.V_ghost;
        ARRAY<T,VECTOR<int,3> > &u=projection.u,&v=projection.v,&w=projection.w;
        ARRAY<bool,VECTOR<int,3> > &psi_N_u_face=projection.elliptic_solver->psi_N_u,&psi_N_v_face=projection.elliptic_solver->psi_N_v,&psi_N_w_face=projection.elliptic_solver->psi_N_w; 

        static int index_lookup[5][2]={{0,0},{0,0},{1,0},{0,1},{1,1}};
        static T one_over_count_lookup[5]={0,1,(T).5,(T)one_third,(T).25};

        for(int i=1;i<=projection.u_grid.m;i++) for(int j=1;j<=projection.u_grid.n;j++) for(int ij=1;ij<=projection.u_grid.mn;ij++) if(!psi_N_u_face(i,j,ij)){
            int count=0;T velocity=0;bool mixed_face=Mixed_X_Face(i,j,ij);
            if(only_use_object_velocity_in_mixed_face && mixed_face){
                for(int c=1;c<=4;c++){VECTOR<int,3> index(i,j+index_lookup[c][0],ij+index_lookup[c][1]);
                    if(phi_objects(index)<0){velocity+=V(index).x;count++;}}
                assert(count);u(i,j,ij)=velocity*one_over_count_lookup[count];}
            else if(average_velocities_in_mixed_face && mixed_face){
                for(int c=1;c<=4;c++){VECTOR<int,3> index(i,j+index_lookup[c][0],ij+index_lookup[c][1]);velocity+=V(index).x;}
                u(i,j,ij)=(T).25*velocity;}
            else if(!fluids_parameters.incompressible.psi_D_old(i-1,j,ij)||!fluids_parameters.incompressible.psi_D_old(i,j,ij)){
                for(int c=1;c<=4;c++){VECTOR<int,3> index(i,j+index_lookup[c][0],ij+index_lookup[c][1]);
                    if(phi_objects(index)>=0){velocity+=V(index).x-V_ghost(index).x;count++;}}
                assert(count);u(i,j,ij)+=velocity*one_over_count_lookup[count];}
            else{
                for(int c=1;c<=4;c++){VECTOR<int,3> index(i,j+index_lookup[c][0],ij+index_lookup[c][1]);
                    if(phi_objects(index)>=0){velocity+=V(index).x;count++;}}
                assert(count);u(i,j,ij)=velocity*one_over_count_lookup[count];}
        }
            
        for(int i=1;i<=projection.v_grid.m;i++) for(int j=1;j<=projection.v_grid.n;j++) for(int ij=1;ij<=projection.v_grid.mn;ij++) if(!psi_N_v_face(i,j,ij)){
            int count=0;T velocity=0;bool mixed_face=Mixed_Y_Face(i,j,ij);
            if(only_use_object_velocity_in_mixed_face && mixed_face){
                for(int c=1;c<=4;c++){VECTOR<int,3> index(i+index_lookup[c][0],j,ij+index_lookup[c][1]);
                    if(phi_objects(index)<0){velocity+=V(index).y;count++;}}
                assert(count);v(i,j,ij)=velocity*one_over_count_lookup[count];}
            else if(average_velocities_in_mixed_face && mixed_face){
                for(int c=1;c<=4;c++){VECTOR<int,3> index(i+index_lookup[c][0],j,ij+index_lookup[c][1]);velocity+=V(index).y;}
                v(i,j,ij)=(T).25*velocity;}
            else if(!fluids_parameters.incompressible.psi_D_old(i,j-1,ij)||!fluids_parameters.incompressible.psi_D_old(i,j,ij)){
                for(int c=1;c<=4;c++){VECTOR<int,3> index(i+index_lookup[c][0],j,ij+index_lookup[c][1]);
                    if(phi_objects(index)>=0){velocity+=V(index).y-V_ghost(index).y;count++;}}
                assert(count);v(i,j,ij)+=velocity*one_over_count_lookup[count];}
            else{
                for(int c=1;c<=4;c++){VECTOR<int,3> index(i+index_lookup[c][0],j,ij+index_lookup[c][1]);
                    if(phi_objects(index)>=0){velocity+=V(index).y;count++;}}
                assert(count);v(i,j,ij)=velocity*one_over_count_lookup[count];}
        }

        for(int i=1;i<=projection.w_grid.m;i++) for(int j=1;j<=projection.w_grid.n;j++) for(int ij=1;ij<=projection.w_grid.mn;ij++) if(!psi_N_w_face(i,j,ij)){
            int count=0;T velocity=0;bool mixed_face=Mixed_Z_Face(i,j,ij);
            if(only_use_object_velocity_in_mixed_face && mixed_face){
                for(int c=1;c<=4;c++){VECTOR<int,3> index(i,j+index_lookup[c][0],ij+index_lookup[c][1]);
                    if(phi_objects(index)<0){velocity+=V(index).z;count++;}}
                assert(count);w(i,j,ij)=velocity*one_over_count_lookup[count];}
            else if(average_velocities_in_mixed_face && mixed_face){
                for(int c=1;c<=4;c++){VECTOR<int,3> index(i,j+index_lookup[c][0],ij+index_lookup[c][1]);velocity+=V(index).z;}
                w(i,j,ij)=(T).25*velocity;}
            else if(!fluids_parameters.incompressible.psi_D_old(i,j,ij-1)||!fluids_parameters.incompressible.psi_D_old(i,j,ij)){
                for(int c=1;c<=4;c++){VECTOR<int,3> index(i+index_lookup[c][0],j+index_lookup[c][1],ij);
                    if(phi_objects(index)>=0){velocity+=V(index).z-V_ghost(index).z;count++;}}
                assert(count);w(i,j,ij)+=velocity*one_over_count_lookup[count];}
            else{
                for(int c=1;c<=4;c++){VECTOR<int,3> index(i+index_lookup[c][0],j+index_lookup[c][1],ij);
                    if(phi_objects(index)>=0){velocity+=V(index).z;count++;}}
                assert(count);w(i,j,ij)=velocity*one_over_count_lookup[count];}
        }
    }
    else WATER_MELTING_EXAMPLE_3D<T,RW>::Average_Node_Velocities_To_Faces();
}
//#####################################################################
// Function Adjust_Velocity_With_Objects
//#####################################################################
void Adjust_Velocity_With_Objects(const T time)
{
    // V_objects should be set with the effective velocity at this point (and phi_objects is in object's initial position for this time step)
    fluids_parameters.Extrapolate_Velocity_Into_Object(phi_objects,V_objects,velocity_extrapolation_band_width,false,time);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time) PHYSBAM_OVERRIDE
{
    if(!set_fluid_boundary_conditions_for_pressure_jump_solve){
        GRID<TV> &grid=fluids_parameters.grid,&u_grid=fluids_parameters.u_grid,&v_grid=fluids_parameters.v_grid,&w_grid=fluids_parameters.w_grid;
        LAPLACE_3D<T>& elliptic_solver=*fluids_parameters.incompressible.projection.elliptic_solver;

        Construct_Velocities_For_Objects(use_object_pseudo_velocity_for_boundary_conditions,false,false,dt,time);

        ARRAY<VECTOR<T,3> ,VECTOR<int,3> >& V=fluids_parameters.incompressible.V;
        if(use_extrapolated_velocity_for_laplace) fluids_parameters.Extrapolate_Velocity_Into_Object(phi_objects,V_objects,velocity_extrapolation_band_width,false,time);
        else for(int i=1;i<=grid.m;i++) for(int j=1;j<=grid.n;j++) for(int ij=1;ij<=grid.mn;ij++) if(phi_objects(i,j,ij)<0) V(i,j,ij)=V_objects(i,j,ij);

        static int index_lookup[5][2]={{0,0},{0,0},{1,0},{0,1},{1,1}};
        static T one_over_count_lookup[5]={0,1,(T).5,(T)one_third,(T).25};

        // Get psi_N velocities by selectively averaging only from object when at least one node is inside object
        for(int i=1;i<=u_grid.m;i++) for(int j=1;j<=u_grid.n;j++) for(int ij=1;ij<=u_grid.mn;ij++) 
            if((!mixed_face_is_solid && phi_objects(i,j,ij)+phi_objects(i,j+1,ij)+phi_objects(i,j,ij+1)+phi_objects(i,j+1,ij+1)<0) || (mixed_face_is_solid && Mixed_X_Face(i,j,ij))){
                int count=0;T velocity=0;
                for(int c=1;c<=4;c++){
                    VECTOR<int,3> index(i,j+index_lookup[c][0],ij+index_lookup[c][1]);
                    if(phi_objects(index)<0){velocity+=V_objects(index).x;count++;}}
                velocity*=one_over_count_lookup[count];
                elliptic_solver.psi_N_u(i,j,ij)=true;fluids_parameters.incompressible.projection.u(i,j,ij)=velocity;}

        for(int i=1;i<=v_grid.m;i++) for(int j=1;j<=v_grid.n;j++) for(int ij=1;ij<=v_grid.mn;ij++) 
            if((!mixed_face_is_solid && phi_objects(i,j,ij)+phi_objects(i+1,j,ij)+phi_objects(i,j,ij+1)+phi_objects(i+1,j,ij+1)<0) || (mixed_face_is_solid && Mixed_Y_Face(i,j,ij))){
                int count=0;T velocity=0;
                for(int c=1;c<=4;c++){
                    VECTOR<int,3> index(i+index_lookup[c][0],j,ij+index_lookup[c][1]);
                    if(phi_objects(index)<0){velocity+=V_objects(index).y;count++;}}
                velocity*=one_over_count_lookup[count];
                elliptic_solver.psi_N_v(i,j,ij)=true;fluids_parameters.incompressible.projection.v(i,j,ij)=velocity;}

        for(int i=1;i<=w_grid.m;i++) for(int j=1;j<=w_grid.n;j++) for(int ij=1;ij<=w_grid.mn;ij++) 
            if((!mixed_face_is_solid && phi_objects(i,j,ij)+phi_objects(i,j+1,ij)+phi_objects(i+1,j,ij)+phi_objects(i+1,j+1,ij)<0) || (mixed_face_is_solid && Mixed_Z_Face(i,j,ij))){
                int count=0;T velocity=0;
                for(int c=1;c<=4;c++){
                    VECTOR<int,3> index(i+index_lookup[c][0],j+index_lookup[c][1],ij);
                    if(phi_objects(index)<0){velocity+=V_objects(index).z;count++;}}
                velocity*=one_over_count_lookup[count];
                elliptic_solver.psi_N_w(i,j,ij)=true;fluids_parameters.incompressible.projection.w(i,j,ij)=velocity;}
    }
    else{
        GRID<TV> &grid=fluids_parameters.grid,&u_grid=fluids_parameters.u_grid,&v_grid=fluids_parameters.v_grid,&w_grid=fluids_parameters.w_grid;
        POISSON_3D<T>* poisson=fluids_parameters.incompressible.projection.poisson;
        assert(poisson && need_to_compute_pressure_coupling && (fluids_parameters.fire || poisson->use_variable_beta));
        ARRAY<T,VECTOR<int,3> >::copy((T)1/fluids_parameters.density,poisson->beta_right);ARRAY<T,VECTOR<int,3> >::copy((T)1/fluids_parameters.density,poisson->beta_top);ARRAY<T,VECTOR<int,3> >::copy((T)1/fluids_parameters.density,poisson->beta_back);
        poisson->Use_Body_Psi_N(true);
        ARRAY<bool,VECTOR<int,3> >::copy(false,poisson->body_psi_N_u);ARRAY<bool,VECTOR<int,3> >::copy(false,poisson->body_psi_N_v);ARRAY<bool,VECTOR<int,3> >::copy(false,poisson->body_psi_N_w);

        Construct_Velocities_For_Objects(use_effective_velocity_for_poisson,subtract_last_pressure_force_from_object,add_gravity_to_object,dt,time);

        if(extrapolate_combined_solid_fluid_velocity){ // extrapolate combined fluid and solid velocity into air
            ARRAY<T,VECTOR<int,3> > combined_solid_fluid_phi(grid);
            FAST_LEVELSET_3D<GRID<TV> > combined_solid_fluid_levelset(grid,combined_solid_fluid_phi);
            combined_solid_fluid_levelset.Set_Band_Width(3);combined_solid_fluid_levelset.Reinitialize();
            for(int i=1;i<=grid.m;i++) for(int j=1;j<=grid.n;j++) for(int ij=1;ij<=grid.mn;ij++){
                combined_solid_fluid_phi(i,j,ij)=min(phi_objects(i,j,ij),fluids_parameters.particle_levelset_evolution.phi(i,j,ij));
                if(phi_objects(i,j,ij)<0) fluids_parameters.incompressible.V(i,j,ij)=V_objects(i,j,ij);}
            EXTRAPOLATION_3D<T,VECTOR<T,3> > extrapolate(grid,combined_solid_fluid_phi,fluids_parameters.incompressible.V);extrapolate.Set_Band_Width(3);extrapolate.Extrapolate();}

        static int index_lookup[5][2]={{0,0},{0,0},{1,0},{0,1},{1,1}};
        static T one_over_count_lookup[5]={0,1,(T).5,(T)one_third,(T).25};

        for(int i=1;i<=u_grid.m;i++) for(int j=1;j<=u_grid.n;j++) for(int ij=1;ij<=u_grid.mn;ij++)
            if((!mixed_face_is_solid && phi_objects(i,j,ij)+phi_objects(i,j+1,ij)+phi_objects(i,j,ij+1)+phi_objects(i,j+1,ij+1)<0) || (mixed_face_is_solid && Mixed_X_Face(i,j,ij))){
                int count=0;T density=0,velocity=0;
                for(int c=1;c<=4;c++){
                    VECTOR<int,3> index(i,j+index_lookup[c][0],ij+index_lookup[c][1]);
                    if(phi_objects(index)<0 && collision_body_affected_by_fluid(id_objects(index))){int body_id=id_objects(index);
                        density+=body_densities(body_id);velocity+=V_objects(index).x;count++;}}
                density*=one_over_count_lookup[count];velocity*=one_over_count_lookup[count];
                poisson->psi_N_u(i,j,ij)=true;fluids_parameters.incompressible.projection.u(i,j,ij)=velocity;
                poisson->body_psi_N_u(i,j,ij)=true;poisson->beta_right(i-1,j,ij)=1/density;}

        for(int i=1;i<=v_grid.m;i++) for(int j=1;j<=v_grid.n;j++) for(int ij=1;ij<=v_grid.mn;ij++)
            if((!mixed_face_is_solid && phi_objects(i,j,ij)+phi_objects(i+1,j,ij)+phi_objects(i,j,ij+1)+phi_objects(i+1,j,ij+1)<0) || (mixed_face_is_solid && Mixed_Y_Face(i,j,ij))){
                int count=0;T density=0,velocity=0;
                for(int c=1;c<=4;c++){
                    VECTOR<int,3> index(i+index_lookup[c][0],j,ij+index_lookup[c][1]);
                    if(phi_objects(index)<0 && collision_body_affected_by_fluid(id_objects(index))){int body_id=id_objects(index);
                        density+=body_densities(body_id);velocity+=V_objects(index).y;count++;}}
                density*=one_over_count_lookup[count];velocity*=one_over_count_lookup[count];
                poisson->psi_N_v(i,j,ij)=true;fluids_parameters.incompressible.projection.v(i,j,ij)=velocity;
                poisson->body_psi_N_v(i,j,ij)=true;poisson->beta_top(i,j-1,ij)=1/density;}

        for(int i=1;i<=w_grid.m;i++) for(int j=1;j<=w_grid.n;j++) for(int ij=1;ij<=w_grid.mn;ij++)
            if((!mixed_face_is_solid && phi_objects(i,j,ij)+phi_objects(i,j+1,ij)+phi_objects(i+1,j,ij)+phi_objects(i+1,j+1,ij)<0) || (mixed_face_is_solid && Mixed_Z_Face(i,j,ij))){
                int count=0;T density=0,velocity=0;
                for(int c=1;c<=4;c++){
                    VECTOR<int,3> index(i+index_lookup[c][0],j+index_lookup[c][1],ij);
                    if(phi_objects(index)<0 && collision_body_affected_by_fluid(id_objects(index))){int body_id=id_objects(index);
                        density+=body_densities(body_id);velocity+=V_objects(index).z;count++;}}
                density*=one_over_count_lookup[count];velocity*=one_over_count_lookup[count];
                poisson->psi_N_w(i,j,ij)=true;fluids_parameters.incompressible.projection.w(i,j,ij)=velocity;
                poisson->body_psi_N_w(i,j,ij)=true;poisson->beta_back(i,j,ij-1)=1/density;}}
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
void Set_Dirichlet_Boundary_Conditions(const T time)
{
    if(!set_fluid_boundary_conditions_for_pressure_jump_solve){
        WATER_MELTING_EXAMPLE_3D<T,RW>::Set_Dirichlet_Boundary_Conditions(time);
        return;}

    GRID<TV>& p_grid=fluids_parameters.incompressible.projection.p_grid;
    ARRAY<T,VECTOR<int,3> >& phi=fluids_parameters.particle_levelset_evolution.phi;
    for(int i=1;i<=p_grid.m;i++) for(int j=1;j<=p_grid.n;j++) for(int ij=1;ij<=p_grid.mn;ij++) 
        if((phi(i,j,ij)+phi(i+1,j,ij)+phi(i,j+1,ij)+phi(i+1,j+1,ij)+phi(i,j,ij+1)+phi(i+1,j,ij+1)+phi(i,j+1,ij+1)+phi(i+1,j+1,ij+1) > 0) &&
           (phi_objects(i,j,ij)+phi_objects(i+1,j,ij)+phi_objects(i,j+1,ij)+phi_objects(i+1,j+1,ij)+phi_objects(i,j,ij+1)+phi_objects(i+1,j,ij+1)+phi_objects(i,j+1,ij+1)+phi_objects(i+1,j+1,ij+1) > 0)){
            fluids_parameters.incompressible.projection.elliptic_solver->psi_D(i,j,ij)=true;fluids_parameters.incompressible.projection.p(i,j,ij)=0;}
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
void Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<T,VECTOR<T,3> >& particles,const typename PARTICLE_LEVELSET<T,VECTOR<T,3> >::PARTICLE_TYPE particle_type,const T time)
{
    LEVELSET_3D<GRID<TV> > levelset_objects(fluids_parameters.grid,phi_objects);
    if(particle_type == PARTICLE_LEVELSET<T,VECTOR<T,3> >::NEGATIVE || particle_type == PARTICLE_LEVELSET<T,VECTOR<T,3> >::REMOVED_NEGATIVE){
        for(int k=particles.array_collection->Size();k>=1;k--) if(levelset_objects.Extended_Phi(particles.X(k))<-fluids_parameters.grid.min_dx_dy_dz) particles.Delete_Particle(k);}
   else for(int k=particles.array_collection->Size();k>=1;k--) if(levelset_objects.Extended_Phi(particles.X(k))<0) particles.Delete_Particle(k);
}
//#####################################################################
// Function Adjust_Particle_For_Objects
//#####################################################################
void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE
{
    const LEVELSET_3D<GRID<TV> > levelset_objects((GRID_3D<T>&)fluids_parameters.grid,(ARRAY<T,VECTOR<int,3> >&)phi_objects);
    FILE_UTILITIES::Write_To_File<RW>(STRING_UTILITIES::string_sprintf("%s/phi_objects.%d",output_directory.c_str(),frame),levelset_objects);
    ARRAY<PAIR<VECTOR<T,3>,VECTOR<T,3> > > rigid_body_forces_and_torques(solids_parameters.rigid_body_parameters.list.Number_Of_Elements());
    for(int i=1;i<=solids_parameters.rigid_body_parameters.list.Number_Of_Elements();i++)if(solids_parameters.rigid_body_parameters.list.Is_Active(i)){
        const RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list(i);
        if(rigid_body.fluid_forces.m){rigid_body_forces_and_torques(i).x=rigid_body.fluid_forces(1)->last_computed_force;rigid_body_forces_and_torques(i).y=rigid_body.fluid_forces(1)->last_computed_torque;}}
    FILE_UTILITIES::Write_To_File<RW>(STRING_UTILITIES::string_sprintf("%s/rigid_body_forces_and_torques.%d",output_directory.c_str(),frame),rigid_body_forces_and_torques);
    WATER_MELTING_EXAMPLE_3D<T,RW>::Write_Output_Files(frame);
}
//#####################################################################
};    
}
#endif  



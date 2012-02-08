//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class THIN_SHELLS_TEST_WATER 
//##################################################################### 
#ifndef __THIN_SHELLS_TEST_WATER__
#define __THIN_SHELLS_TEST_WATER__

#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_2D.h>
namespace PhysBAM{

template<class T,class RW=T>
class THIN_SHELLS_TEST_WATER:public SOLIDS_FLUIDS_EXAMPLE_2D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::solids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::verbose_dt;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::write_time;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::write_frame_title;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::data_directory;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::abort_when_dt_below;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Use_Thin_Shells_Fluid_Coupling_Defaults;

    int example_number;
    bool use_source;
    BOX_2D<T> source;
    MATRIX<T,3> world_to_source;
    VECTOR_2D<T> source_velocity;
    T initial_water_level;
    bool really_use_effective_velocities;

    PARAMETER_LIST parameter_list;
    bool use_source_override;
    T domain_scale;

    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return false;}

    THIN_SHELLS_TEST_WATER(const int example_number_input)
        :SOLIDS_FLUIDS_EXAMPLE_2D<RW>(fluids_parameters.WATER),example_number(example_number_input),domain_scale(1)
    {
        std::cout << "Example number " << example_number << std::endl;
        {std::string filename=STRING_UTILITIES::string_sprintf("Thin_Shells_Test_Water/example_%d.param",example_number);
        if(FILE_UTILITIES::File_Exists(filename)){std::cout << "Reading parameter file '" << filename << "'" << std::endl;parameter_list.Read(filename);}}

        // Two way coupling
        fluids_parameters.use_solid_fluid_coupling=true;
        Set_Parameter(fluids_parameters.use_solid_fluid_coupling,"fluids_parameters.use_solid_fluid_coupling");
        if(fluids_parameters.use_solid_fluid_coupling) Use_Thin_Shells_Fluid_Coupling_Defaults();

        Set_Parameter(domain_scale,"domain_scale");

        // Common parameters
        output_directory="Thin_Shells_Test_Water/output";
        first_frame=0;last_frame=1000;frame_rate=60;
        restart=false;restart_frame=0;

        // Fluids parameters
        fluids_parameters.grid.Initialize(101,101,0,domain_scale,0,domain_scale);
        fluids_parameters.number_particles_per_cell=16;
        fluids_parameters.reseeding_frame_rate=10;
        fluids_parameters.incompressible_iterations=100;
        fluids_parameters.cfl=0.9;
        fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=false;
        fluids_parameters.bias_towards_negative_particles=true;
        fluids_parameters.use_removed_positive_particles=false;fluids_parameters.use_removed_negative_particles=true;
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;fluids_parameters.variable_viscosity=false;fluids_parameters.second_order_cut_cell_method=false;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;
        fluids_parameters.write_particles=true;fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=true;
        fluids_parameters.store_particle_ids=true;
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.enforce_divergence_free_extrapolation=false;
        fluids_parameters.use_old_velocities_for_boundary_conditions=false;

        // Solids parameters
        solids_parameters.gravity=9.8;
        solids_parameters.cfl=(T).9;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies=true;
        solids_parameters.deformable_body_parameters.write=true;
        verbose_dt=true;

        if(example_number==1){
            //initial_water_level=(T)1.021; // some number so the interface won't lie in the center of a cell
            initial_water_level=-1; // no initial height
            use_source=true;
            source=BOX_2D<T>(0.375,0.575,0.8,1);
            world_to_source=MATRIX<T,3>::Identity_Matrix();
            source_velocity=VECTOR_2D<T>(0,-0.16);}
        else if(example_number==2){
            initial_water_level=0.3;
            use_source=false;
            fluids_parameters.simulate=true;
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;}
        else if(example_number==3){
            initial_water_level=-1; // no initial height
            use_source=true;
            source=BOX_2D<T>(0,.1,0.75,0.85);
            world_to_source=MATRIX<T,3>::Identity_Matrix();
            source_velocity=VECTOR_2D<T>((T)2,0);}
        else if(example_number==4){
            initial_water_level=0.3;
            use_source=false;
            fluids_parameters.simulate=true;
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;}
        else if(example_number==5){
            initial_water_level=0.3;
            use_source=false;
            fluids_parameters.simulate=true;
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;}
        else if(example_number==6){
            initial_water_level=0.7;
            use_source=false;
            fluids_parameters.simulate=true;
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;}
        else if(example_number==7){
            initial_water_level=domain_scale*0.6;
            fluids_parameters.simulate=true;
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
            use_source=false;
            source=BOX_2D<T>(domain_scale*.5,domain_scale*.6,domain_scale*0.9,domain_scale*1);
            world_to_source=MATRIX<T,3>::Identity_Matrix();
            source_velocity=VECTOR_2D<T>(0,(T)-domain_scale*1);}
        else if(example_number==8){
            initial_water_level=0.6;
            fluids_parameters.simulate=true;
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
            use_source=false;}
        else if(example_number==9){
            initial_water_level=-1;
            fluids_parameters.simulate=true;
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
            PHYSBAM_FATAL_ERROR("Preroll is gone; this example is likely broken.");
            use_source=true;
            source=BOX_2D<T>(.35-.075,.35+.075,.9,1);
            world_to_source=MATRIX<T,3>::Identity_Matrix();
            source_velocity=VECTOR_2D<T>(0,(T)-1);}
        else if(example_number==10){
            initial_water_level=0.6;
            use_source=false;}
        else if(example_number==11){ // test floating square
            initial_water_level=domain_scale*0.6;
            fluids_parameters.simulate=true;
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
            fluids_parameters.number_particles_per_cell=0;
            fluids_parameters.use_time_filter_for_pressure_jumps=true;
            use_source=false;}
        else if(example_number==12){ // test buoyancy
            initial_water_level=domain_scale*0.6;
            fluids_parameters.simulate=true;
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
            use_source=false;}

        // Debugging
        fluids_parameters.write_debug_data=true;
        abort_when_dt_below=1e-7;
        write_time=true;
        write_frame_title=true;

        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Parameters_From_Parameter_List(*this,parameter_list);
        use_source_override=true;Set_Parameter(use_source_override,"use_source");
        Set_Parameter(source_velocity,"source_velocity");
        really_use_effective_velocities=true;
        Set_Parameter(really_use_effective_velocities,"really_use_effective_velocities");
    }

    ~THIN_SHELLS_TEST_WATER() 
    {}

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
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    int id;

    if(example_number==1){
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies_2D/ground",(T)1/360);
        RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list(id);
        rigid_body.segmented_curve->Set_Density(1);rigid_body.Set_Mass(rigid_body.segmented_curve->Thin_Shell_Mass());
        rigid_body.frame.t=VECTOR_2D<T>(0.5,0.65);
        rigid_body.frame.r=COMPLEX<T>::Unit_Polar(-0.4);
        rigid_body.is_static=true;
        Add_To_Fluid_Simulation(rigid_body,true,false,100);}
    else if(example_number==2){
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies_2D/u",(T)1/3);
        solids_parameters.rigid_body_parameters.list(id)->is_kinematic=true;
        Add_To_Fluid_Simulation(*solids_parameters.rigid_body_parameters.list(id),true,false,100);}
    else if(example_number==3){
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(solids_parameters.deformable_body_parameters.list,30,VECTOR_2D<T>((T).51003,(T).8),VECTOR_2D<T>((T).51003,(T).2));
        T density=1;T pressure_force_scale=1;THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Deformable_Object_Parameters_2D(id,density,pressure_force_scale,parameter_list);
        DEFORMABLE_OBJECT_2D<T>* deformable_object=solids_parameters.deformable_body_parameters.list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Density(*deformable_object->segmented_curve,density);
        deformable_object->Add_Body_Forces(*deformable_object->segmented_curve,solids_parameters.gravity,solids_parameters.gravity_direction);
        deformable_object->Add_Edge_Springs(deformable_object->segmented_curve->segment_mesh,300,(T).9);
        Add_To_Fluid_Simulation(*deformable_object,true,true,pressure_force_scale);}
    else if(example_number==4){
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies_2D/cup",(T)1/3);
        solids_parameters.rigid_body_parameters.list(id)->is_kinematic=true;
        Add_To_Fluid_Simulation(*solids_parameters.rigid_body_parameters.list(id),true,false,100);}
    else if(example_number==5){
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies_2D/cup",(T)1/3);
        solids_parameters.rigid_body_parameters.list(id)->is_kinematic=true;
        Add_To_Fluid_Simulation(*solids_parameters.rigid_body_parameters.list(id),true,false,100);}
    else if(example_number==6){
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies_2D/boat",(T)1/3);
        solids_parameters.rigid_body_parameters.list(id)->frame.t=VECTOR_2D<T>(.5,.4);
        solids_parameters.rigid_body_parameters.list(id)->frame.r=COMPLEX<T>::Unit_Polar(.2);
        solids_parameters.rigid_body_parameters.list(id)->Set_Mass(1);
        solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
        Add_To_Fluid_Simulation(*solids_parameters.rigid_body_parameters.list(id),true,true,6000);}
    else if(example_number==7){
        T mu=0.9,epsilon=0.2;
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies_2D/boat",(T)domain_scale/3);
//        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies_2D/ground_subdivided",(T).01);
        T density=1;T pressure_force_scale=5.75;THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Rigid_Body_Parameters_2D(id,density,pressure_force_scale,parameter_list);
        RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list(id);
        rigid_body.Set_Coefficient_Of_Friction(mu);
        rigid_body.Set_Coefficient_Of_Restitution(epsilon);
        rigid_body.segmented_curve->Set_Density(density);rigid_body.Set_Mass(rigid_body.segmented_curve->Thin_Shell_Mass());
        rigid_body.frame.t=VECTOR_2D<T>(domain_scale*.5,domain_scale*.65);
        //rigid_body.frame.t=VECTOR_2D<T>(domain_scale*.5,domain_scale*.35);
        //rigid_body.frame.r=COMPLEX<T>::Unit_Polar(pi/4);
        rigid_body.Update_Angular_Momentum();
        rigid_body.Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
        Add_To_Fluid_Simulation(rigid_body,true,true,pressure_force_scale);

        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Rigid_Body_Walls(*this,epsilon,mu);

        solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);}
    else if(example_number==8){ // length is .2
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies_2D/ground_subdivided",(T).001);
        RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list(id);
        rigid_body.segmented_curve->Set_Density(5);
        rigid_body.Set_Mass(rigid_body.segmented_curve->Thin_Shell_Mass());
        rigid_body.frame.t=VECTOR_2D<T>(.5,.50);
        rigid_body.Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
        Add_To_Fluid_Simulation(rigid_body,true,true,1);}
    else if(example_number==9){
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(solids_parameters.deformable_body_parameters.list,30,VECTOR_2D<T>((T).1,(T).6),VECTOR_2D<T>((T).7,(T).45));
        DEFORMABLE_OBJECT_2D<T>* deformable_object=solids_parameters.deformable_body_parameters.list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Mass(*deformable_object->segmented_curve,1,true);
        deformable_object->Add_Body_Forces(*deformable_object->segmented_curve,solids_parameters.gravity,solids_parameters.gravity_direction);
        deformable_object->Add_Edge_Springs(deformable_object->segmented_curve->segment_mesh,300,(T).9);
        Add_To_Fluid_Simulation(*deformable_object,true,true,1);}
    else if(example_number==10){
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(solids_parameters.deformable_body_parameters.list,30,VECTOR_2D<T>((T).4,(T).3503),VECTOR_2D<T>((T).6,(T).4503));
        DEFORMABLE_OBJECT_2D<T>* deformable_object=solids_parameters.deformable_body_parameters.list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Mass(*deformable_object->segmented_curve,1,true);
        deformable_object->Add_Body_Forces(*deformable_object->segmented_curve,solids_parameters.gravity,solids_parameters.gravity_direction);
        deformable_object->Add_Edge_Springs(deformable_object->segmented_curve->segment_mesh,300,(T).9);
        Add_To_Fluid_Simulation(*deformable_object,true,true,1);}
    else if(example_number==11){
        T mu=0,epsilon=0.2;

#if 0
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(true,data_directory+"/Rigid_Bodies_2D/Thin_Shells/attached_squares_extra_refined",(T).08123*domain_scale);
        T density=1;T pressure_force_scale=5.75;THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Rigid_Body_Parameters_2D(id,density,pressure_force_scale,parameter_list);
        RIGID_BODY<TV>* rigid_body=solids_parameters.rigid_body_parameters.list(id);
        rigid_body->Set_Coefficient_Of_Friction(mu);
        rigid_body->Set_Coefficient_Of_Restitution(epsilon);
        rigid_body->segmented_curve->Set_Density(density);rigid_body->Set_Mass(rigid_body->segmented_curve->Thin_Shell_Mass());
        //rigid_body->frame.t=VECTOR_2D<T>(domain_scale*.50333,domain_scale*.750333);
        rigid_body->frame.t=VECTOR_2D<T>(domain_scale*.50333,domain_scale*.750333);
        rigid_body->Update_Angular_Momentum();
        rigid_body->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
        Add_To_Fluid_Simulation(*rigid_body,true,true,pressure_force_scale);
#endif
#if 1
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(true,data_directory+"/Rigid_Bodies_2D/Thin_Shells/square_shell_extra_refined",(T).08123*domain_scale);
        T density=1;T pressure_force_scale=5.75;THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Rigid_Body_Parameters_2D(id,density,pressure_force_scale,parameter_list);
        RIGID_BODY<TV>* rigid_body=solids_parameters.rigid_body_parameters.list(id);
        rigid_body->Set_Coefficient_Of_Friction(mu);
        rigid_body->Set_Coefficient_Of_Restitution(epsilon);
        rigid_body->segmented_curve->Set_Density(density);rigid_body->Set_Mass(rigid_body->segmented_curve->Thin_Shell_Mass());
        //rigid_body->frame.t=VECTOR_2D<T>(domain_scale*.50333,domain_scale*.750333);
        //rigid_body->frame.t=VECTOR_2D<T>(domain_scale*.40333,domain_scale*.750333);
        rigid_body->frame.t=VECTOR_2D<T>(domain_scale*.40333,domain_scale*.450333);
        rigid_body->Update_Angular_Momentum();
        rigid_body->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
        Add_To_Fluid_Simulation(*rigid_body,true,true,pressure_force_scale);
#endif

#if 0
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(true,data_directory+"/Rigid_Bodies_2D/Thin_Shells/square_shell_extra_refined",(T).08123*domain_scale);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Rigid_Body_Parameters_2D(id,density,pressure_force_scale,parameter_list);
        rigid_body=solids_parameters.rigid_body_parameters.list(id);
        rigid_body->Set_Coefficient_Of_Friction(mu);
        rigid_body->Set_Coefficient_Of_Restitution(epsilon);
        rigid_body->segmented_curve->Set_Density(density);rigid_body->Set_Mass(rigid_body->segmented_curve->Thin_Shell_Mass());
        //rigid_body->frame.t=VECTOR_2D<T>(domain_scale*.45333,domain_scale*1.00333);
        rigid_body->frame.t=VECTOR_2D<T>(domain_scale*.57333,domain_scale*.750333);
        rigid_body->Update_Angular_Momentum();
        rigid_body->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
        Add_To_Fluid_Simulation(*rigid_body,true,true,pressure_force_scale);
#endif

        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Rigid_Body_Walls(*this,epsilon,mu);

        solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);}

    // Move kinematic bodies to their position at time 0
    solids_parameters.rigid_body_parameters.Reset_Kinematic_Rigid_Bodies(0);
    SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Initialize_Bodies();
}
//#####################################################################
// Function Solids_Example_Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY<VECTOR_2D<T> >& V,const T time)
{
    Zero_Out_Enslaved_Velocity_Nodes(V,time,fragment_id);
}
//#####################################################################
// Function Solids_Example_Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_2D<T> >& V,const T time)
{
    if(example_number==1){
        if(fragment_id==FRAGMENT_ID(1)){V(1)=VECTOR_2D<T>(0,0);V(solids_parameters.deformable_body_parameters.list.deformable_objects(fragment_id)->segmented_curve->particles.array_collection->Size())=VECTOR_2D<T>(0,0);}}
    else if(example_number==3){
        if(fragment_id==FRAGMENT_ID(1)){V(1)=VECTOR_2D<T>(0,0);}}
    else if(example_number==9){
        if(fragment_id==FRAGMENT_ID(1)){V(1)=VECTOR_2D<T>(0,0);V(solids_parameters.deformable_body_parameters.list.deformable_objects(fragment_id)->segmented_curve->particles.array_collection->Size())=VECTOR_2D<T>(0,0);}}
}
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
void Update_Time_Varying_Material_Properties(const T time)
{
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    INTERPOLATION_CURVE<T,VECTOR_2D<T> > position;
    INTERPOLATION_CURVE<T,T> angle;
    if(example_number==2 || example_number==4){
        position.Add_Control_Point(0,VECTOR_2D<T>(.5,.1));
        position.Add_Control_Point(5,VECTOR_2D<T>(.5,.6));
    }
    else if(example_number==5){
        position.Add_Control_Point(0,VECTOR_2D<T>(.3,.1));
        position.Add_Control_Point(5,VECTOR_2D<T>(.3,.6));
        angle.Add_Control_Point(7,0);
        angle.Add_Control_Point(9,-0.5*pi);
    }
    frame.t=position.Value(time);
    frame.r=COMPLEX<T>::Unit_Polar(angle.Value(time));
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time)
{
    SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Update_Fluid_Parameters(dt,time);
    if(example_number==7){use_source=time>1;}
    if(!use_source_override)use_source=false;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    // Not so good to set up a heaviside function here because then the interface will
    // be exactly between the two nodes which can lead to roundoff issues when setting dirichlet cells, etc.
    GRID<TV>& grid=fluids_parameters.grid;
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++)
        fluids_parameters.particle_levelset_evolution.phi(i,j)=grid.y(j)-grid.ymin-initial_water_level;

    if(example_number==11){
        RIGID_BODY<TV>* rigid_body=solids_parameters.rigid_body_parameters.list(1);
        SEGMENTED_CURVE_2D<T>* segmented_curve=rigid_body->segmented_curve;
        for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++){
            T object_phi=segmented_curve->Calculate_Signed_Distance(rigid_body->Object_Space_Point(grid.X(i,j)));
            T& phi=fluids_parameters.particle_levelset_evolution.phi(i,j);
            if(fabs(object_phi)<fabs(phi)){phi=-object_phi;}}}
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
void Adjust_Phi_With_Sources(const T time)
{
    if(!use_source) return;
    SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Adjust_Phi_With_Source(source,world_to_source);
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,2> >*& cell_centered_mask,const T time)
{
    if(!use_source) return;
    SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Get_Source_Reseed_Mask(source,world_to_source,cell_centered_mask,true);
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time)
{
    if(!use_source) return;
    SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Get_Source_Velocities(source,world_to_source,source_velocity);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time)
{
    if(fluids_parameters.use_solid_fluid_coupling){
//        if(!really_use_effective_velocities){std::cout << "not doing effective velocities" << std::endl;use_object_pseudo_velocity_for_boundary_conditions=false;}
        SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Set_Fluid_Boundary_Conditions();
    }
}
//#####################################################################
// Function Adjust_Particle_For_Objects
//#####################################################################
bool Adjust_Particle_For_Objects(PARTICLE_LEVELSET_PARTICLES<T,VECTOR_2D<T> >& particles,const int index,VECTOR_2D<T>& V,const typename PARTICLE_LEVELSET<T,VECTOR_2D<T> >::PARTICLE_TYPE particle_type,const T dt,const T time)
{
    if(fluids_parameters.use_solid_fluid_coupling)
        return fluids_parameters.Adjust_Particle_For_Objects(fluids_parameters.collision_bodies_affecting_fluid,fluids_parameters.particle_levelset_evolution.particle_levelset.Particle_Collision_Distance(particles.quantized_collision_distance(index)),particles,index,V,particle_type,dt,time);
    else
        return true;
}
//#####################################################################
};    
}
#endif  



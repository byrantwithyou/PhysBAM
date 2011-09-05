//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class THIN_SHELLS_TEST_SMOKE_AND_FIRE
//##################################################################### 
#ifndef __THIN_SHELLS_TEST_SMOKE_AND_FIRE__
#define __THIN_SHELLS_TEST_SMOKE_AND_FIRE__

#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_2D.h>

namespace PhysBAM{

template<class T,class RW=T>
class THIN_SHELLS_TEST_SMOKE_AND_FIRE:public SOLIDS_FLUIDS_EXAMPLE_2D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::solids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::write_output_files;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::verbose_dt;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::data_directory;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::abort_when_dt_below;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::write_time;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::write_frame_title;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Use_Thin_Shells_Fluid_Coupling_Defaults;

    int example_number;
    T rho;

    bool use_source;
    BOX_2D<T> source_domain;
    MATRIX<T,3> world_to_source;
    VECTOR_2D<T> source_velocity;

    GRID<TV> body_grid;
    PARAMETER_LIST parameter_list;
    bool use_source_override;

    // EXAMPLES NUMBERED >= 100 ARE FIRE EXAMPLES!
    THIN_SHELLS_TEST_SMOKE_AND_FIRE(const int example_number_input)
        :SOLIDS_FLUIDS_EXAMPLE_2D<RW>((example_number_input>=100)?fluids_parameters.FIRE:fluids_parameters.SMOKE),example_number(example_number_input),
        use_source(false),world_to_source(MATRIX<T,3>::Identity_Matrix())
    {
        example_number=example_number_input;
        std::cout << "Example number " << example_number << std::endl;
        {std::string filename=STRING_UTILITIES::string_sprintf("Thin_Shells_Test_Smoke_And_Fire/example_%d.param",example_number);
        if(FILE_UTILITIES::File_Exists(filename)){std::cout << "Reading parameter file '" << filename << "'" << std::endl;parameter_list.Read(filename);}}

        // Two way coupling
        fluids_parameters.use_solid_fluid_coupling=true;
        Set_Parameter(fluids_parameters.use_solid_fluid_coupling,"fluids_parameters.use_solid_fluid_coupling");
        if(fluids_parameters.use_solid_fluid_coupling) Use_Thin_Shells_Fluid_Coupling_Defaults();

        // Common parameters
        output_directory="Thin_Shells_Test_Smoke_And_Fire/output";
        first_frame=0;last_frame=1000;frame_rate=60;
        restart=false;restart_frame=0;

        // Fluids parameters 
        if(example_number==103) fluids_parameters.grid.Initialize(226,301,0,1.5,0,2);
        else if(example_number>=100) fluids_parameters.grid.Initialize(301,301,0,2,0,2);
        else if(example_number==6) fluids_parameters.grid.Initialize(129,129,0,4,0,4);
        else{
            //fluids_parameters.grid.Initialize(97,65,0,1.5,0,1);
            fluids_parameters.grid.Initialize(65,65,0,1,0,1);
        }

        fluids_parameters.gravity=0;fluids_parameters.cfl=0.9;
        fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][2]=fluids_parameters.domain_walls[2][2]=false;fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.kolmogorov=(T)0;
        fluids_parameters.implicit_viscosity=false;
        fluids_parameters.incompressible_iterations=100;

        if(fluids_parameters.fire){
            fluids_parameters.phi_boundary->Set_Fixed_Boundary(false);
            fluids_parameters.normal_flame_speed=(T).03;
            fluids_parameters.buoyancy_constant=(T)0;
            fluids_parameters.use_vorticity_confinement_fuel=fluids_parameters.use_vorticity_confinement=true;
            fluids_parameters.confinement_parameter_fuel=fluids_parameters.confinement_parameter=(T).6;}
        else{ // smoke defaults
            fluids_parameters.use_vorticity_confinement=true;fluids_parameters.confinement_parameter=(T).3;}

        // density
        if(fluids_parameters.fire){
            fluids_parameters.density_container.Set_Ambient_Density(0);
            fluids_parameters.density_fuel=(T)1;
            fluids_parameters.density=(T).1;}
        else{ // smoke defaults
            rho=1;}

        // temperature
        if(fluids_parameters.fire){
            fluids_parameters.temperature_container.Set_Ambient_Temperature(T(283.15));
            fluids_parameters.temperature_container.Set_Cooling_Constant((T)4000);
            fluids_parameters.temperature_products=3000;fluids_parameters.temperature_fuel=298;}
        else{ // smoke defaults
            fluids_parameters.temperature_container.Set_Cooling_Constant(0);fluids_parameters.temperature_products=3000;}

        // Solids parameters
        solids_parameters.gravity=(T)9.8;
        solids_parameters.cfl=(T)0.9;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies=true;
        solids_parameters.deformable_body_parameters.write=true;
        verbose_dt=true;

        if(example_number==0){
            use_source=true;source_domain=BOX_2D<T>((T)0,(T).1,(T).45,(T).55);source_velocity=VECTOR_2D<T>((T)1,0);}
        else if(example_number==1){
            use_source=true;source_domain=BOX_2D<T>((T)0,(T).1,(T).45,(T).55);source_velocity=VECTOR_2D<T>((T)10,0);}
        else if(example_number==2){
            use_source=true;source_domain=BOX_2D<T>((T)0,(T).1,(T).45,(T).55);source_velocity=VECTOR_2D<T>((T)10,0);}
        else if(example_number==3){
            use_source=true;source_domain=BOX_2D<T>((T)0,(T).1,(T).45,(T).55);source_velocity=VECTOR_2D<T>((T)20,0);}
        else if(example_number==4){
            PHYSBAM_FATAL_ERROR("Preroll is gone; this example is likely broken.");
            use_source=false;}
        else if(example_number==5){
            use_source=false;
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;}
        else if(example_number==6){
            fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[2][1]=true;fluids_parameters.domain_walls[1][2]=fluids_parameters.domain_walls[2][2]=false;
            PHYSBAM_FATAL_ERROR("Preroll is gone; this example is likely broken.");
            use_source=true;source_domain=BOX_2D<T>(0,.5,1.2,1.8);source_velocity=VECTOR_2D<T>(10,0);}
        else if(example_number==100){
            use_source=true;source_domain=BOX_2D<T>((T).85,(T)1.15,(T)-.5,(T).2);source_velocity=VECTOR_2D<T>((T)0,.5);}
        else if(example_number==101){
            use_source=true;source_domain=BOX_2D<T>((T).85,(T)1.15,(T)-.5,(T).2);source_velocity=VECTOR_2D<T>((T)0,.5);}
        else if(example_number==102){
            PHYSBAM_FATAL_ERROR("Preroll is gone; this example is likely broken.");
            fluids_parameters.buoyancy_constant=(T).01;
            use_source=true;source_domain=BOX_2D<T>((T).85,(T)1.15,(T)-.5,(T).2);source_velocity=VECTOR_2D<T>((T)0,.5);}
        else if(example_number==103){
            PHYSBAM_FATAL_ERROR("Preroll is gone; this example is likely broken.");
            fluids_parameters.buoyancy_constant=(T).001;
            use_source=true;source_domain=BOX_2D<T>((T)0,(T).2,(T).35,(T).65);source_velocity=VECTOR_2D<T>((T).5,0);}

        // Debugging
        fluids_parameters.write_debug_data=true;
        fluids_parameters.write_ghost_values=true;
        abort_when_dt_below=1e-7;
        write_time=true;
        write_frame_title=true;
        if(!fluids_parameters.fire){ // need these on for fire
            fluids_parameters.use_density=true;
            fluids_parameters.use_temperature=false;}

        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Parameters_From_Parameter_List(*this,parameter_list);
        use_source_override=true;Set_Parameter(use_source_override,"use_source");
        Set_Parameter(source_velocity,"source_velocity");
    }

    ~THIN_SHELLS_TEST_SMOKE_AND_FIRE()
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
    bool add_bodies=true;Set_Parameter(add_bodies,"add_bodies");
    if(!add_bodies){SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Initialize_Bodies();return;}

    int id;DEFORMABLE_OBJECT_2D<T>* deformable_object=0;

    if(example_number==0){
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(solids_parameters.deformable_body_parameters.list,30,VECTOR_2D<T>((T).533,(T).8),VECTOR_2D<T>((T).533,(T).2));
        deformable_object=solids_parameters.deformable_body_parameters.list.deformable_objects(id);
        T density=1;T pressure_force_scale=1;THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Deformable_Object_Parameters_2D(id,density,pressure_force_scale,parameter_list);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Density(*deformable_object->segmented_curve,density);
        deformable_object->Add_Body_Forces(*deformable_object->segmented_curve,solids_parameters.gravity,solids_parameters.gravity_direction.Vector_2D());
        deformable_object->Add_Edge_Springs(deformable_object->segmented_curve->segment_mesh,3000,(T).9);
        Add_To_Fluid_Simulation(*deformable_object,true,true,pressure_force_scale);}
    else if(example_number==1){
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Circle_Deformable_Object(solids_parameters.deformable_body_parameters.list,30,VECTOR_2D<T>((T).5,(T).5),(T).25);
        deformable_object=solids_parameters.deformable_body_parameters.list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Mass(*deformable_object->segmented_curve,10,true);
        deformable_object->Add_Body_Forces(*deformable_object->segmented_curve,solids_parameters.gravity,solids_parameters.gravity_direction.Vector_2D());
        deformable_object->Add_Edge_Springs(deformable_object->segmented_curve->segment_mesh,3000,(T).9);
        Add_To_Fluid_Simulation(*deformable_object,true,true,100);}
    else if(example_number==2){
        body_grid=GRID_2D<T>(2,2,.2,.5,.4,.7);
//        body_grid=GRID_2D<T>(4,4,.2,.5,.4,.7);
//        body_grid=GRID_2D<T>(2,2,.4,.44,.4,.42);
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Grid_Deformable_Object(solids_parameters.deformable_body_parameters.list,body_grid,10);
        deformable_object=solids_parameters.deformable_body_parameters.list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Mass(*deformable_object->segmented_curve,10,true);
        deformable_object->Add_Body_Forces(*deformable_object->segmented_curve,solids_parameters.gravity,solids_parameters.gravity_direction.Vector_2D());
        deformable_object->Add_Edge_Springs(deformable_object->segmented_curve->segment_mesh,3000,(T).9);
        Add_To_Fluid_Simulation(*deformable_object,true,true,1);}
    else if(example_number==3){
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(solids_parameters.deformable_body_parameters.list,30,VECTOR_2D<T>((T).533,(T).8),VECTOR_2D<T>((T).533,(T).2));
        deformable_object=solids_parameters.deformable_body_parameters.list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Mass(*deformable_object->segmented_curve,1,true);
        deformable_object->Add_Body_Forces(*deformable_object->segmented_curve,solids_parameters.gravity,solids_parameters.gravity_direction.Vector_2D());
        deformable_object->Add_Edge_Springs(deformable_object->segmented_curve->segment_mesh,30,(T).9);
        Add_To_Fluid_Simulation(*deformable_object,true,true,1);}
    else if(example_number==4){
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(solids_parameters.deformable_body_parameters.list,30,VECTOR_2D<T>((T).433,(T).833),VECTOR_2D<T>((T).633,(T).863));
        deformable_object=solids_parameters.deformable_body_parameters.list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Mass(*deformable_object->segmented_curve,1,true);
        deformable_object->segmented_curve->Set_Density((T).25);deformable_object->segmented_curve->Set_Mass_Of_Particles();
        deformable_object->Add_Body_Forces(*deformable_object->segmented_curve,solids_parameters.gravity,solids_parameters.gravity_direction.Vector_2D());
        deformable_object->Add_Edge_Springs(deformable_object->segmented_curve->segment_mesh,30,(T).9);
        Add_To_Fluid_Simulation(*deformable_object,true,true,1);

        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies_2D/ground",(T).005);
        solids_parameters.rigid_body_parameters.list(id)->position=VECTOR_2D<T>(0.5,0);
        solids_parameters.rigid_body_parameters.list(id)->is_static=true;

        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies_2D/ground",(T).005);
        solids_parameters.rigid_body_parameters.list(id)->position=VECTOR_2D<T>(1,0.5);
        solids_parameters.rigid_body_parameters.list(id)->orientation=pi/2;
        solids_parameters.rigid_body_parameters.list(id)->is_static=true;

        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies_2D/ground",(T).005);
        solids_parameters.rigid_body_parameters.list(id)->position=VECTOR_2D<T>(0,0.5);
        solids_parameters.rigid_body_parameters.list(id)->orientation=-pi/2;
        solids_parameters.rigid_body_parameters.list(id)->is_static=true;

        solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);}
    else if(example_number==5){
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies_2D/ground_subdivided",(T).001);
        RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list(id);
        rigid_body.position=VECTOR_2D<T>(0.533,.833);
        rigid_body.orientation=.1;
        rigid_body.segmented_curve->Set_Density(1);rigid_body.Set_Mass(rigid_body.segmented_curve->Thin_Shell_Mass());
        rigid_body.Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction.Vector_2D(),solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
        Add_To_Fluid_Simulation(rigid_body,true,true,10);}
    if(example_number==6){
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(solids_parameters.deformable_body_parameters.list,30,VECTOR_2D<T>((T)2.003,(T)3.8),VECTOR_2D<T>((T)2.003,(T)1));
        deformable_object=solids_parameters.deformable_body_parameters.list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Density(*deformable_object->segmented_curve,1);
        deformable_object->Add_Body_Forces(*deformable_object->segmented_curve,solids_parameters.gravity,solids_parameters.gravity_direction.Vector_2D());
        deformable_object->Add_Edge_Springs(deformable_object->segmented_curve->segment_mesh,3000,(T).9);
        Add_To_Fluid_Simulation(*deformable_object,true,true,1);}
    else if(example_number==100){
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies_2D/ground",(T).005);
        RIGID_BODY<TV>* rigid_body=solids_parameters.rigid_body_parameters.list(id);
        rigid_body->position=VECTOR_2D<T>(0.9,0.7123);
        rigid_body->is_static=true;
        Add_To_Fluid_Simulation(*rigid_body,true,false);}
    else if(example_number==101){
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(solids_parameters.deformable_body_parameters.list,30,VECTOR_2D<T>((T)1.003,(T)1.8),VECTOR_2D<T>((T)1.993,(T)1.7));
        deformable_object=solids_parameters.deformable_body_parameters.list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Mass(*deformable_object->segmented_curve,1,true);
        deformable_object->Add_Body_Forces(*deformable_object->segmented_curve,solids_parameters.gravity,solids_parameters.gravity_direction.Vector_2D());
        deformable_object->Add_Edge_Springs(deformable_object->segmented_curve->segment_mesh,3000,(T).9);
        Add_To_Fluid_Simulation(*deformable_object,true,false);}
    else if(example_number==102){
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(solids_parameters.deformable_body_parameters.list,20,VECTOR_2D<T>((T)0.503,(T)0.603),VECTOR_2D<T>((T)1.503,(T)0.603));
        T density=1;T pressure_force_scale=1;THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Deformable_Object_Parameters_2D(id,density,pressure_force_scale,parameter_list);
        bool affected_by_fluid=(pressure_force_scale!=0);
        deformable_object=solids_parameters.deformable_body_parameters.list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Density(*deformable_object->segmented_curve,density);
        deformable_object->Add_Body_Forces(*deformable_object->segmented_curve,solids_parameters.gravity,solids_parameters.gravity_direction.Vector_2D());
        deformable_object->Add_Edge_Springs(deformable_object->segmented_curve->segment_mesh,3000,(T)1.1);
        Add_To_Fluid_Simulation(*deformable_object,true,affected_by_fluid,pressure_force_scale);}
    else if(example_number==103){
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(solids_parameters.deformable_body_parameters.list,20,VECTOR_2D<T>((T)0.503,(T)1.303),VECTOR_2D<T>((T).503,(T).303));
        T density=1;T pressure_force_scale=1;THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Deformable_Object_Parameters_2D(id,density,pressure_force_scale,parameter_list);
        deformable_object=solids_parameters.deformable_body_parameters.list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Density(*deformable_object->segmented_curve,density);
        deformable_object->Add_Body_Forces(*deformable_object->segmented_curve,solids_parameters.gravity,solids_parameters.gravity_direction.Vector_2D());
        deformable_object->Add_Edge_Springs(deformable_object->segmented_curve->segment_mesh,3000,(T)1.1);
        Add_To_Fluid_Simulation(*deformable_object,true,true,pressure_force_scale);}

    SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Initialize_Bodies();
}
//#####################################################################
// Function Solids_Example_Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY<VECTOR_2D<T> >& V,const T time)
{
    Zero_Out_Enslaved_Velocity_Nodes(V,time,id_number);
}
//#####################################################################
// Function Solids_Example_Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_2D<T> >& V,const T time)
{
    if(example_number==0) V(1)=VECTOR_2D<T>(0,0);
    else if(example_number==1) for(int i=-5;i<=5;i++) V(((i+V.m)%V.m)+1)=VECTOR_2D<T>(0,0);
    else if(example_number==2) for(int i=1;i<=body_grid.m;i++) V(i*body_grid.n)=VECTOR_2D<T>(0,0);
    else if(example_number==3){V(1)=V(V.m)=VECTOR_2D<T>(0,0);}
    else if(example_number==6){V(1)=VECTOR_2D<T>(0,0);}
    else if(example_number==101){V(1)=VECTOR_2D<T>();if(time<1)V(V.m)=VECTOR_2D<T>();}
    else if(example_number==102){V(1)=V(V.m)=VECTOR_2D<T>();}
    else if(example_number==103){V(1)=VECTOR_2D<T>();}
}
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
void Update_Time_Varying_Material_Properties(const T time)
{
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time)
{
    SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Update_Fluid_Parameters(dt,time);
    if(!use_source_override)use_source=false;
}
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Adjust_Density_And_Temperature_With_Sources(const T time)
{
    if(fluids_parameters.fire || !use_source) return;
    SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Adjust_Density_And_Temperature_With_Sources(source_domain,world_to_source,rho,fluids_parameters.temperature_products,time);
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time)
{
    if(!use_source) return;
    SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Get_Source_Velocities(source_domain,world_to_source,source_velocity,time);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time)
{
    SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Set_Fluid_Boundary_Conditions();
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    assert(fluids_parameters.fire);
    if(!use_source) return;
    for(int i=1;i<=fluids_parameters.grid.m;i++)for(int j=1;j<=fluids_parameters.grid.n;j++)
        if(source_domain.Lazy_Inside(fluids_parameters.grid.X(i,j))) fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j)=-fluids_parameters.grid.dx;
        else fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j)=fluids_parameters.grid.dx;
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
void Adjust_Phi_With_Sources(const T time)
{
    assert(fluids_parameters.fire);
    if(!use_source) return;
    for(int i=1;i<=fluids_parameters.grid.m;i++)for(int j=1;j<=fluids_parameters.grid.n;j++)
        if(source_domain.Lazy_Inside(fluids_parameters.grid.X(i,j))&&fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j)>0)
            fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j)=-fluids_parameters.grid.dx;
}
//#####################################################################
};    
}
#endif

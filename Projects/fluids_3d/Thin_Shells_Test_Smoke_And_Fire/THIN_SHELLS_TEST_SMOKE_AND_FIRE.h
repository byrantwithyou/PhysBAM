//#####################################################################
// Copyright 2004, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class THIN_SHELLS_TEST_SMOKE_AND_FIRE
//##################################################################### 
#ifndef __THIN_SHELLS_TEST_SMOKE_AND_FIRE__
#define __THIN_SHELLS_TEST_SMOKE_AND_FIRE__

#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>

namespace PhysBAM{

template<class T,class RW=T>
class THIN_SHELLS_TEST_SMOKE_AND_FIRE:public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::solids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_output_files;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::verbose_dt;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::abort_when_dt_below;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::data_directory;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_time;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_frame_title;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Use_Thin_Shells_Fluid_Coupling_Defaults;

    enum SOURCE_TYPE {BOX_SOURCE,CYLINDER_SOURCE};
    int example_number;
    T rho;
    bool use_source;
    SOURCE_TYPE source_type;
    BOX_3D<T> box_source;
    CYLINDER<T> cylinder_source;
    MATRIX<T,4> world_to_source;
    VECTOR<T,3> source_velocity;
    ARRAY<ARRAY<int> > deformable_object_enslaved_nodes;
    PARAMETER_LIST parameter_list;

    // EXAMPLES NUMBERED >= 100 ARE FIRE EXAMPLES!
    THIN_SHELLS_TEST_SMOKE_AND_FIRE(const int example_number_input)
        :SOLIDS_FLUIDS_EXAMPLE_3D<RW>((example_number_input>=100)?fluids_parameters.FIRE:fluids_parameters.SMOKE),example_number(example_number_input),use_source(false),source_type(BOX_SOURCE),
        world_to_source(MATRIX<T,4>::Identity_Matrix())
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
        if(example_number==4) fluids_parameters.grid.Initialize(211,141,141,-1,5,0,4,0,4);
        else if(example_number==100) fluids_parameters.grid.Initialize(101,151,101,T(0),T(1),0,(T)1.5,T(0),T(1));
        else if(example_number==101) fluids_parameters.grid.Initialize(101,151,101,T(0),T(1),0,(T)1.5,T(0),T(1));
        else if(example_number==102) fluids_parameters.grid.Initialize(81,121,81,T(0),T(1),0,(T)1.5,T(0),T(1));
        else if(example_number==103) fluids_parameters.grid.Initialize(81,121,81,T(0),T(1),0,(T)1.5,T(0),T(1));
        else if(example_number==104) fluids_parameters.grid.Initialize(161,81,81,-1,7,0,(T)4,T(0),T(4));
        else if(example_number==105) fluids_parameters.grid.Initialize(161,81,81,-1,7,0,(T)4,T(0),T(4));
        else fluids_parameters.grid.Initialize(33,33,33,0,1,0,1,0,1);

        fluids_parameters.gravity=0;fluids_parameters.cfl=0.9;
        fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][2]=fluids_parameters.domain_walls[2][2]=fluids_parameters.domain_walls[3][1]=fluids_parameters.domain_walls[3][2]=false;
        fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.kolmogorov=(T)0;
        fluids_parameters.implicit_viscosity=false;
        fluids_parameters.incompressible_iterations=100;

        if(fluids_parameters.fire){
            fluids_parameters.normal_flame_speed=T(.1);
            fluids_parameters.buoyancy_constant=T(.001);
            fluids_parameters.use_vorticity_confinement_fuel=fluids_parameters.use_vorticity_confinement=true;
            fluids_parameters.confinement_parameter_fuel=fluids_parameters.confinement_parameter=(T).6;}
        else{ // smoke defaults
            fluids_parameters.use_vorticity_confinement=true;fluids_parameters.confinement_parameter=(T).3;}

        // density
        if(fluids_parameters.fire){
            fluids_parameters.density_container.Set_Ambient_Density(0);
            fluids_parameters.density_fuel=1;
            fluids_parameters.density=T(0.1);}
        else{ // smoke defaults
            rho=1;}

        // temperature
        if(fluids_parameters.fire){
            fluids_parameters.temperature_container.Set_Ambient_Temperature((T)283.15);
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
        solids_parameters.perform_self_collision=true;
        verbose_dt=false;

        if(example_number==1){
            use_source=true;source_type=BOX_SOURCE;box_source=BOX_3D<T>((T)0,(T).1,(T).50,(T).60,(T).45,(T).55);source_velocity=VECTOR<T,3>((T)3,0,0);}
        else if(example_number==2){
            use_source=true;source_type=BOX_SOURCE;box_source=BOX_3D<T>((T)0,(T).1,(T).50,(T).60,(T).45,(T).55);source_velocity=VECTOR<T,3>((T)3,0,0);}
        else if(example_number==3){
            PHYSBAM_FATAL_ERROR("Preroll is gone; this example is likely broken.");
            use_source=true;source_type=CYLINDER_SOURCE;cylinder_source=CYLINDER<T>(VECTOR<T,3>(0,.55,.5),VECTOR<T,3>(.1,.55,.5),.05);world_to_source=MATRIX<T,4>::Identity_Matrix();source_velocity=VECTOR<T,3>(20,0,0);
            solids_parameters.perform_self_collision=false;}
        else if(example_number==4){
            fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[2][1]=true;
            fluids_parameters.domain_walls[1][2]=fluids_parameters.domain_walls[2][2]=fluids_parameters.domain_walls[3][1]=fluids_parameters.domain_walls[3][2]=false;
            PHYSBAM_FATAL_ERROR("Preroll is gone; this example is likely broken.");
            use_source=true;source_type=CYLINDER_SOURCE;cylinder_source=CYLINDER<T>(VECTOR<T,3>(0,1.5,2),VECTOR<T,3>(.5,1.5,2),.3);world_to_source=MATRIX<T,4>::Identity_Matrix();source_velocity=VECTOR<T,3>(10,0,0);
            solids_parameters.cfl=(T)5.9;
            solids_parameters.perform_self_collision=true;}
        else if(example_number==100){
            use_source=true;source_type=CYLINDER_SOURCE;cylinder_source=CYLINDER<T>(VECTOR<T,3>((T).5,(T)-.25,(T).5),VECTOR<T,3>((T).5,(T).25,(T).5),(T).15);source_velocity=VECTOR<T,3>(0,1,0);}
        else if(example_number==101){
            PHYSBAM_FATAL_ERROR("Preroll is gone; this example is likely broken.");
            solids_parameters.cfl=(T)3;
            use_source=true;source_type=CYLINDER_SOURCE;cylinder_source=CYLINDER<T>(VECTOR<T,3>((T).5,(T)-.25,(T).5),VECTOR<T,3>((T).5,(T).25,(T).5),(T).15);source_velocity=VECTOR<T,3>(0,1,0);}
        else if(example_number==102){
            PHYSBAM_FATAL_ERROR("Preroll is gone; this example is likely broken.");
            fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][2]=fluids_parameters.domain_walls[3][1]=fluids_parameters.domain_walls[3][2]=fluids_parameters.domain_walls[2][1]=true;
            fluids_parameters.domain_walls[2][2]=false;
            solids_parameters.cfl=(T)3;
            use_source=true;source_type=CYLINDER_SOURCE;cylinder_source=CYLINDER<T>(VECTOR<T,3>((T).5,(T)-.25,(T).5),VECTOR<T,3>((T).5,(T).25,(T).5),(T).15);source_velocity=VECTOR<T,3>(0,1,0);}
        else if(example_number==103){
            PHYSBAM_FATAL_ERROR("Preroll is gone; this example is likely broken.");
            solids_parameters.cfl=(T)3;
            use_source=true;source_type=CYLINDER_SOURCE;cylinder_source=CYLINDER<T>(VECTOR<T,3>((T)0,(T).5,(T).5),VECTOR<T,3>((T).2,(T).5,(T).5),(T).15);source_velocity=VECTOR<T,3>(1,0,0);}
        else if(example_number==104){
            fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[2][1]=true;
            fluids_parameters.domain_walls[1][2]=fluids_parameters.domain_walls[2][2]=fluids_parameters.domain_walls[3][1]=fluids_parameters.domain_walls[3][2]=false;
            fluids_parameters.buoyancy_constant=0.15;
            fluids_parameters.confinement_parameter=3; // based on epsilon=60 and h=.05 (in fire paper)
            fluids_parameters.confinement_parameter_fuel=.8; // based on epsilon=16
            fluids_parameters.density_fuel=1;
            fluids_parameters.density=0.01;
            use_source=true;source_type=CYLINDER_SOURCE;cylinder_source=CYLINDER<T>(VECTOR<T,3>((T)0,(T)2,(T)2),VECTOR<T,3>((T).25,(T)2,(T)2),(T).2);source_velocity=VECTOR<T,3>(30,0,0);}
        else if(example_number==105){
            fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[2][1]=true;
            fluids_parameters.domain_walls[1][2]=fluids_parameters.domain_walls[2][2]=fluids_parameters.domain_walls[3][1]=fluids_parameters.domain_walls[3][2]=false;
            fluids_parameters.buoyancy_constant=0.15;
            fluids_parameters.confinement_parameter=3; // based on epsilon=60 and h=.05 (in fire paper)
            fluids_parameters.confinement_parameter_fuel=.8; // based on epsilon=16
            fluids_parameters.density_fuel=1;
            fluids_parameters.density=0.01;
            PHYSBAM_FATAL_ERROR("Preroll is gone; this example is likely broken.");
            solids_parameters.cfl=(T)2.9;
            T source_y=1;Set_Parameter(source_y,"source_y");
            use_source=true;source_type=CYLINDER_SOURCE;cylinder_source=CYLINDER<T>(VECTOR<T,3>((T)0,(T)source_y,(T)2),VECTOR<T,3>((T).25,(T)source_y,(T)2),(T).2);source_velocity=VECTOR<T,3>(30,0,0);}

        // Debugging
        abort_when_dt_below=1e-7;
        write_time=true;
        write_frame_title=true;
        fluids_parameters.write_debug_data=false;
        fluids_parameters.write_thin_shells_advection_cache=false;
        if(!fluids_parameters.fire){ // need these on for fire
            fluids_parameters.use_density=true;
            fluids_parameters.use_temperature=false;}

        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Parameters_From_Parameter_List(*this,parameter_list);
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
    if(!add_bodies){SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Initialize_Bodies();return;}

    T edge_stiffness_scaling=2;
    T altitude_stiffness_scaling=2;

    DEFORMABLE_OBJECT_LIST_3D<T>& deformable_object_list=solids_parameters.deformable_body_parameters.list;
    int id;DEFORMABLE_OBJECT_3D<T>* deformable_object=0;

    if(example_number==1){
        MATRIX<T,4> transform=MATRIX<T,4>::Translation_Matrix(VECTOR<T,3>(0.5,0.25,0.75))*MATRIX<T,4>::Rotation_Matrix_Y_Axis((T).5*pi);
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(deformable_object_list,deformable_object_enslaved_nodes,GRID<TV>(10,15,0,.5,0,.75),transform);
        deformable_object=deformable_object_list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Mass(*deformable_object->triangulated_surface,1);
        deformable_object->Add_Body_Forces(*deformable_object->triangulated_surface);
        deformable_object->Add_Edge_Springs(deformable_object->triangulated_surface->triangle_mesh,edge_stiffness_scaling/(1+sqrt((T)2)),4);
        deformable_object->Add_Altitude_Springs(deformable_object->triangulated_surface->triangle_mesh,altitude_stiffness_scaling*4/(1+sqrt((T)2)),8);
        deformable_object->Add_Bending_Elements(deformable_object->triangulated_surface->triangle_mesh);
        Add_To_Fluid_Simulation(*deformable_object,true,true,40);}
    else if(example_number==2){
#if 0
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies/box",(T)0.2);
        solids_parameters.rigid_body_parameters.list(id)->position=VECTOR<T,3>(0.6003,0.5003,0.5003);
        solids_parameters.rigid_body_parameters.list(id)->is_static=true;
        Add_To_Fluid_Simulation(*solids_parameters.rigid_body_parameters.list(id),true,false,100);
#else
        MATRIX<T,4> transform=MATRIX<T,4>::Translation_Matrix(VECTOR<T,3>(0.6,0.5,0.5))*MATRIX<T,4>::Scale_Matrix(0.2);
        PLANE<T> halfplane(VECTOR<T,3>(0,-1,0),VECTOR<T,3>(0,0.3,0));
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object_From_File(deformable_object_list,deformable_object_enslaved_nodes,
                                                                                          data_directory+"/Rigid_Bodies/sphere.tri",transform,&halfplane);
        deformable_object=deformable_object_list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Mass(*deformable_object->triangulated_surface,1);
        deformable_object->Add_Body_Forces(*deformable_object->triangulated_surface);
        deformable_object->Add_Edge_Springs(deformable_object->triangulated_surface->triangle_mesh,edge_stiffness_scaling/(1+sqrt((T)2)),4);
        deformable_object->Add_Altitude_Springs(deformable_object->triangulated_surface->triangle_mesh,altitude_stiffness_scaling*4/(1+sqrt((T)2)),8);
//        deformable_object->Add_Bending_Elements(deformable_object->triangulated_surface->triangle_mesh);
        Add_To_Fluid_Simulation(*deformable_object,true,true,4);

        transform=MATRIX<T,4>::Translation_Matrix(VECTOR<T,3>(0.4,0.5,0.5))*MATRIX<T,4>::Scale_Matrix(0.2);
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object_From_File(deformable_object_list,deformable_object_enslaved_nodes,
                                                                                          data_directory+"/Rigid_Bodies/sphere.tri",transform,&halfplane);
        deformable_object=deformable_object_list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Mass(*deformable_object->triangulated_surface,1);
        deformable_object->Add_Body_Forces(*deformable_object->triangulated_surface);
        deformable_object->Add_Edge_Springs(deformable_object->triangulated_surface->triangle_mesh,edge_stiffness_scaling/(1+sqrt((T)2)),4);
        deformable_object->Add_Altitude_Springs(deformable_object->triangulated_surface->triangle_mesh,altitude_stiffness_scaling*4/(1+sqrt((T)2)),8);
//        deformable_object->Add_Bending_Elements(deformable_object->triangulated_surface->triangle_mesh);
        Add_To_Fluid_Simulation(*deformable_object,true,true,4);

        solids_parameters.perform_self_collision=false;
#endif
    }
    else if(example_number==3){
        edge_stiffness_scaling=20;altitude_stiffness_scaling=2;
        MATRIX<T,4> transform=MATRIX<T,4>::Translation_Matrix(VECTOR<T,3>(0.4,0,0.75))*MATRIX<T,4>::Rotation_Matrix_Y_Axis((T).5*pi);
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(deformable_object_list,deformable_object_enslaved_nodes,GRID<TV>(20,30,0,.5,.25,1),transform);
        deformable_object=deformable_object_list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Density(*deformable_object->triangulated_surface,.05);
        deformable_object->Add_Body_Forces(*deformable_object->triangulated_surface);
        deformable_object->Add_Edge_Springs(deformable_object->triangulated_surface->triangle_mesh,edge_stiffness_scaling/(1+sqrt((T)2)),4);
        deformable_object->Add_Altitude_Springs(deformable_object->triangulated_surface->triangle_mesh,altitude_stiffness_scaling*4/(1+sqrt((T)2)),8);
        deformable_object->Add_Bending_Elements(deformable_object->triangulated_surface->triangle_mesh);
        Add_To_Fluid_Simulation(*deformable_object,true,true,1);}
    else if(example_number==4){
        //MATRIX<T,4> transform=MATRIX<T,4>::Translation_Matrix(VECTOR<T,3>(2,4,3))*MATRIX<T,4>::Rotation_Matrix_Y_Axis((T).5*pi);
        //id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(deformable_object_list,deformable_object_enslaved_nodes,GRID<TV>(101,151,0,2,-3,0),transform);
        //MATRIX<T,4> transform=MATRIX<T,4>::Translation_Matrix(VECTOR<T,3>(2,4,3))*MATRIX<T,4>::Rotation_Matrix_Y_Axis((T).5*pi);
        //id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(deformable_object_list,deformable_object_enslaved_nodes,GRID<TV>(51,76,0,2,-3,0),transform);
        MATRIX<T,4> transform=MATRIX<T,4>::Translation_Matrix(VECTOR<T,3>(2,4,3))*MATRIX<T,4>::Rotation_Matrix_Z_Axis((T).3)*MATRIX<T,4>::Rotation_Matrix_Y_Axis((T).5*pi);
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(deformable_object_list,deformable_object_enslaved_nodes,GRID<TV>(11,16,0,2,-3,0),transform);
        edge_stiffness_scaling=20;altitude_stiffness_scaling=2;T density=1;T pressure_force_scale=1;
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Deformable_Object_Parameters_3D(id,edge_stiffness_scaling,altitude_stiffness_scaling,density,pressure_force_scale,parameter_list);

        deformable_object=deformable_object_list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Density(*deformable_object->triangulated_surface,density);
        deformable_object->Add_Body_Forces(*deformable_object->triangulated_surface);
        deformable_object->Add_Edge_Springs(deformable_object->triangulated_surface->triangle_mesh,edge_stiffness_scaling/(1+sqrt((T)2)),4);
        deformable_object->Add_Altitude_Springs(deformable_object->triangulated_surface->triangle_mesh,altitude_stiffness_scaling*4/(1+sqrt((T)2)),8);
        deformable_object->Add_Bending_Elements(deformable_object->triangulated_surface->triangle_mesh);
        Add_To_Fluid_Simulation(*deformable_object,true,true,pressure_force_scale);}
    else if(example_number==100){
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies/ground",(T).003);
        RIGID_BODY<TV>* rigid_body=solids_parameters.rigid_body_parameters.list(id);
        rigid_body->position=VECTOR<T,3>(0.5003,0.9003,0.5003);
        rigid_body->orientation=QUATERNION<T>(.7,VECTOR<T,3>(1,2,3));
        rigid_body->is_static=true;
        Add_To_Fluid_Simulation(*rigid_body,true,false);}
    else if(example_number==101){
        MATRIX<T,4> transform=MATRIX<T,4>::Translation_Matrix(VECTOR<T,3>(.5,1.25,.75))*MATRIX<T,4>::Rotation_Matrix_Z_Axis((T).6)*MATRIX<T,4>::Rotation_Matrix_Y_Axis((T).5*pi);
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(deformable_object_list,deformable_object_enslaved_nodes,GRID<TV>(81,121,0,.5,-.75,0),transform,2);
        edge_stiffness_scaling=20;altitude_stiffness_scaling=2;T density=1;T pressure_force_scale=1;
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Deformable_Object_Parameters_3D(id,edge_stiffness_scaling,altitude_stiffness_scaling,density,pressure_force_scale,parameter_list);
        deformable_object=deformable_object_list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Density(*deformable_object->triangulated_surface,density);
        deformable_object->Add_Body_Forces(*deformable_object->triangulated_surface);
        deformable_object->Add_Edge_Springs(deformable_object->triangulated_surface->triangle_mesh,edge_stiffness_scaling/(1+sqrt((T)2)),4);
        deformable_object->Add_Altitude_Springs(deformable_object->triangulated_surface->triangle_mesh,altitude_stiffness_scaling*4/(1+sqrt((T)2)),8);
        deformable_object->Add_Bending_Elements(deformable_object->triangulated_surface->triangle_mesh);
        Add_To_Fluid_Simulation(*deformable_object,true,false);}
    else if(example_number==102){
        MATRIX<T,4> transform=MATRIX<T,4>::Translation_Matrix(VECTOR<T,3>(.2,.7,.8))*MATRIX<T,4>::Rotation_Matrix_Z_Axis((T).5*pi)*MATRIX<T,4>::Rotation_Matrix_Y_Axis((T).5*pi);
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(deformable_object_list,deformable_object_enslaved_nodes,GRID<TV>(61,61,0,.6,-.6,0),transform,2);
        edge_stiffness_scaling=20;altitude_stiffness_scaling=2;T density=1;T pressure_force_scale=1;
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Deformable_Object_Parameters_3D(id,edge_stiffness_scaling,altitude_stiffness_scaling,density,pressure_force_scale,parameter_list);
        deformable_object=deformable_object_list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Density(*deformable_object->triangulated_surface,density);
        deformable_object->Add_Body_Forces(*deformable_object->triangulated_surface);
        deformable_object->Add_Edge_Springs(deformable_object->triangulated_surface->triangle_mesh,edge_stiffness_scaling/(1+sqrt((T)2)),4);
        deformable_object->Add_Altitude_Springs(deformable_object->triangulated_surface->triangle_mesh,altitude_stiffness_scaling*4/(1+sqrt((T)2)),8);
        deformable_object->Add_Bending_Elements(deformable_object->triangulated_surface->triangle_mesh);
        Add_To_Fluid_Simulation(*deformable_object,true,true,pressure_force_scale);}
    else if(example_number==103){
        MATRIX<T,4> transform=MATRIX<T,4>::Translation_Matrix(VECTOR<T,3>(.5,1.15,.75))*MATRIX<T,4>::Rotation_Matrix_Y_Axis((T).5*pi);
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(deformable_object_list,deformable_object_enslaved_nodes,GRID<TV>(61,91,0,.5,-.75,0),transform,2);
        edge_stiffness_scaling=20;altitude_stiffness_scaling=2;T density=1;T pressure_force_scale=1;
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Deformable_Object_Parameters_3D(id,edge_stiffness_scaling,altitude_stiffness_scaling,density,pressure_force_scale,parameter_list);
        deformable_object=deformable_object_list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Density(*deformable_object->triangulated_surface,density);
        deformable_object->Add_Body_Forces(*deformable_object->triangulated_surface);
        deformable_object->Add_Edge_Springs(deformable_object->triangulated_surface->triangle_mesh,edge_stiffness_scaling/(1+sqrt((T)2)),4);
        deformable_object->Add_Altitude_Springs(deformable_object->triangulated_surface->triangle_mesh,altitude_stiffness_scaling*4/(1+sqrt((T)2)),8);
        deformable_object->Add_Bending_Elements(deformable_object->triangulated_surface->triangle_mesh);
        Add_To_Fluid_Simulation(*deformable_object,true,true,pressure_force_scale);}
    else if(example_number==104){
        MATRIX<T,4> transform=MATRIX<T,4>::Translation_Matrix(VECTOR<T,3>(2,4,3))*MATRIX<T,4>::Rotation_Matrix_Y_Axis((T).5*pi);
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(deformable_object_list,deformable_object_enslaved_nodes,GRID<TV>(3,4,0,2,-3,0),transform);
        deformable_object=deformable_object_list.deformable_objects(id);
        deformable_object->simulate=false;
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Density(*deformable_object->triangulated_surface,1);
        Add_To_Fluid_Simulation(*deformable_object,true,false);}
    else if(example_number==105){
        VECTOR<T,3> cloth_corner(2.5,4,3);Set_Parameter(cloth_corner,"cloth_corner");
        VECTOR_2D<T> cloth_size(2,3);Set_Parameter(cloth_size,"cloth_size");
        VECTOR_2D<int> cloth_resolution(101,151);Set_Parameter(cloth_resolution,"cloth_resolution");
        int constraint_mode=1;Set_Parameter(constraint_mode,"constraint_mode");
        MATRIX<T,4> transform=MATRIX<T,4>::Translation_Matrix(cloth_corner)*MATRIX<T,4>::Rotation_Matrix_Y_Axis((T).5*pi);
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(deformable_object_list,deformable_object_enslaved_nodes,GRID<TV>(cloth_resolution.x,cloth_resolution.y,0,cloth_size.x,-cloth_size.y,0),transform,constraint_mode);
        edge_stiffness_scaling=20;altitude_stiffness_scaling=2;T density=1;T pressure_force_scale=1;
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Deformable_Object_Parameters_3D(id,edge_stiffness_scaling,altitude_stiffness_scaling,density,pressure_force_scale,parameter_list);

        deformable_object=deformable_object_list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Density(*deformable_object->triangulated_surface,density);
        deformable_object->Add_Body_Forces(*deformable_object->triangulated_surface);
        deformable_object->Add_Edge_Springs(deformable_object->triangulated_surface->triangle_mesh,edge_stiffness_scaling/(1+sqrt((T)2)),4);
        deformable_object->Add_Altitude_Springs(deformable_object->triangulated_surface->triangle_mesh,altitude_stiffness_scaling*4/(1+sqrt((T)2)),8);
        deformable_object->Add_Bending_Elements(deformable_object->triangulated_surface->triangle_mesh);
        Add_To_Fluid_Simulation(*deformable_object,true,true,pressure_force_scale);
    
        for(int i=1;i<=deformable_object_enslaved_nodes(id).m;i++){
            std::cout << "SETTING PARTICLES " << deformable_object_enslaved_nodes(id)(i) << " to large mass" << std::endl;
            deformable_object->triangulated_surface->particles.mass(deformable_object_enslaved_nodes(id)(i))=1e30;}
    }

    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Initialize_Bodies();
}
//#####################################################################
// Function Solids_Example_Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY<VECTOR<T,3> >& V,const T time)
{
    Zero_Out_Enslaved_Velocity_Nodes(V,time,id_number);
}
//#####################################################################
// Function Solids_Example_Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR<T,3> >& V,const T time)
{
    assert(id_number<=deformable_object_enslaved_nodes.m);
    if(example_number==101){
        if(time<3){for(int i=1;i<=deformable_object_enslaved_nodes(id_number).m;i++) V(deformable_object_enslaved_nodes(id_number)(i))=VECTOR<T,3>();}
        else{for(int i=1;i<=2;i++) V(deformable_object_enslaved_nodes(id_number)(i))=VECTOR<T,3>();}}
    else{for(int i=1;i<=deformable_object_enslaved_nodes(id_number).m;i++) V(deformable_object_enslaved_nodes(id_number)(i))=VECTOR<T,3>();}
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
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Update_Fluid_Parameters(dt,time);
    if(example_number==4) use_source=time<5;
    if(example_number==104 || example_number==105) use_source=time<5;
}
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Adjust_Density_And_Temperature_With_Sources(const T time)
{
    if(fluids_parameters.fire || !use_source) return;
    if(source_type==BOX_SOURCE) SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Adjust_Density_And_Temperature_With_Sources(box_source,world_to_source,rho,fluids_parameters.temperature_products,time);
    else if(source_type==CYLINDER_SOURCE) SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Adjust_Density_And_Temperature_With_Sources(cylinder_source,world_to_source,rho,fluids_parameters.temperature_products,time);
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time)
{
    if(!use_source) return;
    else if(source_type==BOX_SOURCE) SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Get_Source_Velocities(box_source,world_to_source,source_velocity,time);
    else if(source_type==CYLINDER_SOURCE) SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Get_Source_Velocities(cylinder_source,world_to_source,source_velocity,time);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time)
{
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Set_Fluid_Boundary_Conditions();
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
template<class SOURCE> void Initialize_Phi(const SOURCE& source)
{
    for(int i=1;i<=fluids_parameters.grid.m;i++) for(int j=1;j<=fluids_parameters.grid.n;j++) for(int ij=1;ij<=fluids_parameters.grid.mn;ij++) 
        if(source.Lazy_Inside(world_to_source*fluids_parameters.grid.X(i,j,ij)))
            fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j,ij)=-fluids_parameters.grid.dx;
        else
            fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j,ij)=fluids_parameters.grid.dx;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    assert(fluids_parameters.fire);
    if(!use_source) return;
    else if(source_type==BOX_SOURCE) Initialize_Phi(box_source);
    else if(source_type==CYLINDER_SOURCE) Initialize_Phi(cylinder_source);
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
template<class SOURCE> void Adjust_Phi_With_Sources(const SOURCE& source,const T time)
{
    for(int i=1;i<=fluids_parameters.grid.m;i++) for(int j=1;j<=fluids_parameters.grid.n;j++) for(int ij=1;ij<=fluids_parameters.grid.mn;ij++) 
//        if(source.Lazy_Inside(world_to_source*fluids_parameters.grid.X(i,j,ij)))
        if(source.Lazy_Inside(world_to_source*fluids_parameters.grid.X(i,j,ij))&&fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j,ij)>0)
            fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j,ij)=-fluids_parameters.grid.dx;
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
void Adjust_Phi_With_Sources(const T time)
{
    assert(fluids_parameters.fire);
    if(!use_source) return;
    else if(source_type==BOX_SOURCE) Adjust_Phi_With_Sources(box_source,time);
    else if(source_type==CYLINDER_SOURCE) Adjust_Phi_With_Sources(cylinder_source,time);
}
//#####################################################################
};    
}
#endif

//#####################################################################
// Copyright 2004, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class THIN_SHELLS_TEST_WATER 
//##################################################################### 
#ifndef __THIN_SHELLS_TEST_WATER__
#define __THIN_SHELLS_TEST_WATER__

#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>
namespace PhysBAM{

template<class T,class RW=T>
class THIN_SHELLS_TEST_WATER:public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::solids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::verbose_dt;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_time;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_frame_title;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::data_directory;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::abort_when_dt_below;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Use_Thin_Shells_Fluid_Coupling_Defaults;

    enum SOURCE_TYPE {BOX_SOURCE,CYLINDER_SOURCE};
    int example_number;
    bool use_source;
    SOURCE_TYPE source_type;
    BOX_3D<T> box_source;
    CYLINDER<T> cylinder_source;
    MATRIX<T,4> world_to_source;
    VECTOR<T,3> source_velocity;
    T initial_water_level;
    ARRAY<ARRAY<int> > deformable_object_enslaved_nodes;
    T domain_scale;
    PARAMETER_LIST parameter_list;
    bool use_source_override;
    T source_stop_time;
    T enslaving_stop_time;

    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return false;}

    THIN_SHELLS_TEST_WATER(const int example_number_input)
        :SOLIDS_FLUIDS_EXAMPLE_3D<RW>(fluids_parameters.WATER),source_type(BOX_SOURCE),domain_scale(1)
    {
        parameter_list.Set_Verbose(true);
        example_number=example_number_input;
        std::cout << "Example number " << example_number << std::endl;
        {std::string filename=STRING_UTILITIES::string_sprintf("Thin_Shells_Test_Water/example_%d.param",example_number);
        if(FILE_UTILITIES::File_Exists(filename)){std::cout << "Reading parameter file '" << filename << "'" << std::endl;parameter_list.Read(filename);}}

        // Two way coupling
        fluids_parameters.use_solid_fluid_coupling=true;
        Set_Parameter(fluids_parameters.use_solid_fluid_coupling,"fluids_parameters.use_solid_fluid_coupling");
        if(fluids_parameters.use_solid_fluid_coupling) Use_Thin_Shells_Fluid_Coupling_Defaults();

        // Common parameters
        output_directory="Thin_Shells_Test_Water/output";
        first_frame=0;last_frame=1000;frame_rate=60;
        restart=false;restart_frame=0;

        if(example_number==2)domain_scale=6; // different default
        Set_Parameter(domain_scale,"domain_scale");

        // Fluids parameters
        if(example_number==8||example_number==9){
            fluids_parameters.grid.Initialize(101,101,76,0,domain_scale,0,domain_scale,domain_scale*.125,domain_scale*.875);}
        else if(example_number==2){
            fluids_parameters.grid.Initialize(101,101,81,0,domain_scale,0,domain_scale,.1*domain_scale,.9*domain_scale);}
        else if(example_number==6){
            fluids_parameters.grid.Initialize(101,101,81,0,1,0,1,.1,.9);}
        else if(example_number==12){
            fluids_parameters.grid.Initialize(101,101,81,0,4,0,4,.4,3.6);}
        else if(example_number==14)
            fluids_parameters.grid.Initialize(201,241,161,0,1,0,1.2,0,0.8);
        else{
            fluids_parameters.grid.Initialize(101,101,101,0,1,0,1,0,1);}
        fluids_parameters.number_particles_per_cell=8;
        fluids_parameters.reseeding_frame_rate=10;
        fluids_parameters.incompressible_iterations=20;
        fluids_parameters.cfl=0.9;
        fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][2]=fluids_parameters.domain_walls[2][1]=fluids_parameters.domain_walls[3][1]=fluids_parameters.domain_walls[3][2]=true;fluids_parameters.domain_walls[2][2]=false;
        fluids_parameters.bias_towards_negative_particles=true;
        fluids_parameters.use_removed_positive_particles=false;fluids_parameters.use_removed_negative_particles=true;
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;fluids_parameters.variable_viscosity=false;fluids_parameters.second_order_cut_cell_method=false;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;
        fluids_parameters.write_particles=true;fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=true;
        fluids_parameters.store_particle_ids=true;
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.enforce_divergence_free_extrapolation=false;

        // Solids parameters
        solids_parameters.gravity=9.8;
        solids_parameters.cfl=(T).9;
        solids_parameters.perform_self_collision=true;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies=true;
        solids_parameters.deformable_body_parameters.write=true;
        verbose_dt=true;

        if(example_number==1){
            //initial_water_level=(T)1.021; // some number so the interface won't lie in the center of a cell
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
            initial_water_level=-1; // no initial height
            use_source=true;box_source=BOX_3D<T>((T).45,(T).55,(T).9,1,(T).45,(T).55);world_to_source=MATRIX<T,4>::Identity_Matrix();source_velocity=VECTOR<T,3>(0,-1,0);}
        else if(example_number==2){
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
            initial_water_level=(T).33333*domain_scale; // no initial height
            use_source=false;}
        else if(example_number==3){
            initial_water_level=-1;
            use_source=true;box_source=BOX_3D<T>((T)0,(T).1,(T).45,(T).55,(T).44,(T).55);world_to_source=MATRIX<T,4>::Identity_Matrix();source_velocity=VECTOR<T,3>(3,0,0);}
        else if(example_number==4){
            initial_water_level=-1;
            use_source=true;box_source=BOX_3D<T>((T)0.4,(T).6,(T).9,(T)1,(T).4,(T).6);world_to_source=MATRIX<T,4>::Identity_Matrix();source_velocity=VECTOR<T,3>(0,-2,0);}
        else if(example_number==5){
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
            initial_water_level=(T).40333; // no initial height
            use_source=false;}
        else if(example_number==6){ // water on curtain
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
            PHYSBAM_FATAL_ERROR("Preroll is gone; this example is likely broken.");
            initial_water_level=(T)-1; // no initial height
            use_source=true;source_type=CYLINDER_SOURCE;cylinder_source=CYLINDER<T>(VECTOR<T,3>(0,.65,.5),VECTOR<T,3>(.1,.65,.5),.05);world_to_source=MATRIX<T,4>::Identity_Matrix();source_velocity=VECTOR<T,3>(3,0,0);}
        else if(example_number==7){ // water over cloth funnel
            initial_water_level=(T)-1; // no initial height
            use_source=true;source_type=CYLINDER_SOURCE;cylinder_source=CYLINDER<T>(VECTOR<T,3>(.25,.9,.5),VECTOR<T,3>(.25,1,.5),.075);world_to_source=MATRIX<T,4>::Identity_Matrix();source_velocity=VECTOR<T,3>(0,-.25,0);}
        else if(example_number==8){
            initial_water_level=domain_scale*0.6;
            fluids_parameters.simulate=true;
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
            use_source=false;source_type=CYLINDER_SOURCE;cylinder_source=CYLINDER<T>(domain_scale*VECTOR<T,3>(.55,.9,.5),domain_scale*VECTOR<T,3>(.55,1,.5),domain_scale*.05);world_to_source=MATRIX<T,4>::Identity_Matrix();source_velocity=VECTOR<T,3>(0,-domain_scale*1,0);}
        else if(example_number==9){
            initial_water_level=0.6;
            fluids_parameters.simulate=true;
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
            use_source=false;source_type=CYLINDER_SOURCE;cylinder_source=CYLINDER<T>(VECTOR<T,3>(.55,.9,.5),VECTOR<T,3>(.55,1,.5),.05);world_to_source=MATRIX<T,4>::Identity_Matrix();source_velocity=VECTOR<T,3>(0,-1,0);}
        else if(example_number==10){ // buddha cup profile (start cup in water)
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
            initial_water_level=(T).40333; // no initial height
            use_source=false;}
        else if(example_number==11){ // cup (like buddha example 5)
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
            initial_water_level=(T).40333; // no initial height
            use_source=false;}
        else if(example_number==12){ // water on curtain
            PHYSBAM_FATAL_ERROR("Preroll is gone; this example is likely broken.");
            initial_water_level=(T)-1; // no initial height
            bool longer_cylinder=false;Set_Parameter(longer_cylinder,"longer_cylinder");
            use_source=true;source_type=CYLINDER_SOURCE;cylinder_source=CYLINDER<T>(VECTOR<T,3>(longer_cylinder?-1:0,2.6,2),VECTOR<T,3>(.4,2.6,2),.25);world_to_source=MATRIX<T,4>::Identity_Matrix();source_velocity=VECTOR<T,3>(4,0,0);}
        else if(example_number==14){ // like 11 but with taller domain
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
            initial_water_level=(T).40333; // no initial height
            use_source=false;}

        // Debugging
        fluids_parameters.write_debug_data=true;
        abort_when_dt_below=1e-7;
        write_time=true;
        write_frame_title=true;
        fluids_parameters.simulate=true;
        fluids_parameters.write_ghost_values=false;

        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Parameters_From_Parameter_List(*this,parameter_list);
        use_source_override=true;Set_Parameter(use_source_override,"use_source");
        Set_Parameter(source_velocity,"source_velocity");

        source_stop_time=-1;Set_Parameter(source_stop_time,"source_stop_time");
        enslaving_stop_time=-1;Set_Parameter(enslaving_stop_time,"enslaving_stop_time");
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
// Function Translate_Body_Simplicial_Object
//#####################################################################
void Translate_Body_Simplicial_Object(int id,const VECTOR<T,3>& v)
{
    solids_parameters.rigid_body_parameters.list(id)->simplicial_object->particles.X+=v;
    solids_parameters.rigid_body_parameters.list(id)->simplicial_object->Update_Triangle_List();
    solids_parameters.rigid_body_parameters.list(id)->simplicial_object->Initialize_Hierarchy();
    solids_parameters.rigid_body_parameters.list(id)->simplicial_object->hierarchy->Update_Box_Radii();
    solids_parameters.rigid_body_parameters.list(id)->simplicial_object->Update_Bounding_Box();
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    bool add_bodies=true;Set_Parameter(add_bodies,"add_bodies");
    if(!add_bodies){SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Initialize_Bodies();return;}
    bool simulate_bodies=true;Set_Parameter(simulate_bodies,"simulate_bodies");

    DEFORMABLE_OBJECT_LIST_3D<T>& deformable_object_list=solids_parameters.deformable_body_parameters.list;
    int id;DEFORMABLE_OBJECT_3D<T>* deformable_object=0;

    if(example_number==1){
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies/Thin_Shells/cup",(T)1/100);
        solids_parameters.rigid_body_parameters.list(id)->is_static=false;solids_parameters.rigid_body_parameters.list(id)->is_kinematic=true;
        Add_To_Fluid_Simulation(*solids_parameters.rigid_body_parameters.list(id),true,false,100);}
    else if(example_number==2){
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies/Thin_Shells/cup",(T)domain_scale/100);
        solids_parameters.rigid_body_parameters.list(id)->is_static=false;solids_parameters.rigid_body_parameters.list(id)->is_kinematic=true;
        Add_To_Fluid_Simulation(*solids_parameters.rigid_body_parameters.list(id),true,false,100);}
    else if(example_number==3){
        MATRIX<T,4> transform=MATRIX<T,4>::Translation_Matrix(VECTOR<T,3>(0.5,0.33,0.75))*MATRIX<T,4>::Rotation_Matrix_Y_Axis((T).5*pi);
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(deformable_object_list,deformable_object_enslaved_nodes,GRID<TV>(40,60,0,.5,0,.75),transform);
        deformable_object=deformable_object_list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Mass(*deformable_object->triangulated_surface,1);
        deformable_object->Add_Body_Forces(*deformable_object->triangulated_surface);
        deformable_object->Add_Edge_Springs(deformable_object->triangulated_surface->triangle_mesh,2/(1+sqrt((T)2)),2);
        deformable_object->Add_Altitude_Springs(deformable_object->triangulated_surface->triangle_mesh,2*4/(1+sqrt((T)2)),4);
        deformable_object->Add_Bending_Elements(deformable_object->triangulated_surface->triangle_mesh);
        Add_To_Fluid_Simulation(*deformable_object,true,true,30);}
    else if(example_number==4){
        T edge_stiffness_scaling=2,altitude_stiffness_scaling=2;
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
        Add_To_Fluid_Simulation(*deformable_object,true,false,4);

        transform=MATRIX<T,4>::Translation_Matrix(VECTOR<T,3>(0.4,0.5,0.5))*MATRIX<T,4>::Scale_Matrix(0.2);
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object_From_File(deformable_object_list,deformable_object_enslaved_nodes,
                                                                                       data_directory+"/Rigid_Bodies/sphere.tri",transform,&halfplane);
        deformable_object=deformable_object_list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Mass(*deformable_object->triangulated_surface,1);
        deformable_object->Add_Body_Forces(*deformable_object->triangulated_surface);
        deformable_object->Add_Edge_Springs(deformable_object->triangulated_surface->triangle_mesh,edge_stiffness_scaling/(1+sqrt((T)2)),4);
        deformable_object->Add_Altitude_Springs(deformable_object->triangulated_surface->triangle_mesh,altitude_stiffness_scaling*4/(1+sqrt((T)2)),8);
//        deformable_object->Add_Bending_Elements(deformable_object->triangulated_surface->triangle_mesh);
        Add_To_Fluid_Simulation(*deformable_object,true,false,4);
        solids_parameters.perform_self_collision=false;}
    else if(example_number==5){
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies/Thin_Shells/buddha_tophalf",(T).4);
        Translate_Body_Simplicial_Object(id,VECTOR<T,3>(0,-0.2,0));
        solids_parameters.rigid_body_parameters.list(id)->is_static=false;solids_parameters.rigid_body_parameters.list(id)->is_kinematic=true;
        Add_To_Fluid_Simulation(*solids_parameters.rigid_body_parameters.list(id),true,false,100);}
    else if(example_number==6){
        MATRIX<T,4> transform=MATRIX<T,4>::Translation_Matrix(VECTOR<T,3>(0.4,0,0.75))*MATRIX<T,4>::Rotation_Matrix_Y_Axis((T).5*pi);
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(deformable_object_list,deformable_object_enslaved_nodes,GRID<TV>(101,151,0,.5,.25,1),transform);
        T edge_stiffness_scaling=20,altitude_stiffness_scaling=2;T density=1;T pressure_force_scale=1;
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Deformable_Object_Parameters_3D(id,edge_stiffness_scaling,altitude_stiffness_scaling,density,pressure_force_scale,parameter_list);
        deformable_object=deformable_object_list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Density(*deformable_object->triangulated_surface,density);
        deformable_object->Add_Body_Forces(*deformable_object->triangulated_surface);
        deformable_object->Add_Edge_Springs(deformable_object->triangulated_surface->triangle_mesh,edge_stiffness_scaling/(1+sqrt((T)2)),4);
        deformable_object->Add_Altitude_Springs(deformable_object->triangulated_surface->triangle_mesh,altitude_stiffness_scaling*4/(1+sqrt((T)2)),8);
        deformable_object->Add_Bending_Elements(deformable_object->triangulated_surface->triangle_mesh);
        Add_To_Fluid_Simulation(*deformable_object,true,true,pressure_force_scale);}
    else if(example_number==7){
        T edge_stiffness_scaling=20,altitude_stiffness_scaling=2;
        MATRIX<T,4> transform=MATRIX<T,4>::Translation_Matrix(VECTOR<T,3>(0.4,0.5,0.5))*MATRIX<T,4>::Rotation_Matrix_Z_Axis((T)-.12*pi)*MATRIX<T,4>::Rotation_Matrix_X_Axis((T).5*pi);
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(deformable_object_list,deformable_object_enslaved_nodes,GRID<TV>(30,20,-.375,.375,-.25,.25),transform,2);
        deformable_object=deformable_object_list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Mass(*deformable_object->triangulated_surface,1);
        deformable_object->Add_Body_Forces(*deformable_object->triangulated_surface);
        deformable_object->Add_Edge_Springs(deformable_object->triangulated_surface->triangle_mesh,edge_stiffness_scaling/(1+sqrt((T)2)),4);
        deformable_object->Add_Altitude_Springs(deformable_object->triangulated_surface->triangle_mesh,altitude_stiffness_scaling*4/(1+sqrt((T)2)),8);
        deformable_object->Add_Bending_Elements(deformable_object->triangulated_surface->triangle_mesh);
        Add_To_Fluid_Simulation(*deformable_object,true,true,40);}
    else if(example_number==8){
        T mu=0.9,epsilon=0.2;
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies/Thin_Shells/canoe_high",(T)domain_scale*.15);
        T density=1,pressure_force_scale=3;THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Rigid_Body_Parameters_2D(id,density,pressure_force_scale,parameter_list);
        RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list(id);
        rigid_body.frame.t=domain_scale*VECTOR<T,3>(.5,.69,.5);
        rigid_body.frame.r=QUATERNION<T>(.5*pi,VECTOR<T,3>(0,1,0));
        rigid_body.triangulated_surface->Set_Density(density);rigid_body.Set_Mass(rigid_body.triangulated_surface->Thin_Shell_Mass());
        rigid_body.Set_Coefficient_Of_Friction(mu);
        rigid_body.Set_Coefficient_Of_Restitution(epsilon);
        rigid_body.Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
        Add_To_Fluid_Simulation(rigid_body,true,true,pressure_force_scale);

        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Rigid_Body_Walls(*this,epsilon,mu);

        solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);}
    else if(example_number==9){
        T mu=0.9,epsilon=0.2;
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies/Thin_Shells/boat_hires",(T).20);
        RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list(id);
        rigid_body.frame.t=VECTOR<T,3>(.5,.65,.5);
        rigid_body.frame.r=QUATERNION<T>(-.1,VECTOR<T,3>(0,0,1)) * QUATERNION<T>(.5*pi,VECTOR<T,3>(0,1,0));
        rigid_body.triangulated_surface->Set_Density(1);rigid_body.Set_Mass(rigid_body.triangulated_surface->Thin_Shell_Mass());
        rigid_body.Set_Coefficient_Of_Friction(mu);
        rigid_body.Set_Coefficient_Of_Restitution(epsilon);
        rigid_body.Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
        Add_To_Fluid_Simulation(rigid_body,true,true,1.75);

        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Rigid_Body_Walls(*this,epsilon,mu);

        solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);}
    else if(example_number==10){
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies/Thin_Shells/buddha_tophalf",(T).4);
        Translate_Body_Simplicial_Object(id,VECTOR<T,3>(0,-0.2,0));
        solids_parameters.rigid_body_parameters.list(id)->is_static=false;solids_parameters.rigid_body_parameters.list(id)->is_kinematic=true;
        Add_To_Fluid_Simulation(*solids_parameters.rigid_body_parameters.list(id),true,false,100);}
    else if(example_number==11){
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies/Thin_Shells/cup",(T).01);
        solids_parameters.rigid_body_parameters.list(id)->is_static=false;solids_parameters.rigid_body_parameters.list(id)->is_kinematic=true;
        Add_To_Fluid_Simulation(*solids_parameters.rigid_body_parameters.list(id),true,false,100);}
    else if(example_number==12){ // water curtain
        VECTOR<T,3> cloth_translation(1.2,4,3);Set_Parameter(cloth_translation,"cloth_translation");
        MATRIX<T,4> transform=MATRIX<T,4>::Translation_Matrix(cloth_translation)*MATRIX<T,4>::Rotation_Matrix_Y_Axis((T).5*pi);
        VECTOR_2D<int> cloth_res(41,51);Set_Parameter(cloth_res,"cloth_resolution");
        id=THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Deformable_Object(deformable_object_list,deformable_object_enslaved_nodes,GRID<TV>(cloth_res.x,cloth_res.y,0,2,-2.5,0),transform);
        T edge_stiffness_scaling=20,altitude_stiffness_scaling=2;T density=1;T pressure_force_scale=1;
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Deformable_Object_Parameters_3D(id,edge_stiffness_scaling,altitude_stiffness_scaling,density,pressure_force_scale,parameter_list);
        T bending_stiffness=1e-3,bending_damping=1e-3,bending_max_strain_per_time_step=.1;bool use_constant_mass=false;
        Set_Parameter(bending_stiffness,"deformable_object_1.bending_stiffness");
        Set_Parameter(bending_damping,"deformable_object_1.bending_damping");
        Set_Parameter(bending_max_strain_per_time_step,"deformable_object_1.bending_max_strain_per_time_step");
        Set_Parameter(use_constant_mass,"deformable_object_1.use_constant_mass");
        deformable_object=deformable_object_list.deformable_objects(id);
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Density(*deformable_object->triangulated_surface,density,use_constant_mass);
        if(simulate_bodies){
            deformable_object->Add_Body_Forces(*deformable_object->triangulated_surface);
            deformable_object->Add_Edge_Springs(deformable_object->triangulated_surface->triangle_mesh,edge_stiffness_scaling/(1+sqrt((T)2)),4);
            deformable_object->Add_Altitude_Springs(deformable_object->triangulated_surface->triangle_mesh,altitude_stiffness_scaling*4/(1+sqrt((T)2)),8);
            deformable_object->Add_Bending_Elements(deformable_object->triangulated_surface->triangle_mesh,bending_stiffness,bending_damping,true,bending_max_strain_per_time_step);
            Add_To_Fluid_Simulation(*deformable_object,true,true,pressure_force_scale);}
        else{
            Add_To_Fluid_Simulation(*deformable_object,true,false);}}
    else if(example_number==14){
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies/Thin_Shells/cup",(T).01);
        solids_parameters.rigid_body_parameters.list(id)->is_static=false;solids_parameters.rigid_body_parameters.list(id)->is_kinematic=true;
        Add_To_Fluid_Simulation(*solids_parameters.rigid_body_parameters.list(id),true,false,100);}

    // Move kinematic bodies to their position at time 0
    solids_parameters.rigid_body_parameters.Reset_Kinematic_Rigid_Bodies(0);

    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Initialize_Bodies();
}
//#####################################################################
// Function Solids_Example_Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY<VECTOR<T,3> >& V,const T time)
{  
    Zero_Out_Enslaved_Velocity_Nodes(V,time,fragment_id);
}
//#####################################################################
// Function Solids_Example_Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR<T,3> >& V,const T time)
{
    assert(fragment_id<=deformable_object_enslaved_nodes.m);
    if(enslaving_stop_time!=-1 && time>enslaving_stop_time) return;
    for(int i=0;i<deformable_object_enslaved_nodes(fragment_id).m;i++) V(deformable_object_enslaved_nodes(fragment_id)(i))=VECTOR<T,3>();
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
    T time_scale=1;Set_Parameter(time_scale,"time_scale");
    if(example_number==1){
        INTERPOLATION_CURVE<T,T> angle;angle.Add_Control_Point(0,-pi/2);angle.Add_Control_Point(1,-pi/2);angle.Add_Control_Point(2,pi/4);
        frame.r=QUATERNION<T>(angle.Value(time),VECTOR<T,3>(1,0,0));frame.t=VECTOR<T,3>((T).5,(T).5,(T).5);}
    else if(example_number==2){
        INTERPOLATION_CURVE<T,VECTOR<T,3> > position;position.Add_Control_Point(0,VECTOR<T,3>((T).5*domain_scale,(T).05*domain_scale,(T).5*domain_scale));position.Add_Control_Point(1.5,VECTOR<T,3>((T).5*domain_scale,(T).65*domain_scale,(T).5*domain_scale));
        INTERPOLATION_CURVE<T,T> angle;angle.Add_Control_Point(3,0);angle.Add_Control_Point(5,-(T).55*pi);
        frame.r=QUATERNION<T>(angle.Value(time),VECTOR<T,3>(0,0,1))*QUATERNION<T>(-pi/2,VECTOR<T,3>(1,0,0));
        //frame.r=QUATERNION<T>(-pi/2,VECTOR<T,3>(1,0,0));
        frame.t=position.Value(time);}
    else if(example_number==5){
        INTERPOLATION_CURVE<T,VECTOR<T,3> > position;
        INTERPOLATION_CURVE<T,T> angle;
        position.Add_Control_Point(0,VECTOR<T,3>((T).5,(T).7,(T).5));
        position.Add_Control_Point(0.5,VECTOR<T,3>((T).5,(T).22,(T).5));
        position.Add_Control_Point(5,VECTOR<T,3>((T).5,(T).22,(T).5));
        position.Add_Control_Point(6,VECTOR<T,3>((T).5,(T).7,(T).5));
        angle.Add_Control_Point(3.5,(T)3*pi/4);
        angle.Add_Control_Point(4.5,(T)pi);
        angle.Add_Control_Point(6.5,(T)pi);
        angle.Add_Control_Point(8.0,(T)3*pi/8);
        frame.r=QUATERNION<T>(angle.Value(time),VECTOR<T,3>(0,0,1));
        frame.t=position.Value(time);}
    else if(example_number==10){
        INTERPOLATION_CURVE<T,VECTOR<T,3> > position;
        INTERPOLATION_CURVE<T,T> angle;
        position.Add_Control_Point(0,VECTOR<T,3>((T).5,(T).22,(T).5));
        position.Add_Control_Point(2,VECTOR<T,3>((T).5,(T).7,(T).5));
        frame.r=QUATERNION<T>(pi,VECTOR<T,3>(0,0,1));
        frame.t=position.Value(time);}
    else if(example_number==11){
        INTERPOLATION_CURVE<T,VECTOR<T,3> > position;
        INTERPOLATION_CURVE<T,T> angle;
        position.Add_Control_Point(0,VECTOR<T,3>((T).5,(T).65,(T).5));
        position.Add_Control_Point(0.5,VECTOR<T,3>((T).5,(T).17,(T).5));
        position.Add_Control_Point(3,VECTOR<T,3>((T).5,(T).17,(T).5));
        position.Add_Control_Point(4,VECTOR<T,3>((T).5,(T).65,(T).5));
        angle.Add_Control_Point(1.5,(T)3*pi/4);
        angle.Add_Control_Point(2.5,(T)pi);
        angle.Add_Control_Point(5,(T)pi);
        angle.Add_Control_Point(6,(T)3*pi/8);
        frame.r=QUATERNION<T>(angle.Value(time),VECTOR<T,3>(0,0,1))*QUATERNION<T>(pi/2,VECTOR<T,3>(1,0,0));
        frame.t=position.Value(time);}
    else if(example_number==14){
        INTERPOLATION_CURVE<T,VECTOR<T,3> > position;
        INTERPOLATION_CURVE<T,T> angle;
        T center_offset=.5;Set_Parameter(center_offset,"center_offset");
        position.Add_Control_Point(0,VECTOR<T,3>((T)center_offset,(T).65,(T).5));
        position.Add_Control_Point(0.5,VECTOR<T,3>((T)center_offset,(T).17,(T).5));
        position.Add_Control_Point(3,VECTOR<T,3>((T)center_offset,(T).17,(T).5));
        position.Add_Control_Point(4,VECTOR<T,3>((T)center_offset,(T).8,(T).5));
        angle.Add_Control_Point(1.5,(T)3*pi/4);
        angle.Add_Control_Point(2.5,(T)pi);
        angle.Add_Control_Point(4.5,(T)pi);
        angle.Add_Control_Point(5.5,(T)3*pi/8);
        frame.r=QUATERNION<T>(angle.Value(time_scale*time),VECTOR<T,3>(0,0,1))*QUATERNION<T>(pi/2,VECTOR<T,3>(1,0,0));
        frame.t=position.Value(time_scale*time);}
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time)
{
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Update_Fluid_Parameters(dt,time);
    if(example_number==1 && time>=.583) use_source=false;
    if(example_number==8||example_number==9){use_source=time>4;}
    if(!use_source_override)use_source=false;
    if(source_stop_time!=-1 && time>source_stop_time)use_source=false;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    // Not so good to set up a heaviside function here because then the interface will
    // be exactly between the two nodes which can lead to roundoff issues when setting dirichlet cells, etc.
    GRID<TV>& grid=fluids_parameters.grid;
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) for(int ij=0;ij<grid.mn;ij++)
        fluids_parameters.particle_levelset_evolution.phi(i,j,ij)=grid.y(j)-grid.ymin-initial_water_level;
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
void Adjust_Phi_With_Sources()
{
    if(!use_source) return;
    else if(source_type==BOX_SOURCE) SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Adjust_Phi_With_Source(box_source,world_to_source);
    else if(source_type==CYLINDER_SOURCE) SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Adjust_Phi_With_Source(cylinder_source,world_to_source);
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,3> >*& cell_centered_mask)
{
    if(!use_source) return;
    else if(source_type==BOX_SOURCE) SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Get_Source_Reseed_Mask(box_source,world_to_source,cell_centered_mask,true);
    else if(source_type==CYLINDER_SOURCE) SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Get_Source_Reseed_Mask(cylinder_source,world_to_source,cell_centered_mask,true);
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities()
{
    if(!use_source) return;
    else if(source_type==BOX_SOURCE) SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Get_Source_Velocities(box_source,world_to_source,source_velocity);
    else if(source_type==CYLINDER_SOURCE) SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Get_Source_Velocities(cylinder_source,world_to_source,source_velocity);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time)
{
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Set_Fluid_Boundary_Conditions();
}
//#####################################################################
// Function Adjust_Particle_For_Objects
//#####################################################################
bool Adjust_Particle_For_Objects(PARTICLE_LEVELSET_PARTICLES<T,VECTOR<T,3> >& particles,const int index,VECTOR<T,3>& V,const typename PARTICLE_LEVELSET<T,VECTOR<T,3> >::PARTICLE_TYPE particle_type,const T dt,const T time)
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

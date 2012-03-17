//#####################################################################
// Copyright 2004, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARTICULATED_RIGID_BODIES_2D_EXAMPLE
//##################################################################### 
#ifndef __ARTICULATED_RIGID_BODIES_2D_EXAMPLE__
#define __ARTICULATED_RIGID_BODIES_2D_EXAMPLE__

#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RIGID_BODY_HINTS.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RGB_COLORS.h>
#include <fstream>
//#include "RANDOM_PLACEMENT_2D.h"
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NR0.h>
#include <Rigid_Bodies/RIGID_BODY_INTERSECTIONS_2D.h>
#include <Rigid_Bodies/RIGID_BODY_LIST_2D.h>
#include <Solids_And_Fluids/FLUIDS_PARAMETERS_2D.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_2D.h>
namespace PhysBAM{

template<class T,class RW>
class ARTICULATED_RIGID_BODIES_2D_EXAMPLE : public SOLIDS_FLUIDS_EXAMPLE_2D<T,RW>
{
public:
    RIGID_BODY_LIST_2D<T>* rigid_body_list;
    ARRAY<OPENGL_RIGID_BODY_HINTS> opengl_hints;
    ARTICULATED_RIGID_BODY<TV> arb;
    RANDOM_NR0 random_numbers;
    T ether_viscosity;
    T artificial_maximum_speed; // artificial limit on the maximum possible speed
    bool spatial_partition_based_on_scene_size;
    int spatial_partition_number_of_cells;
    bool spatial_partition_based_on_object_size,spatial_partition_with_max_size;
    bool use_particle_partition,use_particle_partition_center_phi_test;
    int particle_partition_size;
    bool use_triangle_hierarchy,use_triangle_hierarchy_center_phi_test;
    bool use_edge_intersection;
    bool extra_verbose,verbose_dt;
    bool is_arb;
    T max_rotation_per_time_step;
    T max_linear_movement_fraction_per_time_step;

    int parameter;

    ARTICULATED_RIGID_BODIES_2D_EXAMPLE(int parameter=0)
        :SOLIDS_FLUIDS_EXAMPLE_2D<T,RW>(FLUIDS_PARAMETERS_2D<T>::NONE),ether_viscosity(0),artificial_maximum_speed(0),
        spatial_partition_based_on_scene_size(false),spatial_partition_number_of_cells(100),
        spatial_partition_based_on_object_size(true),spatial_partition_with_max_size(true),
        use_particle_partition(true),use_particle_partition_center_phi_test(true),
        particle_partition_size(2),
        use_triangle_hierarchy(false),use_triangle_hierarchy_center_phi_test(false),
        use_edge_intersection(false),extra_verbose(false),verbose_dt(true),is_arb(false),
        max_rotation_per_time_step(0.1*pi),max_linear_movement_fraction_per_time_step(0.1),
        parameter(parameter)
    {
        solids_parameters.cfl=(T).5; // To maintain the old rigid bodies default of .5 instead of the new SIMULATION_EXAMPLE default of .9
        assert(spatial_partition_based_on_scene_size || spatial_partition_based_on_object_size);
        random_numbers.Set_Seed(12345);
        rigid_body_list=&arb.rigid_bodies_list;
    }

    virtual ~ARTICULATED_RIGID_BODIES_2D_EXAMPLE()
    {}

    virtual void Initialize_Rigid_Bodies() {}
    //virtual void Initialize_Bodies() {SOLIDS_FLUIDS_EXAMPLE_2D<T,RW>::Initialize_Bodies();}
    virtual void Update_Animated_Parameters(const T time){}
    virtual void Calculate_ARB_Impulse() {}

//#####################################################################
// Function Initialize_Rigid_Body_Forces
//#####################################################################
virtual void Initialize_Rigid_Body_Forces()
{
    for(int i=0;i<rigid_body_list->rigid_bodies.m;i++)
        if(!rigid_body_list->rigid_bodies(i)->is_static) 
//            rigid_body_list->rigid_bodies(i)->Add_Basic_Forces(gravity,VECTOR_2D<T>(gravity_direction.x,gravity_direction.y),ether_viscosity,0);
            rigid_body_list->rigid_bodies(i)->Add_Basic_Forces(solids_parameters.gravity,VECTOR_2D<T>(solids_parameters.gravity_direction.x,solids_parameters.gravity_direction.y),ether_viscosity,0);
    //temporarily remove gravity
}
//#####################################################################
// Function Load_Restart_Data
//#####################################################################
void Load_Restart_Data(const int frame)
{
    std::cout << "Restarting from frame " << frame << std::endl;
    rigid_body_list->template Read<T>(output_directory,frame);
    for(int i=0;i<rigid_body_list->rigid_bodies.m;i++) rigid_body_list->rigid_bodies(i)->Update_Angular_Velocity();
}
//#####################################################################
// Function Initialize_Rigid_Body
//#####################################################################
RIGID_BODY<TV>* Initialize_Rigid_Body(const std::string& basename,T scaling_factor=1)
{
    std::string fullbasename=data_directory+"/Rigid_Bodies_2D/"+basename;
    int index = rigid_body_list->template Add_Rigid_Body<T>(fullbasename, scaling_factor);

    RIGID_BODY<TV> *rigid_body = rigid_body_list->rigid_bodies(index);

    // If this is a ground plane (based on filename which is kind of a hack) then we fix the bounding box here to only be extended 
    // downwards (rather than ymin-=1, ymax+=1 as is done in ARTICULATED_RIGID_BODIES_2D_DRIVER).  This reduces bounding volume collisions between 
    // objects near the plane and the plane.
    if(basename=="ground"){
        BOX_2D<T>& box=*rigid_body->segmented_curve->bounding_box;
        if(box.ymin == box.ymax){ 
            std::cout << "Fixing ground bounding box" << std::endl;
            std::cout << "before: " << box << std::endl;
            box.ymin-=2;
            std::cout << "after: " << box << std::endl;}}

    // hints
    Append_Hints(basename,scaling_factor);

    return rigid_body;
}
//#####################################################################
// Function Append_Hints
//#####################################################################
void Append_Hints(const std::string& name,const T scale)
{
    VECTOR_3D<double> default_color(1,1,1);
    opengl_hints.Append(OPENGL_RIGID_BODY_HINTS(OPENGL_MATERIAL::Plastic(OPENGL_COLOR(default_color)),true));
}
//#####################################################################
// Function Append_Diffuse_Hints
//#####################################################################
void Append_Diffuse_Hints(const VECTOR_3D<double>& rgb_color,const bool include_bounding_box=true)
{
    opengl_hints(opengl_hints.m).material=OPENGL_MATERIAL::Plastic(OPENGL_COLOR(rgb_color));
    opengl_hints(opengl_hints.m).include_bounding_box = include_bounding_box;
}
//#####################################################################
// Function Write_Data_Files
//#####################################################################
virtual void Write_Data_Files(const int frame)
{
    if(!FILE_UTILITIES::Directory_Exists(output_directory)) FILE_UTILITIES::Create_Directory(output_directory);
    rigid_body_list->template Write<RW>(output_directory,frame);
    Write_Last_Valid_Frame(frame);
}
//#####################################################################
// Function Write_OpenGL_Hints
//#####################################################################
void Write_OpenGL_Hints()
{
    if(opengl_hints.m){
        char filename[256];sprintf(filename,"%s/opengl_hints",output_directory.c_str());std::ofstream output(filename,std::ios::binary);
        opengl_hints.template Write<RW>(output);output.close();}
}
//#####################################################################
// Function Write_Last_Valid_Frame
//#####################################################################
void Write_Last_Valid_Frame(int frame)
{
    char filename[256];sprintf(filename,"%s/last_frame",output_directory.c_str());std::ofstream output(filename,std::ios::trunc);
    output << frame << std::endl;
}
//#####################################################################
// Function Write_Articulated_Points
//#####################################################################
void Write_Articulation_Points(int frame)
{
    arb.Output_Articulation_Points(output_directory,frame);
}
//#####################################################################
};
}
#endif

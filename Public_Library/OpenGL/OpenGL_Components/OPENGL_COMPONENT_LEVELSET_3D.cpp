//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <OpenGL/OpenGL/OPENGL_GRID_3D.h>
#include <OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_3D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_LEVELSET_3D<T>::
OPENGL_COMPONENT_LEVELSET_3D(const VIEWER_DIR& viewer_dir,
                             const std::string& levelset_filename_input,
                             const std::string& triangulated_surface_filename_input,
                             const std::string& filename_set_input,
                             const std::string& filename_triangulated_surface_set_input,
                             bool write_generated_triangulated_surface_input,
                             bool check_triangulated_surface_file_time_input)
    :OPENGL_COMPONENT<T>(viewer_dir,"Levelset 3D"), opengl_levelset_multiview(0), 
      levelset_filename(levelset_filename_input),triangulated_surface_filename(triangulated_surface_filename_input),
      filename_set(filename_set_input),filename_triangulated_surface_set(filename_triangulated_surface_set_input),
      write_generated_triangulated_surface(write_generated_triangulated_surface_input),
      check_triangulated_surface_file_time(check_triangulated_surface_file_time_input),
set(0),set_loaded(-1),use_sets(true),draw_multiple_levelsets(true),ghost_cells(3)
{
    viewer_callbacks.Set("toggle_display_overlay",{[this](){Toggle_Display_Overlay();},"Toggle display overlay (in slice mode)"});
    viewer_callbacks.Set("toggle_slice_color_mode",{[this](){Toggle_Slice_Color_Mode();},"Toggle solid/gradient slice colors"});
    viewer_callbacks.Set("toggle_smooth_slice",{[this](){Toggle_Smooth_Slice();},"Toggle smooth levelset draw"});
    viewer_callbacks.Set("next_set",{[this](){Next_Set();},"Switch to next set"});
    viewer_callbacks.Set("previous_set",{[this](){Previous_Set();},"Switch to previous set"});
    viewer_callbacks.Set("toggle_draw_multiple_levelsets",{[this](){Toggle_Draw_Multiple_Levelsets();},"Toggle mutliple/single levelset draw"});

    int number_of_sets=0;
    while(filename_set!=""){
        std::string filename=LOG::sprintf(filename_set.c_str(),number_of_sets);
        LOG::cout<<"Checking "<<filename<<std::endl;
        if(File_Exists(filename)) number_of_sets++;
        else break;}
    LOG::cout<<"Found "<<number_of_sets<<" levelsets for multiphase"<<std::endl;

    if(number_of_sets==0){use_sets=false;draw_multiple_levelsets=false;number_of_sets=1;}

    opengl_levelset_multiviews.Resize(number_of_sets);
    OPENGL_INDEXED_COLOR_MAP* color_map=OPENGL_INDEXED_COLOR_MAP::Levelset_Multiple_Color_Map();
    for(int i=0;i<opengl_levelset_multiviews.m;i++){
        opengl_levelset_multiviews(i)=new OPENGL_LEVELSET_MULTIVIEW<T>;
        if(use_sets){
            OPENGL_COLOR color=color_map->Lookup(i);
            opengl_levelset_multiviews(i)->Set_Slice_Color(color,OPENGL_COLOR::Transparent());
            opengl_levelset_multiviews(i)->Set_Surface_Material(color,color);
            opengl_levelset_multiviews(i)->Set_Two_Sided(false);}}
    opengl_levelset_multiview=opengl_levelset_multiviews(0);
    delete color_map;

    if(triangulated_surface_filename.length()==0) triangulated_surface_filename="";

    Reinitialize();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT_LEVELSET_3D<T>::
~OPENGL_COMPONENT_LEVELSET_3D()
{
    opengl_levelset_multiviews.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Set_Surface_Material
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Set_Surface_Material(const OPENGL_MATERIAL &front_surface_mat,
                     const OPENGL_MATERIAL &back_surface_mat)
{
    for(int i=0;i<opengl_levelset_multiviews.m;i++)
        opengl_levelset_multiviews(i)->Set_Surface_Material(front_surface_mat, back_surface_mat);
}
//#####################################################################
// Function Set_Overlayed_Surface_Material
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Set_Overlayed_Surface_Material(const OPENGL_MATERIAL &overlayed_surface_mat)
{
    for(int i=0;i<opengl_levelset_multiviews.m;i++)
        opengl_levelset_multiviews(i)->Set_Overlayed_Surface_Material(overlayed_surface_mat);
}
//#####################################################################
// Function Set_Slice_Color
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Set_Slice_Color(const OPENGL_COLOR &inside_slice_color, const OPENGL_COLOR &outside_slice_color)
{
    for(int i=0;i<opengl_levelset_multiviews.m;i++)
        opengl_levelset_multiviews(i)->Set_Slice_Color(inside_slice_color, outside_slice_color);
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Display() const
{
    if(draw){
        if(draw_multiple_levelsets) for(int i=0;i<opengl_levelset_multiviews.m;i++) opengl_levelset_multiviews(i)->Display();
        else opengl_levelset_multiview->Display();}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_LEVELSET_3D<T>::
Bounding_Box() const
{
    if(draw) return opengl_levelset_multiview->Bounding_Box();
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Print_Cell_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Print_Cell_Selection_Info(std::ostream& stream,const TV_INT& cell) const
{
    if(!draw) return;
    bool is_MAC=true;
    if(opengl_levelset_multiviews.Size() && opengl_levelset_multiviews(0)->Levelset() && !opengl_levelset_multiviews(0)->Levelset()->grid.Is_MAC_Grid()) is_MAC=false;
    if(is_MAC){
        TV_INT index=cell;
        opengl_levelset_multiviews(0)->Levelset()->grid.Clamp(index,ghost_cells);
        for(int i=0;i<opengl_levelset_multiviews.m;i++){
            const LEVELSET<TV>& levelset=*opengl_levelset_multiviews(i)->Levelset();
            stream<<component_name<<": phi["<<i<<"]="<<levelset.phi(index)
                  <<" curvature["<<i<<"]="<<levelset.Compute_Curvature(levelset.grid.Center(index))<<std::endl;}}
}
//#####################################################################
// Function Print_Node_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Print_Node_Selection_Info(std::ostream& stream,const TV_INT& node) const
{
    if(!draw) return;
    bool is_MAC=true;
    if(opengl_levelset_multiviews.Size() && opengl_levelset_multiviews(0)->Levelset() && !opengl_levelset_multiviews(0)->Levelset()->grid.Is_MAC_Grid()) is_MAC=false;
    if(!is_MAC){
        TV_INT index=node;
        opengl_levelset_multiviews(0)->Levelset()->grid.Clamp(index,ghost_cells);
        for(int i=0;i<opengl_levelset_multiviews.m;i++)
            if(opengl_levelset_multiviews(i)->Levelset())
                stream<<component_name<<": phi["<<i<<"]="<<(*opengl_levelset_multiviews(i)->Levelset()).phi(index)<<std::endl;}
}
//#####################################################################
// Function Turn_Smooth_Shading_On
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Turn_Smooth_Shading_On()
{
    for(int i=0;i<opengl_levelset_multiviews.m;i++)
        opengl_levelset_multiviews(i)->Turn_Smooth_Shading_On();
}
//#####################################################################
// Function Turn_Smooth_Shading_Off
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Turn_Smooth_Shading_Off()
{
    for(int i=0;i<opengl_levelset_multiviews.m;i++)
        opengl_levelset_multiviews(i)->Turn_Smooth_Shading_Off();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Reinitialize()
{
    if(!draw) return;
    if(use_sets){
        for(int i=0;i<opengl_levelset_multiviews.m;i++)
            Reinitialize_Levelset(viewer_dir.current_directory+"/"+LOG::sprintf(filename_set.c_str(),i),viewer_dir.current_directory+"/"+LOG::sprintf(filename_triangulated_surface_set.c_str(),i),opengl_levelset_multiviews(i));
        set_loaded=set;}
    else Reinitialize_Levelset(viewer_dir.current_directory+"/"+levelset_filename.c_str(),
        viewer_dir.current_directory+"/"+triangulated_surface_filename.c_str(), 
        opengl_levelset_multiview);
}
//#####################################################################
// Function Reinitialize_Levelset
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Reinitialize_Levelset(const std::string& levelset_filename, const std::string& triangulated_surface_filename, OPENGL_LEVELSET_MULTIVIEW<T>* levelset_multiview)
{
    if(File_Exists(levelset_filename)) levelset_multiview->Read_Levelset(levelset_filename);
    else return;
    if(!triangulated_surface_filename.empty() && File_Exists(triangulated_surface_filename) &&
       (!check_triangulated_surface_file_time || Compare_File_Times(triangulated_surface_filename, levelset_filename)>=0)){
        levelset_multiview->Read_Triangulated_Surface(triangulated_surface_filename);}
    else levelset_multiview->Generate_Triangulated_Surface(write_generated_triangulated_surface,triangulated_surface_filename);
    levelset_multiview->Initialize();
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Set_Frame()
{
    
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Toggle_Display_Overlay
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Toggle_Display_Overlay()
{
    for(int i=0;i<opengl_levelset_multiviews.m;i++)
        opengl_levelset_multiviews(i)->Toggle_Display_Overlay();
}
//#####################################################################
// Function Toggle_Slice_Color_Mode
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Toggle_Slice_Color_Mode()
{
    for(int i=0;i<opengl_levelset_multiviews.m;i++)
        opengl_levelset_multiviews(i)->Toggle_Slice_Color_Mode();
}
//#####################################################################
// Function Toggle_Smooth_Slice
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Toggle_Smooth_Slice()
{
    for(int i=0;i<opengl_levelset_multiviews.m;i++)
        opengl_levelset_multiviews(i)->Toggle_Smooth_Slice_Texture();
}
//#####################################################################
// Function Next_Set
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Next_Set()
{
    if(!use_sets) return;
    set=min(set+1,opengl_levelset_multiviews.m-1);
    opengl_levelset_multiview=opengl_levelset_multiviews(set);
    Reinitialize();
}
//#####################################################################
// Function Previous_Set
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Previous_Set()
{
    if(!use_sets) return;
    set=max(set-1,0);
    opengl_levelset_multiview=opengl_levelset_multiviews(set);
    Reinitialize();
}
//#####################################################################
// Function Toggle_Draw_Multiple_Levelsets
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Toggle_Draw_Multiple_Levelsets()
{
    draw_multiple_levelsets=!draw_multiple_levelsets;
}
//#####################################################################
// Function Slice_Has_Changed
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Slice_Has_Changed()
{
    for(int i=0;i<opengl_levelset_multiviews.m;i++)
        opengl_levelset_multiviews(i)->Set_Slice(slice);
}
namespace PhysBAM{
template class OPENGL_COMPONENT_LEVELSET_3D<float>;
template class OPENGL_COMPONENT_LEVELSET_3D<double>;
}

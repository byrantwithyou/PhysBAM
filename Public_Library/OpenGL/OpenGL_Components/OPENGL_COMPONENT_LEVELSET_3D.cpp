//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/LOG.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
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
OPENGL_COMPONENT_LEVELSET_3D(STREAM_TYPE stream_type,const std::string& levelset_filename_input,
                             const std::string& triangulated_surface_filename_input,
                             const std::string& filename_set_input,
                             const std::string& filename_triangulated_surface_set_input,
                             bool write_generated_triangulated_surface_input,
                             bool check_triangulated_surface_file_time_input)
    :OPENGL_COMPONENT<T>(stream_type,"Levelset 3D"), opengl_levelset_multiview(0), 
      levelset_filename(levelset_filename_input),triangulated_surface_filename(triangulated_surface_filename_input),
      filename_set(filename_set_input),filename_triangulated_surface_set(filename_triangulated_surface_set_input),
      write_generated_triangulated_surface(write_generated_triangulated_surface_input),
      frame_loaded(-1),check_triangulated_surface_file_time(check_triangulated_surface_file_time_input),
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
        std::string filename=LOG::sprintf(filename_set.c_str(),frame,number_of_sets);
        LOG::cout<<"Checking "<<filename<<std::endl;
        if(FILE_UTILITIES::File_Exists(filename)) number_of_sets++;
        else break;}
    LOG::cout<<"Found "<<number_of_sets<<" levelsets for multiphase"<<std::endl;

    if(number_of_sets==0){use_sets=false;draw_multiple_levelsets=false;number_of_sets=1;}

    opengl_levelset_multiviews.Resize(number_of_sets);
    OPENGL_INDEXED_COLOR_MAP* color_map=OPENGL_INDEXED_COLOR_MAP::Levelset_Multiple_Color_Map();
    for(int i=0;i<opengl_levelset_multiviews.m;i++){
        opengl_levelset_multiviews(i)=new OPENGL_LEVELSET_MULTIVIEW<T>(stream_type);
        if(use_sets){
            OPENGL_COLOR color=color_map->Lookup(i);
            opengl_levelset_multiviews(i)->Set_Slice_Color(color,OPENGL_COLOR::Transparent());
            opengl_levelset_multiviews(i)->Set_Surface_Material(color,color);
            opengl_levelset_multiviews(i)->Set_Two_Sided(false);}}
    opengl_levelset_multiview=opengl_levelset_multiviews(0);
    delete color_map;

    if(triangulated_surface_filename.length()==0) triangulated_surface_filename="";

    is_animation=levelset_filename.find("%d")!=std::string::npos;
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
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_LEVELSET_3D<T>::
Valid_Frame(int frame_input) const
{
    if(use_sets) return FILE_UTILITIES::File_Exists(LOG::sprintf(filename_set.c_str(),frame,set));
    else return FILE_UTILITIES::File_Exists(is_animation?LOG::sprintf(levelset_filename.c_str(),frame_input):levelset_filename);
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
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION<T>* current_selection) const
{
    if(Is_Up_To_Date(frame)){
        bool is_MAC=true;
        if(opengl_levelset_multiviews.Size() && opengl_levelset_multiviews(0)->Levelset() && !opengl_levelset_multiviews(0)->Levelset()->grid.Is_MAC_Grid()) is_MAC=false;
        if(current_selection && current_selection->type==OPENGL_SELECTION<T>::GRID_CELL_3D && is_MAC){
            VECTOR<int,3> index=((OPENGL_SELECTION_GRID_CELL_3D<T>*)current_selection)->index;
            opengl_levelset_multiviews(0)->Levelset()->grid.Clamp(index,ghost_cells);
            for(int i=0;i<opengl_levelset_multiviews.m;i++){
                const LEVELSET<TV>& levelset=*opengl_levelset_multiviews(i)->Levelset();
                output_stream<<component_name<<": phi["<<i<<"]="<<levelset.phi(index)
                             <<" curvature["<<i<<"]="<<levelset.Compute_Curvature(levelset.grid.Center(index))<<std::endl;}}
        if(current_selection && current_selection->type==OPENGL_SELECTION<T>::GRID_NODE_3D && !is_MAC){
            VECTOR<int,3> index=((OPENGL_SELECTION_GRID_NODE_3D<T>*)current_selection)->index;
            opengl_levelset_multiviews(0)->Levelset()->grid.Clamp(index,ghost_cells);
            for(int i=0;i<opengl_levelset_multiviews.m;i++)  if(opengl_levelset_multiviews(i)->Levelset())
                output_stream<<component_name<<": phi["<<i<<"]="<<(*opengl_levelset_multiviews(i)->Levelset()).phi(index)<<std::endl;}
        if(current_selection && current_selection->type==OPENGL_SELECTION<T>::COMPONENT_PARTICLES_3D){
            OPENGL_SELECTION_COMPONENT_PARTICLES_3D<T> *selection=(OPENGL_SELECTION_COMPONENT_PARTICLES_3D<T>*)current_selection;
            VECTOR<T,3> location=selection->location;
            for(int i=0;i<opengl_levelset_multiviews.m;i++) if(opengl_levelset_multiviews(i)->Levelset())
                output_stream<<component_name<<": phi["<<i<<"] @ particle="<<opengl_levelset_multiviews(i)->Levelset()->Phi(location)<<std::endl;}}
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
    if(draw){
        if((is_animation && (frame_loaded != frame || set_loaded != set)) || (!is_animation && frame_loaded < 0)){
            if(use_sets){
                for(int i=0;i<opengl_levelset_multiviews.m;i++)
                    Reinitialize_Levelset(LOG::sprintf(filename_set.c_str(),frame,i),LOG::sprintf(filename_triangulated_surface_set.c_str(),i,frame),opengl_levelset_multiviews(i));
                set_loaded=set;}
            else Reinitialize_Levelset(FILE_UTILITIES::Get_Frame_Filename(levelset_filename.c_str(),frame), FILE_UTILITIES::Get_Frame_Filename(triangulated_surface_filename.c_str(),frame), opengl_levelset_multiview);
            frame_loaded=frame;
        }
    }
}
//#####################################################################
// Function Reinitialize_Levelset
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Reinitialize_Levelset(const std::string& levelset_filename, const std::string& triangulated_surface_filename, OPENGL_LEVELSET_MULTIVIEW<T>* levelset_multiview)
{
    if(FILE_UTILITIES::File_Exists(levelset_filename)) levelset_multiview->Read_Levelset(levelset_filename);
    else return;
    if(!triangulated_surface_filename.empty() && FILE_UTILITIES::File_Exists(triangulated_surface_filename) &&
       (!check_triangulated_surface_file_time || FILE_UTILITIES::Compare_File_Times(triangulated_surface_filename, levelset_filename)>=0)){
        levelset_multiview->Read_Triangulated_Surface(triangulated_surface_filename);}
    else levelset_multiview->Generate_Triangulated_Surface(write_generated_triangulated_surface,triangulated_surface_filename);
    levelset_multiview->Initialize();
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_3D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
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

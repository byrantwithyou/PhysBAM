//#####################################################################
// Copyright 2004-2005, Eran Guendelman, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL/OPENGL_GRID_2D.h>
#include <OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_LEVELSET_2D<T>::
OPENGL_COMPONENT_LEVELSET_2D(STREAM_TYPE stream_type,const std::string& levelset_filename_input,const std::string filename_set_input)
    :OPENGL_COMPONENT<T>(stream_type,"Levelset 2D"),opengl_levelset(0),levelset_filename(levelset_filename_input),filename_set(filename_set_input),
    frame_loaded(-1),set(0),use_sets(true),set_loaded(-1),valid(false),draw_multiple_levelsets(false)
{
    viewer_callbacks.Set("toggle_color_mode",{[this](){Toggle_Color_Mode();},"Toggle color mode"});
    viewer_callbacks.Set("toggle_smooth",{[this](){Toggle_Smooth();},"Toggle smooth levelset draw"});
    viewer_callbacks.Set("toggle_normals",{[this](){Toggle_Normals();},"Toggle levelset normals draw"});
    viewer_callbacks.Set("toggle_draw_mode",{[this](){Toggle_Draw_Mode();},"Toggle levelset contour/cellview"});
    viewer_callbacks.Set("toggle_draw_sign",{[this](){Toggle_Draw_Sign();},"Toggle levelset sign direction"});
    viewer_callbacks.Set("next_set",{[this](){Next_Set();},"Switch to next set"});
    viewer_callbacks.Set("previous_set",{[this](){Previous_Set();},"Switch to previous set"});
    viewer_callbacks.Set("toggle_draw_multiple_levelsets",{[this](){Toggle_Draw_Multiple_Levelsets();},"Toggle mutliple/single levelset draw"});
    viewer_callbacks.Set("toggle_draw_ghost_values",{[this](){Toggle_Draw_Ghost_Values();},"Toggle draw ghost values"});

    int number_of_sets=0;
    while(filename_set!=""){
        std::string filename=LOG::sprintf(filename_set.c_str(),frame,number_of_sets);
        LOG::cout<<"Checking "<<filename<<std::endl;
        if(File_Exists(filename)) number_of_sets++;else break;}
    LOG::cout<<"Found "<<number_of_sets<<" levelsets for multiphase"<<std::endl;
    if(number_of_sets==0){use_sets=false;number_of_sets=1;}else draw_multiple_levelsets=true;

    opengl_levelsets.Resize(number_of_sets);
    OPENGL_INDEXED_COLOR_MAP* color_map=OPENGL_INDEXED_COLOR_MAP::Levelset_Multiple_Color_Map();
    for(int j=0;j<opengl_levelsets.m;j++)
        opengl_levelsets(j)=new OPENGL_LEVELSET_2D<T>(stream_type,*(new LEVELSET<TV>(*(new GRID<TV>),*(new ARRAY<T,VECTOR<int,2> >))),color_map->Lookup(j),OPENGL_COLOR::Transparent());
    opengl_levelset=opengl_levelsets(0);
    delete color_map;

    is_animation=Is_Animated(levelset_filename);
    Reinitialize();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT_LEVELSET_2D<T>::
~OPENGL_COMPONENT_LEVELSET_2D()
{
    for(int j=0;j<opengl_levelsets.m;j++){
        delete &opengl_levelsets(j)->levelset.grid;
        delete &opengl_levelsets(j)->levelset.phi;
        delete &opengl_levelsets(j)->levelset;
        delete opengl_levelsets(j);}
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_LEVELSET_2D<T>::
Valid_Frame(int frame_input) const
{
    if(use_sets) return File_Exists(LOG::sprintf(filename_set.c_str(),frame_input,set));
    else return Frame_File_Exists(levelset_filename,frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_2D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_2D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_2D<T>::
Display() const
{
    if(valid && draw){
        if(draw_multiple_levelsets) for(int j=0;j<opengl_levelsets.m;j++) opengl_levelsets(j)->Display();
        else opengl_levelsets(set)->Display();}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_LEVELSET_2D<T>::
Bounding_Box() const
{
    if(valid && draw) return opengl_levelset->Bounding_Box();
    return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Print_Cell_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_2D<T>::
Print_Cell_Selection_Info(std::ostream& stream,const TV_INT& cell) const
{
    if(!valid) return;
    if(Is_Up_To_Date(frame)){
        if(opengl_levelsets(0)->levelset.grid.Is_MAC_Grid()){
            for(int i=0;i<opengl_levelsets.m;i++) 
                stream<<component_name<<": phi["<<i<<"]="<<opengl_levelsets(i)->levelset.phi(cell)
                             <<" curvature["<<i<<"]="<<opengl_levelsets(i)->levelset.Compute_Curvature(opengl_levelsets(i)->levelset.grid.Center(cell))<<std::endl;}}
        // if(current_selection && current_selection->type==OPENGL_SELECTION::COMPONENT_PARTICLES_2D){
        //     OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T> *selection=(OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>*)current_selection;
        //     VECTOR<T,2> location=selection->location;
        //     for(int i=0;i<opengl_levelsets.m;i++) 
        //         output_stream<<component_name<<": phi["<<i<<"] @ particle="<<opengl_levelsets(i)->levelset.Phi(location)<<std::endl;}}
}
//#####################################################################
// Function Print_Node_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_2D<T>::
Print_Node_Selection_Info(std::ostream& stream,const TV_INT& node) const
{
    if(!valid) return;
    if(Is_Up_To_Date(frame)){
        if(!opengl_levelsets(0)->levelset.grid.Is_MAC_Grid()){
            for(int i=0;i<opengl_levelsets.m;i++) 
                stream<<component_name<<": phi["<<i<<"]="<<opengl_levelsets(i)->levelset.phi(node)<<std::endl;}}
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_2D<T>::
Reinitialize(const bool force_even_if_not_drawn)
{
    if(draw||force_even_if_not_drawn){
        if((is_animation && (frame_loaded!=frame || set_loaded!=set)) || (!is_animation && frame_loaded<0)){
            valid=false;std::string filename;
            if(use_sets)
                for(int i=0;i<opengl_levelsets.m;i++){
                    filename=LOG::sprintf(filename_set.c_str(),frame,i);
                    if(File_Exists(filename)) Read_From_File(stream_type,filename.c_str(),opengl_levelsets(i)->levelset);
                    else return;
                    opengl_levelsets(i)->Update();}
            else{
                filename=Get_Frame_Filename(levelset_filename,frame);
                if(File_Exists(filename)) Read_From_File(stream_type,filename.c_str(),opengl_levelset->levelset);
                else return;
                opengl_levelset->Update();}
            frame_loaded=frame;set_loaded=set;valid=true;}}
}
//#####################################################################
// Function Toggle_Color_Mode
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_2D<T>::
Toggle_Color_Mode()
{
    for(int j=0;j<opengl_levelsets.m;j++) opengl_levelsets(j)->Toggle_Color_Map();
}
//#####################################################################
// Function Toggle_Smooth
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_2D<T>::
Toggle_Smooth()
{
    for(int j=0;j<opengl_levelsets.m;j++) opengl_levelsets(j)->Toggle_Smooth_Texture();
}
//#####################################################################
// Function Toggle_Normals
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_2D<T>::
Toggle_Normals()
{
    for(int j=0;j<opengl_levelsets.m;j++) opengl_levelsets(j)->Toggle_Normals();
}
//#####################################################################
// Function Toggle_Draw_Mode
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_2D<T>::
Toggle_Draw_Mode()
{
    for(int j=0;j<opengl_levelsets.m;j++){
        int a=opengl_levelsets(j)->draw_area,c=opengl_levelsets(j)->draw_curve,e=opengl_levelsets(j)->draw_cells;
        int mask=c+(a+e*2)*2,newmask=(mask+1)%6;
        opengl_levelsets(j)->draw_curve=(newmask&1)!=0;
        opengl_levelsets(j)->draw_area=(newmask/2)==1;
        opengl_levelsets(j)->draw_cells=(newmask/2)==2;
        opengl_levelsets(j)->Update();}
}
//#####################################################################
// Function Toggle_Draw_Sign
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_2D<T>::
Toggle_Draw_Sign()
{
    for(int j=0;j<opengl_levelsets.m;j++){
        opengl_levelsets(j)->dominant_sign=(opengl_levelsets(j)->dominant_sign==1)?-1:1;
        opengl_levelsets(j)->Update();}
}
//#####################################################################
// Function Next_Set
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_2D<T>::
Next_Set()
{
    set=min(set+1,opengl_levelsets.m-1);
    LOG::cout<<"viewing levelset set "<<set<<std::endl;
    opengl_levelset=opengl_levelsets(set);
    Reinitialize();
}
//#####################################################################
// Function Previous_Set
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_2D<T>::
Previous_Set()
{
    set=max(set-1,0);
    LOG::cout<<"viewing levelset set "<<set<<std::endl;
    opengl_levelset=opengl_levelsets(set);
    Reinitialize();
}
//#####################################################################
// Function Toggle_Draw_Multiple_Levelsets
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_2D<T>::
Toggle_Draw_Multiple_Levelsets()
{
    draw_multiple_levelsets=!draw_multiple_levelsets;
}
//#####################################################################
// Function Toggle_Draw_Ghost_Values
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_2D<T>::
Toggle_Draw_Ghost_Values()
{
    for(int j=0;j<opengl_levelsets.m;j++) opengl_levelsets(j)->Toggle_Draw_Ghost_Values();
}
namespace PhysBAM{
template class OPENGL_COMPONENT_LEVELSET_2D<float>;
template class OPENGL_COMPONENT_LEVELSET_2D<double>;
}

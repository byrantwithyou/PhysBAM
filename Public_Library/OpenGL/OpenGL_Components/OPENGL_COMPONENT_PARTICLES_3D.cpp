//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_3D.h>
#include <sstream>
#include <stdexcept>
using namespace PhysBAM;
//#####################################################################
// Function OPENGL_COMPONENT_PARTICLES_3D
//#####################################################################
template<class T> OPENGL_COMPONENT_PARTICLES_3D<T>::
OPENGL_COMPONENT_PARTICLES_3D(const VIEWER_DIR& viewer_dir,const std::string &filename_input, const std::string &filename_set_input, bool use_ids_input, bool particles_stored_per_cell_uniform_input)
    :OPENGL_COMPONENT<T>(viewer_dir,"Particles 3D"), particles(new GEOMETRY_PARTICLES<TV>),opengl_points(new OPENGL_POINTS_3D<T>(*(new ARRAY<TV>))),
    opengl_vector_field(*(new ARRAY<TV>),opengl_points->points,OPENGL_COLOR::Cyan()),
    filename(filename_input), filename_set(filename_set_input), set(0), set_loaded(-1),number_of_sets(0),use_sets(false),valid(false),
    draw_velocities(false),have_velocities(false),use_ids(use_ids_input),
    particles_stored_per_cell_uniform(particles_stored_per_cell_uniform_input),
    draw_multiple_particle_sets(false)
{
    viewer_callbacks.Set("toggle_draw_point_numbers",{[this](){Toggle_Draw_Point_Numbers();},"Toggle draw point numbers"});
    viewer_callbacks.Set("toggle_draw_velocities",{[this](){Toggle_Draw_Velocities();},"Toggle draw velocities"});
    viewer_callbacks.Set("command_prompt",{[this](){Command_Prompt();},"Command prompt"});
    viewer_callbacks.Set("increase_vector_size",{[this](){Increase_Vector_Size();},"Increase vector size"});
    viewer_callbacks.Set("decrease_vector_size",{[this](){Decrease_Vector_Size();},"Decrease vector size"});
    viewer_callbacks.Set("toggle_arrowhead",{[this](){Toggle_Arrowhead();},"Toggle arrow head style"});
    viewer_callbacks.Set("next_set",{[this](){Next_Set();},"Switch to next set"});
    viewer_callbacks.Set("previous_set",{[this](){Previous_Set();},"Switch to previous set"});
    viewer_callbacks.Set("toggle_draw_multiple_particle_sets",{[this](){Toggle_Draw_Multiple_Particle_Sets();},"Toggle drawing multiple particle sets"});

    number_of_sets=0;
    while(filename_set!=""){
        std::string filename=viewer_dir.current_directory+"/"+LOG::sprintf(filename_set.c_str(),number_of_sets);
        LOG::cout<<"Checking "<<filename<<std::endl;
        if(File_Exists(filename)) number_of_sets++;
        else break;}
    if(number_of_sets>0){use_sets=true;draw_multiple_particle_sets=true;}
    else number_of_sets=1;

    particles_multiple.Resize(number_of_sets);opengl_points_multiple.Resize(number_of_sets);
    particles_multiple(0)=particles;opengl_points_multiple(0)=opengl_points;
    OPENGL_INDEXED_COLOR_MAP* color_map=OPENGL_INDEXED_COLOR_MAP::Particle_Multiple_Color_Map();
    opengl_points->color=color_map->Lookup(0);
    for(int i=1;i<number_of_sets;i++){
        particles_multiple(i)=new GEOMETRY_PARTICLES<TV>;
        opengl_points_multiple(i)=new OPENGL_POINTS_3D<T>(*(new ARRAY<TV>),color_map->Lookup(i));}
    delete color_map;
    // Don't need to call Reinitialize here because it will be called in first call to Set_Frame
}
//#####################################################################
// Function ~OPENGL_COMPONENT_PARTICLES_3D
//#####################################################################
template<class T> OPENGL_COMPONENT_PARTICLES_3D<T>::
~OPENGL_COMPONENT_PARTICLES_3D()
{
    for(int i=0;i<number_of_sets;i++){delete particles_multiple(i);delete &opengl_points_multiple(i)->points;delete opengl_points_multiple(i);}
    //delete &opengl_points->points; TODO: get rid of opengl_points altogether, but this is safe as opengl_points_multiple(0) is always an alias of this.
    delete &opengl_vector_field.vector_field;
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_3D<T>::
Set_Frame()
{
    
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_3D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_3D<T>::
Display() const
{
    if(!valid || !draw) return;
    if(slice && slice->Is_Slice_Mode()){glPushAttrib(GL_ENABLE_BIT);slice->Enable_Clip_Planes();}

    GLint mode;
    glGetIntegerv(GL_RENDER_MODE, &mode);
    if(mode==GL_SELECT){
        glPushName(0);
        glPushName(0);
        if(draw_multiple_particle_sets)
            for(int i=0;i<number_of_sets;i++){
                glLoadName(i);
                opengl_points_multiple(i)->Display();}
        else{
            glLoadName(set);
            opengl_points->Display();}
        glPopName();
        if(draw_velocities && have_velocities){
            glLoadName(1);
            opengl_vector_field.Display();}
        glPopName();}
    else{
        if(draw_multiple_particle_sets) for(int i=0;i<number_of_sets;i++) opengl_points_multiple(i)->Display();
        else opengl_points->Display();
        if(draw_velocities && have_velocities) opengl_vector_field.Display();}

    if(slice && slice->Is_Slice_Mode()) glPopAttrib();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_PARTICLES_3D<T>::
Bounding_Box() const
{
    if(valid && draw) return opengl_points->Bounding_Box();
    return RANGE<TV>::Centered_Box();
}
//#####################################################################
// Function Get_Selection_Priority
//#####################################################################
template<class T> int OPENGL_COMPONENT_PARTICLES_3D<T>::
Get_Selection_Priority(ARRAY_VIEW<GLuint> indices)
{
    if(!indices.m) return -1;
    return 100;
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> bool OPENGL_COMPONENT_PARTICLES_3D<T>::
Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers)
{
    if(indices(0)==0){
        selected_set=indices(1);
        return opengl_points_multiple(selected_set)->Set_Selection(indices.Array_View(2,indices.m-2),modifiers);}
    else if(indices(0)==1)
        return opengl_vector_field.Set_Selection(indices.Array_View(1,indices.m-1),modifiers);
    return false;
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_3D<T>::
Clear_Selection()
{
    opengl_points_multiple(selected_set)->Clear_Selection();
    selected_set=-1;
    opengl_vector_field.Clear_Selection();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_3D<T>::
Print_Selection_Info(std::ostream& output_stream) const
{
    if(!draw_multiple_particle_sets && selected_set!=set) return;
    
    output_stream<<"Selected particle in ["<<component_name<<"("<<selected_set<<")] (total number = "<<particles->Size()<<")"<<std::endl;
    opengl_points_multiple(selected_set)->Print_Selection_Info(output_stream);
}
//#####################################################################
// Destroy_Selection_After_Frame_Change
//#####################################################################
template<class T> bool OPENGL_COMPONENT_PARTICLES_3D<T>::
Destroy_Selection_After_Frame_Change()
{
    if(opengl_points_multiple(selected_set)->Destroy_Selection_After_Frame_Change()){
        selected_set=-1;
        return true;}
    return false;
}
//#####################################################################
// Function Selection_Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_PARTICLES_3D<T>::
Selection_Bounding_Box() const
{
    return opengl_points_multiple(selected_set)->Selection_Bounding_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_3D<T>::
Reinitialize()
{
    if(!draw) return;
    valid=true;
    have_velocities=false;

    for(int i=0;i<number_of_sets;i++){
        std::string frame_filename;
        if(use_sets) frame_filename=viewer_dir.current_directory+"/"+LOG::sprintf(filename_set.c_str(),i);
        else frame_filename=viewer_dir.current_directory+"/"+filename;
        
        try{
            FILE_ISTREAM input;
            Safe_Open_Input(input,frame_filename);
            if(particles_stored_per_cell_uniform){
                ARRAY<GEOMETRY_PARTICLES<TV>*,VECTOR<int,3> > particles_per_cell;
                Read_Binary(input,particles_per_cell);
                ARRAY<PARTICLES<TV>*> initialization_array(particles_per_cell.array.Size());
                for(int j=0;j<particles_per_cell.array.Size();j++){
                    if(particles_per_cell.array(j)) initialization_array(j)=particles_per_cell.array(j);
                    else initialization_array(j)=0;}
                ARRAY_VIEW<const PARTICLES<TV>* const> initialization_array_view(initialization_array.Get_Array_Pointer(),initialization_array.Size());
                particles_multiple(i)->Initialize(initialization_array_view);
                particles_per_cell.Delete_Pointers_And_Clean_Memory();}
            else{
                Read_Binary(input,*particles_multiple(i));}
            opengl_points_multiple(i)->Set_Points_From_Particles(*particles_multiple(i),true);}
        catch(FILESYSTEM_ERROR&){valid=false;}}
#if 0
    if(draw_velocities && particles->store_velocity){
        have_velocities=true;
        opengl_vector_field.vector_field.Resize(particles->Size());
        for(int i=0;i<particles->Size();i++)opengl_vector_field.vector_field(i)=particles->V(i);}
#endif
}
//#####################################################################
// Function Toggle_Draw_Point_Numbers
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_3D<T>::
Toggle_Draw_Point_Numbers()
{
    for(int i=0;i<opengl_points_multiple.m;i++)
        opengl_points_multiple(i)->draw_point_numbers=!opengl_points_multiple(i)->draw_point_numbers;
}
//#####################################################################
// Function Toggle_Draw_Velocities
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_3D<T>::
Toggle_Draw_Velocities()
{
    draw_velocities=!draw_velocities;
    if(draw_velocities) Reinitialize();
}
//#####################################################################
// Function Command_Prompt_Response
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_3D<T>::
Command_Prompt_Response()
{
    if(!OPENGL_WORLD<T>::Singleton()->prompt_response.empty()){
        std::string command;
        std::istringstream sstream(OPENGL_WORLD<T>::Singleton()->prompt_response);
        sstream>>command;
        if(command=="s"){
            int index;
            if(sstream>>index){
                if(ARRAY_VIEW<int>* ids=particles_multiple(set)->template Get_Array<int>("id"))
                    if(!ids->Find(index,index)) return;
                opengl_points->Clear_Selection();
                opengl_points->Select_Point(index);}}}
}
//#####################################################################
// Function Command_Prompt
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_3D<T>::
Command_Prompt()
{
    OPENGL_WORLD<T>::Singleton()->Prompt_User("Command: ",{[this](){Command_Prompt_Response();},""});
}
//#####################################################################
// Function Next_Set
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_3D<T>::
Next_Set()
{
    set=min(set+1,number_of_sets-1);
    LOG::cout<<"viewing particle set "<<set<<std::endl;
    particles=particles_multiple(set);opengl_points=opengl_points_multiple(set);
    Reinitialize();
}
//#####################################################################
// Function Previous_Set
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_3D<T>::
Previous_Set()
{
    set=max(set-1,0);
    LOG::cout<<"viewing particle set "<<set<<std::endl;
    particles=particles_multiple(set);opengl_points=opengl_points_multiple(set);
    Reinitialize();
}
//#####################################################################
// Function Toggle_Draw_Multiple_Particle_Sets
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_3D<T>::
Toggle_Draw_Multiple_Particle_Sets()
{
    draw_multiple_particle_sets=!draw_multiple_particle_sets;
}
//#####################################################################
// Function Get_Particles_Id_Array
//#####################################################################
template<class T> ARRAY_VIEW<int>* OPENGL_COMPONENT_PARTICLES_3D<T>::
Get_Particles_Id_Array(int set_number) const
{
    if(set_number<0) set_number=set;
    ARRAY_VIEW<int>* ids=particles->template Get_Array<int>("id");
    if(ids && ids->Size() && (*ids)(0)) return ids; // A hack to ignore ids if the first one equals zero
    return 0;
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_3D<T>::
Increase_Vector_Size()
{
    opengl_vector_field.Scale_Vector_Size(1.1);
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_3D<T>::
Decrease_Vector_Size()
{
    opengl_vector_field.Scale_Vector_Size(1/1.1);
}
//#####################################################################
// Function Toggle_Arrowhead
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_3D<T>::
Toggle_Arrowhead()
{
    opengl_vector_field.Toggle_Arrowhead_Mode();
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_COMPONENT_PARTICLES_3D<float>;
template class OPENGL_COMPONENT_PARTICLES_3D<double>;
}

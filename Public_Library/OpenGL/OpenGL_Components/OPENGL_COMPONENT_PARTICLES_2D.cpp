//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_2D.h>
#include <sstream>
#include <stdexcept>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_PARTICLES_2D<T>::
OPENGL_COMPONENT_PARTICLES_2D(STREAM_TYPE stream_type,const std::string &filename_input, const std::string &filename_set_input, bool use_ids_input, bool particles_stored_per_cell_uniform_input)
    :OPENGL_COMPONENT<T>(stream_type,"Particles 2D"), particles(new GEOMETRY_PARTICLES<TV>),opengl_points(new OPENGL_POINTS_2D<T>(stream_type,*(new ARRAY<TV>))),
    opengl_vector_field(stream_type,*(new ARRAY<TV>), opengl_points->points, OPENGL_COLOR::Cyan()), 
      filename(filename_input), filename_set(filename_set_input), frame_loaded(-1), set(0), set_loaded(-1),number_of_sets(0),use_sets(false),valid(false),
      draw_velocities(false), have_velocities(false), use_ids(use_ids_input), 
      particles_stored_per_cell_uniform(particles_stored_per_cell_uniform_input),
draw_multiple_particle_sets(false),selected_set(-1)

{
    viewer_callbacks.Set("toggle_draw_point_numbers",{[this](){Toggle_Draw_Point_Numbers();},"Toggle draw point numbers"});
    viewer_callbacks.Set("toggle_draw_radii",{[this](){Toggle_Draw_Radii();},"Toggle draw radii"});
    viewer_callbacks.Set("toggle_draw_velocities",{[this](){Toggle_Draw_Velocities();},"Toggle draw velocities"});
    viewer_callbacks.Set("command_prompt",{[this](){Command_Prompt();},"Command prompt"});
    viewer_callbacks.Set("next_set",{[this](){Next_Set();},"Switch to next set"});
    viewer_callbacks.Set("previous_set",{[this](){Previous_Set();},"Switch to previous set"});
    viewer_callbacks.Set("toggle_draw_multiple_particle_sets",{[this](){Toggle_Draw_Multiple_Particle_Sets();},"Toggle drawing multiple particle sets"});

    number_of_sets=0;
    while(filename_set!=""){
        std::string filename=LOG::sprintf(filename_set.c_str(),frame,number_of_sets);
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
        opengl_points_multiple(i)=new OPENGL_POINTS_2D<T>(stream_type,*(new ARRAY<TV>),color_map->Lookup(i));}
    delete color_map;
        
    is_animation=Is_Animated(filename);
    // Don't need to call Reinitialize here because it will be called in first call to Set_Frame
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT_PARTICLES_2D<T>::
~OPENGL_COMPONENT_PARTICLES_2D()
{
    for(int i=0;i<number_of_sets;i++){delete particles_multiple(i);delete &opengl_points_multiple(i)->points;delete opengl_points_multiple(i);}
    delete &opengl_vector_field.vector_field;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_PARTICLES_2D<T>::
Valid_Frame(int frame_input) const
{
    return Frame_File_Exists(filename,frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Display() const
{
    if(!valid || !draw) return;
    GLint mode=0;
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
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_PARTICLES_2D<T>::
Bounding_Box() const
{
    if(valid && draw) return opengl_points->Bounding_Box();
    return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Get_Selection_Priority
//#####################################################################
template<class T> int OPENGL_COMPONENT_PARTICLES_2D<T>::
Get_Selection_Priority(ARRAY_VIEW<GLuint> indices)
{
    if(!indices.m) return -1;
    return 100;
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> bool OPENGL_COMPONENT_PARTICLES_2D<T>::
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
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Clear_Selection()
{
    opengl_points_multiple(selected_set)->Clear_Selection();
    selected_set=-1;
    opengl_vector_field.Clear_Selection();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Print_Selection_Info(std::ostream &output_stream) const
{
    if(!draw_multiple_particle_sets && selected_set!=set) return;
    
    output_stream<<"Selected particle in ["<<component_name<<"("<<selected_set<<")] (total number = "<<particles->Size()<<")"<<std::endl;
    opengl_points_multiple(selected_set)->Print_Selection_Info(output_stream);
}
//#####################################################################
// Function Destroy_Selection_After_Frame_Change
//#####################################################################
template<class T> bool OPENGL_COMPONENT_PARTICLES_2D<T>::
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
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_PARTICLES_2D<T>::
Selection_Bounding_Box() const
{
    return opengl_points_multiple(selected_set)->Selection_Bounding_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Reinitialize(bool force)
{
    if(!draw || !(force || !valid || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded<0))) return;
    valid=true;
    have_velocities=false;

    for(int i=0;i<number_of_sets;i++){
        std::string frame_filename;
        if(use_sets) frame_filename=LOG::sprintf(filename_set.c_str(),frame,i);
        else frame_filename=Get_Frame_Filename(filename,frame);
        
        try{
            std::istream* input_file=Safe_Open_Input(frame_filename);
            TYPED_ISTREAM typed_input(*input_file,stream_type);
            if(particles_stored_per_cell_uniform){
                ARRAY<GEOMETRY_PARTICLES<TV>*,VECTOR<int,2> > particles_per_cell;
                Read_Binary(typed_input,particles_per_cell);
                ARRAY<PARTICLES<TV>*> initialization_array(particles_per_cell.array.Size());
                for(int j=0;j<particles_per_cell.array.Size();j++){
                    if(particles_per_cell.array(j)) initialization_array(j)=particles_per_cell.array(j);
                    else initialization_array(j)=0;}
                ARRAY_VIEW<const PARTICLES<TV>* const> initialization_array_view(initialization_array.Size(),initialization_array.Get_Array_Pointer());
                particles_multiple(i)->Initialize(initialization_array_view);
                particles_per_cell.Delete_Pointers_And_Clean_Memory();}
            else{
                Read_Binary(typed_input,*particles_multiple(i));}
            delete input_file;
            opengl_points_multiple(i)->Set_Points_From_Particles(*particles_multiple(i),true);}
        catch(FILESYSTEM_ERROR&){valid=false;}}
    frame_loaded=frame;
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
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Toggle_Draw_Point_Numbers()
{
    for(int i=0;i<opengl_points_multiple.m;i++)
        opengl_points_multiple(i)->draw_point_numbers=!opengl_points_multiple(i)->draw_point_numbers;
}
//#####################################################################
// Function Toggle_Draw_Radii
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Toggle_Draw_Radii()
{
    for(int i=0;i<opengl_points_multiple.m;i++)
        opengl_points_multiple(i)->draw_radii=!opengl_points_multiple(i)->draw_radii;
}
//#####################################################################
// Function Toggle_Draw_Velocities
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Toggle_Draw_Velocities()
{
    draw_velocities=!draw_velocities;
    if(draw_velocities) Reinitialize(true);
}
//#####################################################################
// Function Command_Prompt_Response
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
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
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Command_Prompt()
{
    OPENGL_WORLD<T>::Singleton()->Prompt_User("Command: ",{[this](){Command_Prompt_Response();},""});
}
//#####################################################################
// Function Next_Set
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
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
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
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
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Toggle_Draw_Multiple_Particle_Sets()
{
    draw_multiple_particle_sets=!draw_multiple_particle_sets;
}
//#####################################################################
// Function Get_Particles_Id_Array
//#####################################################################
template<class T> ARRAY_VIEW<int>* OPENGL_COMPONENT_PARTICLES_2D<T>::
Get_Particles_Id_Array(int set_number) const
{
    if(set_number<0) set_number=set;
    ARRAY_VIEW<int>* ids=particles_multiple(set_number)->template Get_Array<int>("id");
    if(ids && ids->Size() && (*ids)(0)) return ids; // A hack to ignore ids if the first one equals zero
    return 0;
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_COMPONENT_PARTICLES_2D<float>;
template class OPENGL_COMPONENT_PARTICLES_2D<double>;
}

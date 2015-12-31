//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Log/LOG.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_DEBUG_PARTICLES_3D.h>
#include <sstream>
#include <stdexcept>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
OPENGL_COMPONENT_DEBUG_PARTICLES_3D(STREAM_TYPE stream_type,const std::string &filename_input)
    :OPENGL_COMPONENT<T>(stream_type,"Particles 3D"),particles(*new GEOMETRY_PARTICLES<TV>),debug_objects(*new ARRAY<DEBUG_OBJECT<TV> >),opengl_particles(*new OPENGL_DEBUG_PARTICLES_3D<T>(stream_type,particles,debug_objects)),
    filename(filename_input),frame_loaded(-1),set(0),set_loaded(-1),valid(false),draw_multiple_particle_sets(false)
{
    viewer_callbacks.Set("show_colored_wireframe",{[this](){Show_Colored_Wireframe();},"Show colored wireframe"});
    viewer_callbacks.Set("toggle_draw_velocities",{[this](){Toggle_Draw_Velocities();},"Toggle draw velocities"});
    viewer_callbacks.Set("increase_vector_size",{[this](){Increase_Vector_Size();},"Increase vector size"});
    viewer_callbacks.Set("decrease_vector_size",{[this](){Decrease_Vector_Size();},"Decrease vector size"});
    viewer_callbacks.Set("command_prompt",{[this](){Command_Prompt();},"Command prompt"});

    is_animation=FILE_UTILITIES::Is_Animated(filename);
    // Don't need to call Reinitialize here because it will be called in first call to Set_Frame
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
~OPENGL_COMPONENT_DEBUG_PARTICLES_3D()
{
    delete &particles;
    delete &debug_objects;
    delete &opengl_particles;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename,frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Display() const
{
    if(valid && draw) opengl_particles.Display();
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
template<class T> bool OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Use_Bounding_Box() const
{
    return draw && valid && opengl_particles.Use_Bounding_Box();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Bounding_Box() const
{
    if(valid && draw) return opengl_particles.Bounding_Box();
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    int point_index;
    if(buffer_size==1) point_index=buffer[0];
    else if(buffer_size==2) point_index=buffer[1];
    else return 0;

    OPENGL_SELECTION_COMPONENT_DEBUG_PARTICLES_3D<T> *selection=new OPENGL_SELECTION_COMPONENT_DEBUG_PARTICLES_3D<T>(this);

    // We have the OPENGL_PARTICLES_3D index but need to find the particle index
    selection->index=point_index;
    selection->location=particles.X(point_index);
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Highlight_Selection(OPENGL_SELECTION<T>* selection)
{
    if(selection->type != OPENGL_SELECTION<T>::DEBUG_PARTICLES_3D) return;
    OPENGL_SELECTION_COMPONENT_DEBUG_PARTICLES_3D<T> *real_selection=(OPENGL_SELECTION_COMPONENT_DEBUG_PARTICLES_3D<T>*)selection;
    int particle_index=real_selection->index;
    int point_index=0;
    for(int i=0;i<particle_index;i++)
        point_index++;
    opengl_particles.Select_Point(point_index);
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Clear_Highlight()
{
    opengl_particles.Clear_Selection();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION<T>* selection) const
{
    if(selection && selection->type == OPENGL_SELECTION<T>::DEBUG_PARTICLES_3D && selection->object == this){
        OPENGL_SELECTION_COMPONENT_DEBUG_PARTICLES_3D<T> *real_selection=(OPENGL_SELECTION_COMPONENT_DEBUG_PARTICLES_3D<T>*)selection;
        
        if(!draw_multiple_particle_sets && set!=0) return;
    
        output_stream<<"Selected particle in ["<<component_name<<"("<<1<<")] (total number = "<<particles.Size()<<")"<<std::endl;
    
        int current_index=-1;
        if(real_selection->index<particles.Size()){
            output_stream<<"  Selected by index "<<real_selection->index<<std::endl;
            current_index=real_selection->index;}

        if(current_index>=0){
            // real_selection->index is index into particles array at time of selection.  Not very useful. current_index is more useful
            output_stream<<"current index = "<<current_index<<std::endl;
            particles.Print(output_stream,current_index);}}
}
//#####################################################################
// Function Create_Or_Destroy_Selection_After_Frame_Change
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION<T>* old_selection,bool& delete_selection)
{
    // TODO: reimplement transfering particles between objects.
    if(old_selection && old_selection->object == this){
        delete_selection=true;
    }
    return 0;
}
//#####################################################################
// Function Selection_Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Selection_Bounding_Box(OPENGL_SELECTION<T>* selection) const
{
    int current_index=dynamic_cast<OPENGL_SELECTION_COMPONENT_DEBUG_PARTICLES_3D<T>&>(*selection).index;
    return World_Space_Box(RANGE<TV>(particles.X(current_index)));
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Reinitialize(bool force)
{
    if(!draw || !(force || !valid || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded<0))) return;
    valid=true;

    std::string frame_filename;
    frame_filename=FILE_UTILITIES::Get_Frame_Filename(filename,frame);
        
    try{
        std::istream* input_file=FILE_UTILITIES::Safe_Open_Input(frame_filename);
        TYPED_ISTREAM typed_input(*input_file,stream_type);
        Read_Binary(typed_input,particles,debug_objects);
        delete input_file;}
    catch(FILESYSTEM_ERROR&){valid=false;}
    frame_loaded=frame;
}
//#####################################################################
// Function Show_Colored_Wireframe
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Show_Colored_Wireframe()
{
    opengl_particles.wireframe_only=!opengl_particles.wireframe_only;
}
//#####################################################################
// Function Toggle_Draw_Velocities
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Toggle_Draw_Velocities()
{
    opengl_particles.draw_velocities=!opengl_particles.draw_velocities;
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Increase_Vector_Size()
{
    opengl_particles.scale_velocities*=(T)1.1;
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Decrease_Vector_Size()
{
    opengl_particles.scale_velocities/=(T)1.1;
}
//#####################################################################
// Function Command_Prompt_Response
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Command_Prompt_Response()
{
    if(!OPENGL_WORLD<T>::Singleton()->prompt_response.empty()){
        std::string command;
        std::istringstream sstream(OPENGL_WORLD<T>::Singleton()->prompt_response);
        sstream>>command;
        if(command=="s"){
            ARRAY<int> indices;
            int index;
            while(sstream>>index) indices.Append(index);
            opengl_particles.Clear_Selection();
            opengl_particles.Select_Points(indices);}}
}
//#####################################################################
// Function Set_Slice
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Set_Slice(OPENGL_SLICE *slice_input)
{
    slice=slice_input;
    Slice_Has_Changed();
    opengl_particles.Set_Slice(slice_input);
}
//#####################################################################
// Function Command_Prompt
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Command_Prompt()
{
    OPENGL_WORLD<T>::Singleton()->Prompt_User("Command: ",{[this](){Command_Prompt_Response();},""});
}
//#####################################################################
// Selection object functions
//#####################################################################
template<class T>
RANGE<VECTOR<T,3> > OPENGL_SELECTION_COMPONENT_DEBUG_PARTICLES_3D<T>::
Bounding_Box() const
{
    return object->Selection_Bounding_Box((OPENGL_SELECTION<T>*)this);
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_COMPONENT_DEBUG_PARTICLES_3D<float>;
template class OPENGL_COMPONENT_DEBUG_PARTICLES_3D<double>;
}

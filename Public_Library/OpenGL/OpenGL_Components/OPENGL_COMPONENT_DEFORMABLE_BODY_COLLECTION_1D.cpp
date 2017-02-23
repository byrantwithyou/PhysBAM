//#####################################################################
// Copyright 2009, Nipun Kwatra
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <OpenGL/OpenGL/OPENGL_FREE_PARTICLES.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <OpenGL/OpenGL/OPENGL_TETRAHEDRALIZED_VOLUME.h>
#include <OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D.h>
using namespace PhysBAM;
//#####################################################################
// Function OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D
//#####################################################################
template<class T> OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D(STREAM_TYPE stream_type,const std::string& prefix,const int start_frame)
    :OPENGL_COMPONENT<T>(stream_type,"Deformable Object List"),prefix(prefix),frame_loaded(-1),valid(false),use_active_list(false),display_mode(0),
    incremented_active_object(0),smooth_shading(false),selected_vertex(-1),
    deformable_body_collection(*new DEFORMABLE_BODY_COLLECTION<TV>(0,0)),
    color_map(OPENGL_INDEXED_COLOR_MAP::Basic_16_Color_Map())
{
    viewer_callbacks.Set("toggle_active_value",{[this](){Toggle_Active_Value();},"Toggle viewing of elements"});
    viewer_callbacks.Set("toggle_use_active_list",{[this](){Toggle_Use_Active_List();},"Toggle drawing subset of the deformable objects in the list"});
    viewer_callbacks.Set("toggle_selection_mode",{[this](){Toggle_Selection_Mode();},"Toggle selecting a whole segment or just one part"});
    viewer_callbacks.Set("increment_active_object",{[this](){Increment_Active_Object();},"Increment deformable object being drawn"});
    viewer_callbacks.Set("cycle_display_mode",{[this](){Cycle_Display_Mode();},"Cycle embedded display mode"});
    viewer_callbacks.Set("show_only_first",{[this](){Show_Only_First();},"Show only first deformable object"});
    viewer_callbacks.Set("highlight_particle",{[this](){Highlight_Particle();},"Highlight a particle"});
    viewer_callbacks.Set("toggle_velocity_vectors",{[this](){Toggle_Velocity_Vectors();},"Toggle particle velocity vectors"});
    viewer_callbacks.Set("increase_vector_size",{[this](){Increase_Vector_Size();},"Increase vector size"});
    viewer_callbacks.Set("decrease_vector_size",{[this](){Decrease_Vector_Size();},"Decrease vector size"});
    viewer_callbacks.Set("toggle_active_value_response",{[this](){Toggle_Active_Value_Response();},""});
    viewer_callbacks.Set("highlight_particle_response",{[this](){Highlight_Particle_Response();},""});

    // check for per frame particles
    if(FILE_UTILITIES::File_Exists(LOG::sprintf("%s/%d/deformable_object_structures",prefix.c_str(),start_frame)))
        invalidate_deformable_objects_selection_each_frame=true;
    else invalidate_deformable_objects_selection_each_frame=false;

    is_animation=true;
    Initialize();
}
//#####################################################################
// Function ~OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D
//#####################################################################
template<class T> OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
~OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D()
{
    delete color_map;
    delete &deformable_body_collection;
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Initialize()
{
    draw_velocity_vectors=false;
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Reinitialize(bool force,bool read_geometry)
{
    if(!(draw && (force || (is_animation && frame_loaded!=frame) || (!is_animation && frame_loaded<0)))) return;
    static bool first_time=true;
    std::string frame_string=LOG::sprintf("%s/%d/",prefix.c_str(),frame);
    std::string static_frame_string=frame_string;
    int static_frame=FILE_UTILITIES::File_Exists(frame_string+"deformable_object_structures")?frame:-1;
    bool read_static_variables=static_frame!=-1 || first_time || !deformable_body_collection.structures.m;
    if(read_geometry) deformable_body_collection.Read(stream_type,prefix,prefix,frame,static_frame,read_static_variables,true);
    if(read_static_variables){
        int m=deformable_body_collection.structures.m;active_list.Resize(m,true,true,true);
        point_simplices_1d_objects.Delete_Pointers_And_Clean_Memory();point_simplices_1d_objects.Resize(m);
        for(int i=0;i<m;i++){
            STRUCTURE<TV>* structure=deformable_body_collection.structures(i);
            if(POINT_SIMPLICES_1D<T>* point_simplices_1d=dynamic_cast<POINT_SIMPLICES_1D<T>*>(structure)){
                if(first_time) LOG::cout<<"object "<<i<<": point simplices 1d\n";
                point_simplices_1d_objects(i)=new OPENGL_POINT_SIMPLICES_1D<T>(stream_type,*point_simplices_1d,OPENGL_COLOR((T).5,(T).25,0));
                point_simplices_1d_objects(i)->draw_vertices=true;}
            else{if(first_time) LOG::cout<<"object "<<i<<": object unrecognized at geometry level\n";}}}

    Update_Velocity_Field();
    frame_loaded=frame;
    valid=true;
    first_time=false;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(LOG::sprintf("%s/%d/deformable_object_particles",prefix.c_str(),frame_input));
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    if(draw_input) active_list.Fill(true);
    Reinitialize();
}
//#####################################################################
// Function Draw_All_Objects
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Draw_All_Objects()
{
    Set_Draw(true);
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Display() const
{
    if(!draw || !valid) return;

    glPushName(0);
    glPushName(0);
    for(int i=0;i<point_simplices_1d_objects.m;i++){
        if(!active_list(i)) continue;
        glLoadName(i);
        if(point_simplices_1d_objects(i))
            point_simplices_1d_objects(i)->Display();}
    glPopName();
    glPopName();
    
    if(selected_vertex>=0) OPENGL_SELECTION::Draw_Highlighted_Vertex(deformable_body_collection.particles.X(selected_vertex),selected_vertex);

    //if(draw_velocity_vectors) velocity_field.Display();
}
//#####################################################################
// Function Cycle_Display_Mode
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Cycle_Display_Mode()
{
    display_mode=(display_mode+1)%4;
}
//#####################################################################
// Function Show_Only_First
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Show_Only_First()
{
    active_list.Fill(false);
    active_list(0)=true;
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
template<class T> bool OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Use_Bounding_Box() const
{
    return draw && valid && deformable_body_collection.structures.m>0;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Bounding_Box() const
{
    RANGE<VECTOR<T,3> > box;
    if(draw && valid && deformable_body_collection.structures.m>0){
        for(int i=0;i<point_simplices_1d_objects.m;i++) if(point_simplices_1d_objects(i)) box.Enlarge_To_Include_Box(point_simplices_1d_objects(i)->Bounding_Box());}
    return box;
}
//#####################################################################
// Function Highlight_Particle
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Highlight_Particle()
{
    OPENGL_WORLD<T>::Singleton()->Prompt_User("Enter Particle Number: ",{[this](){Highlight_Particle_Response();},""});
}
//#####################################################################
// Function Toggle_Active_Value
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Toggle_Active_Value()
{
    OPENGL_WORLD<T>::Singleton()->Prompt_User("Enter Component Number: ",{[this](){Toggle_Active_Value_Response();},""});
}
//#####################################################################
// Function Toggle_Use_Active_List
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Toggle_Use_Active_List()
{
    if(active_list.m<1) return;
    use_active_list=!use_active_list;if(!use_active_list) active_list.Fill(true);
    else{active_list.Fill(false);incremented_active_object=1;active_list(incremented_active_object)=true;}
}
//#####################################################################
// Function Select_Segment
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Toggle_Selection_Mode()
{
}
//#####################################################################
// Function Increment_Active_Object
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Increment_Active_Object()
{
    if(active_list.m<1 || !use_active_list) return;
    active_list(incremented_active_object)=false;incremented_active_object++;
    if(incremented_active_object>active_list.m) incremented_active_object=1;
    active_list(incremented_active_object)=true;
}
//#####################################################################
// Function Get_Selection_Priority
//#####################################################################
template<class T> int OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Get_Selection_Priority(ARRAY_VIEW<GLuint> indices)
{
    if(!indices.m) return -1;
    return 10;
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> bool OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers)
{
    if(modifiers&GLUT_ACTIVE_CTRL){
        active_list(indices(1))=false;
        Update_Velocity_Field();
        return false;}

    if(indices(0)==0){
        selected_vertex=indices(1);
        return point_simplices_1d_objects(indices(1))->Set_Selection(indices.Array_View(2,indices.m-2),modifiers);}
    return false;
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Clear_Selection()
{
    selected_vertex=-1;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Print_Selection_Info(std::ostream &output_stream) const
{
    output_stream<<"Deformable object "<<selected_vertex<<std::endl;
    point_simplices_1d_objects(selected_vertex)->Print_Selection_Info(output_stream);
}
//#####################################################################
// Function Toggle_Active_Value_Response
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Toggle_Active_Value_Response()
{
    bool start_val=true;
    if(!OPENGL_WORLD<T>::Singleton()->prompt_response.empty()){
        int index;
        std::istringstream sstream(OPENGL_WORLD<T>::Singleton()->prompt_response);
        sstream>>index;
        if(index>=0) {
            if(active_list.m<=index) active_list.Resize(index+1,true,true,start_val);
            active_list(index)=!active_list(index);}}
}
//#####################################################################
// Function Highlight_Particle_Response
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Highlight_Particle_Response()
{
    if(!OPENGL_WORLD<T>::Singleton()->prompt_response.empty()){
        int index=-1;std::istringstream sstream(OPENGL_WORLD<T>::Singleton()->prompt_response);sstream>>index;
        if(index>=0 && index<deformable_body_collection.particles.Size()) selected_vertex=index;}
    Reinitialize(true);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Selection_Bounding_Box() const
{
    return World_Space_Box(point_simplices_1d_objects(selected_vertex)->Bounding_Box());
}
//#####################################################################
// Function Toggle_Velocity_Vectors
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Toggle_Velocity_Vectors()
{
    draw_velocity_vectors=!draw_velocity_vectors;
    Reinitialize();
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Increase_Vector_Size()
{
    //velocity_field.Scale_Vector_Size(1.1);
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Decrease_Vector_Size()
{
    //velocity_field.Scale_Vector_Size(1/1.1);
}
//#####################################################################
// Function Update_Object_Labels
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Update_Velocity_Field()
{
    if(draw_velocity_vectors){
        positions=deformable_body_collection.particles.X;
        velocity_vectors=deformable_body_collection.particles.V;}
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> bool OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Destroy_Selection_After_Frame_Change()
{
    if(invalidate_deformable_objects_selection_each_frame){
        Clear_Selection();
        return true;}
    return false;
}
//#####################################################################
// Function Turn_Smooth_Shading_On
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Turn_Smooth_Shading_On()
{
    smooth_shading=true;
    for(int i=0;i<point_simplices_1d_objects.m;i++)
        if(point_simplices_1d_objects(i))
            point_simplices_1d_objects(i)->Turn_Smooth_Shading_On();
}
//#####################################################################
// Function Turn_Smooth_Shading_Off
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>::
Turn_Smooth_Shading_Off()
{
    smooth_shading=false;
    for(int i=0;i<point_simplices_1d_objects.m;i++)
        if(point_simplices_1d_objects(i))
            point_simplices_1d_objects(i)->Turn_Smooth_Shading_Off();
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<float>;
template class OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<double>;
}

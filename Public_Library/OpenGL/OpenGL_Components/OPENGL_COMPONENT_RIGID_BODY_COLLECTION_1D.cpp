//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h> // TODO: remove once MUSCLE.cpp exists (workaround for windows compiler bug)
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <OpenGL/OpenGL/OPENGL_AXES.h>
#include <OpenGL/OpenGL/OPENGL_POINT_SIMPLICES_1D.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>::
OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D(STREAM_TYPE stream_type,const std::string& basedir_input)
    :OPENGL_COMPONENT<T>(stream_type,"Rigid Geometry Collection 1D"),basedir(basedir_input),frame_loaded(-1),valid(false),show_object_names(false),output_positions(true),draw_velocity_vectors(false),
    draw_node_velocity_vectors(false),draw_point_simplices(true),
    rigid_body_collection(*new RIGID_BODY_COLLECTION<TV>(0)),current_selection(0),need_destroy_rigid_body_collection(true)
{
    viewer_callbacks.Set("toggle_output_positions",{[this](){Toggle_Output_Positions();},"Toggle output positions"});
    viewer_callbacks.Set("toggle_show_object_names",{[this](){Toggle_Show_Object_Names();},"Toggle show object names"});
    viewer_callbacks.Set("toggle_draw_mode",{[this](){Toggle_Draw_Mode();},"Toggle draw mode"});

    is_animation=true;
    has_init_destroy_information=true;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>::
OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D(STREAM_TYPE stream_type,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const std::string& basedir_input)
    :OPENGL_COMPONENT<T>(stream_type,"Rigid Geometry Collection 1D"),basedir(basedir_input),frame_loaded(-1),valid(false),show_object_names(false),output_positions(true),draw_velocity_vectors(false),
    draw_node_velocity_vectors(false),draw_point_simplices(true),rigid_body_collection(rigid_body_collection),current_selection(0),need_destroy_rigid_body_collection(false)
{
    viewer_callbacks.Set("toggle_output_positions",{[this](){Toggle_Output_Positions();},"Toggle output positions"});
    viewer_callbacks.Set("toggle_show_object_names",{[this](){Toggle_Show_Object_Names();},"Toggle show object names"});
    viewer_callbacks.Set("toggle_draw_mode",{[this](){Toggle_Draw_Mode();},"Toggle draw mode"});

    is_animation=true;
    has_init_destroy_information=true;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>::
~OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D()
{
    if(need_destroy_rigid_body_collection) delete &rigid_body_collection;
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>::
Reinitialize(const bool force,const bool read_geometry)
{
    if(draw && (force || (is_animation && (frame_loaded!=frame)) || (!is_animation && (frame_loaded<0)))){
        valid=false;
        if(!FILE_UTILITIES::File_Exists(LOG::sprintf("%s/%d/rigid_body_particles",basedir.c_str(),frame))) return;
        rigid_body_collection.Read(stream_type,basedir,frame,&needs_init,&needs_destroy);

        if(has_init_destroy_information) for(int i=0;i<needs_destroy.m;i++) Destroy_Geometry(needs_destroy(i));

        int max_number_of_bodies(max(opengl_point_simplices.Size(),rigid_body_collection.rigid_body_particles.Size()));
        // only enlarge array as we read in more geometry to memory
        opengl_point_simplices.Resize(max_number_of_bodies);
        opengl_axes.Resize(max_number_of_bodies);
        draw_object.Resize(max_number_of_bodies);draw_object.Fill(true);
        use_object_bounding_box.Resize(max_number_of_bodies);use_object_bounding_box.Fill(true);

        // Initialize bodies which have become active
        if(has_init_destroy_information) for(int i=0;i<needs_init.m;i++){
            int id=needs_init(i);PHYSBAM_ASSERT(rigid_body_collection.Is_Active(id));
            Create_Geometry(id);}
        else for(int i=0;i<max_number_of_bodies;i++){if(rigid_body_collection.Is_Active(i)) Create_Geometry(i);} // TODO: can we figure out what bodies need_init

        // Update active bodies / remove inactive bodies
        for(int id=0;id<rigid_body_collection.rigid_body_particles.Size();id++){
            if(rigid_body_collection.Is_Active(id)) Update_Geometry(id);
            else Destroy_Geometry(id);}
        for(int id=rigid_body_collection.rigid_body_particles.Size();id<opengl_point_simplices.Size();id++){
            Destroy_Geometry(id);}

        frame_loaded=frame;
        valid=true;}

    Update_Object_Labels();
}
//#####################################################################
// Function Create_Geometry
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>::
Create_Geometry(const int id)
{
    RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(id);

    if(!opengl_axes(id)) opengl_axes(id)=new OPENGL_AXES<T>(stream_type);
    if(rigid_body.simplicial_object && !opengl_point_simplices(id)){
        opengl_point_simplices(id)=new OPENGL_POINT_SIMPLICES_1D<T>(stream_type,*rigid_body.simplicial_object);
        opengl_point_simplices(id)->draw_vertices=true;
        opengl_point_simplices(id)->Enslave_Transform_To(*opengl_axes(id));}
}
//#####################################################################
// Function Update_Geometry
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>::
Update_Geometry(const int id)
{
    if(opengl_axes(id)) *opengl_axes(id)->frame=FRAME<VECTOR<T,3> >(Convert_1d_To_3d(rigid_body_collection.Rigid_Body(id).Frame()));
}
//#####################################################################
// Function Destroy_Geometry
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>::
Destroy_Geometry(const int id)
{
    if(opengl_point_simplices(id)){delete opengl_point_simplices(id);opengl_point_simplices(id)=0;}
    if(opengl_axes(id)){delete opengl_axes(id);opengl_axes(id)=0;}
    draw_object(id)=false;
}
//#####################################################################
// Function Update_Object_Labels
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>::
Update_Object_Labels()
{
    for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++){
        if(draw_object(i)){
            if(opengl_point_simplices(i)){
                if(output_positions){
                    //rigid_body_collection.Rigid_Body(i).Update_Angular_Velocity();
                    opengl_point_simplices(i)->Set_Name(LOG::sprintf("%s <%.3f>",rigid_body_collection.Rigid_Body(i).name.c_str(),rigid_body_collection.rigid_body_particles.frame(i).t.x));}
                else opengl_point_simplices(i)->Set_Name(rigid_body_collection.Rigid_Body(i).name);}}}
    for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++){
        if(draw_object(i)){
            if(opengl_point_simplices(i)){
                if(output_positions){
                    rigid_body_collection.Rigid_Body(i).Update_Angular_Velocity();}}}}
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(LOG::sprintf("%s/%d/rigid_bodies",basedir.c_str(),frame_input));
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>::
Display() const
{
    if(draw){
        glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);

        if(draw_point_simplices){
            for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++){
                if(draw_object(i) && opengl_point_simplices(i)) opengl_point_simplices(i)->Display();}}

        if(show_object_names){
            glColor3f(1,1,1);
            for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++){
                if(draw_object(i) && rigid_body_collection.Rigid_Body(i).name.length()){
                    OpenGL_String(rigid_body_collection.rigid_body_particles.frame(i).t,LOG::sprintf("%s %f",rigid_body_collection.Rigid_Body(i).name.c_str(),rigid_body_collection.rigid_body_particles.twist(i).linear.x));}}}
        glPopAttrib();}
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
template<class T> bool OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>::
Use_Bounding_Box() const
{
    int num_drawable_objects=0;
    for(int i=0;i<opengl_point_simplices.Size();i++)
        if(draw_object(i) && use_object_bounding_box(i))
            num_drawable_objects++;
    return (draw && num_drawable_objects>0);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>::
Bounding_Box() const
{
    RANGE<VECTOR<T,3> > box;
    if(draw)
        for(int i=0;i<opengl_point_simplices.Size();i++)
            if(draw_object(i) && use_object_bounding_box(i) && opengl_point_simplices(i)){
                box=RANGE<VECTOR<T,3> >::Combine(box,opengl_point_simplices(i)->Bounding_Box());}
    return box;
}
//#####################################################################
// Function Set_Draw_Object
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>::
Set_Draw_Object(int i,bool draw_it)
{
    if(i>draw_object.Size()) draw_object.Resize(i);
    draw_object(i)=draw_it;
}
//#####################################################################
// Function Set_Object_Color
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>::
Set_Object_Color(int i,const OPENGL_COLOR &color_input)
{
    if(!opengl_point_simplices(i)) return;
    opengl_point_simplices(i)->color=color_input;
}
//#####################################################################
// Function Get_Draw_Object
//#####################################################################
template<class T> bool OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>::
Get_Draw_Object(int i) const
{
    return draw_object(i);
}
//#####################################################################
// Function Set_Use_Object_Bounding_Box
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>::
Set_Use_Object_Bounding_Box(int i,bool use_it)
{
    use_object_bounding_box(i)=use_it;
}
//#####################################################################
// Function Toggle_Output_Positions
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>::
Toggle_Output_Positions()
{
    output_positions=!output_positions;
    Update_Object_Labels();
}
//#####################################################################
// Function Toggle_Show_Object_Names
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>::
Toggle_Show_Object_Names()
{
    show_object_names=!show_object_names;
}
//#####################################################################
// Function Toggle_Draw_Mode
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>::
Toggle_Draw_Mode()
{
    draw_point_simplices=!draw_point_simplices;
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<float>;
template class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<double>;
}

//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_DEBUG_PARTICLES_2D.h>
#include <sstream>
#include <stdexcept>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
OPENGL_COMPONENT_DEBUG_PARTICLES_2D(STREAM_TYPE stream_type,const std::string &filename_input)
    :OPENGL_COMPONENT<T>(stream_type,"Particles 2D"),particles(*new GEOMETRY_PARTICLES<TV>),
    debug_objects(*new ARRAY<DEBUG_OBJECT<TV> >),default_color(OPENGL_COLOR::White()),
    velocity_color(OPENGL_COLOR(1,(T).078,(T).576)),
    draw_velocities(false),draw_arrows(true),scale_velocities((T).025),
    wireframe_only(false),filename(filename_input),frame_loaded(-1),valid(false),
    selected_index(-1)
{
    viewer_callbacks.Set("show_colored_wireframe",{[this](){Show_Colored_Wireframe();},"Show colored wireframe"});
    viewer_callbacks.Set("toggle_draw_velocities",{[this](){Toggle_Draw_Velocities();},"Toggle draw velocities"});
    viewer_callbacks.Set("increase_vector_size",{[this](){Increase_Vector_Size();},"Increase vector size"});
    viewer_callbacks.Set("decrease_vector_size",{[this](){Decrease_Vector_Size();},"Decrease vector size"});
    viewer_callbacks.Set("toggle_arrowhead",{[this](){Toggle_Arrowhead();},"Toggle arrow heads"});
    viewer_callbacks.Set("command_prompt",{[this](){Command_Prompt();},"Command prompt"});

    is_animation=FILE_UTILITIES::Is_Animated(filename);
    // Don't need to call Reinitialize here because it will be called in first call to Set_Frame
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
~OPENGL_COMPONENT_DEBUG_PARTICLES_2D()
{
    delete &particles;
    delete &debug_objects;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename,frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
Display() const
{
    if(!valid || !draw) return;
    if(slice && slice->Is_Slice_Mode()){
        glPushAttrib(GL_ENABLE_BIT);
        slice->Enable_Clip_Planes();}

    for(int i=0;i<debug_objects.m;i++)
        if(debug_objects(i).draw_vertices)
            for(int i=0;i<debug_objects(i).type;i++)
                OPENGL_SHAPES::Draw_Dot(debug_objects(i).X(i),OPENGL_COLOR(1,0,1),5);

    if(TV::m==2){
        glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT);
        glDisable(GL_LIGHTING);
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,1);}
    glDisable(GL_CULL_FACE);
    if(TV::m==3) glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    if(wireframe_only){
        glPushAttrib(GL_POLYGON_BIT);
        glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);}

    for(int i=0;i<debug_objects.m;i++)
        if(debug_objects(i).type==DEBUG_OBJECT<TV>::triangle){
            OpenGL_Begin(GL_TRIANGLES);
            if(TRIANGLE_2D<T>::Signed_Area(debug_objects(i).X(0),debug_objects(i).X(1),debug_objects(i).X(2))>0)
                OPENGL_COLOR(debug_objects(i).color).Send_To_GL_Pipeline();
            else OPENGL_COLOR(debug_objects(i).bgcolor).Send_To_GL_Pipeline();
            OpenGL_Triangle(debug_objects(i).X(0),debug_objects(i).X(1),debug_objects(i).X(2));
            OpenGL_End();}

    if(wireframe_only) glPopAttrib();

    for(int i=0;i<debug_objects.m;i++)
        if(debug_objects(i).type==DEBUG_OBJECT<TV>::segment){
            if(debug_objects(i).separation && debug_objects(i).bgcolor!=debug_objects(i).color){
                TV t=debug_objects(i).X(1)-debug_objects(i).X(0),n=t.Unit_Orthogonal_Vector()*debug_objects(i).separation;
                OPENGL_SHAPES::Draw_Segment(debug_objects(i).X(0)+n,debug_objects(i).X(1)+n,OPENGL_COLOR(debug_objects(i).color),2);
                OPENGL_SHAPES::Draw_Segment(debug_objects(i).X(0)-n,debug_objects(i).X(1)-n,OPENGL_COLOR(debug_objects(i).bgcolor),2);}
            else OPENGL_SHAPES::Draw_Segment(debug_objects(i).X(0),debug_objects(i).X(1),OPENGL_COLOR(debug_objects(i).color),2);}

    if(TV::m==2) glPopAttrib();


    GLint mode=0;
    glGetIntegerv(GL_RENDER_MODE,&mode);

    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,0);
    glEnable(GL_CULL_FACE);
    glPopMatrix();

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    glPushAttrib(GL_ENABLE_BIT | GL_POINT_BIT);
    glPointSize(5);
    glDisable(GL_LIGHTING);

    glGetIntegerv(GL_RENDER_MODE,&mode);

    ARRAY_VIEW<VECTOR<T,3> >* colors=particles.template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
    ARRAY_VIEW<T>* sizes=particles.template Get_Array<T>(ATTRIBUTE_ID_DISPLAY_SIZE);
    ARRAY_VIEW<TV>* V=particles.template Get_Array<TV>(ATTRIBUTE_ID_V);

    if(draw_velocities && V && mode!=GL_SELECT){
        glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
        glDisable(GL_LIGHTING);
        velocity_color.Send_To_GL_Pipeline();
        OpenGL_Begin(GL_LINES);
        for(int i=0;i<particles.X.m;i++){
            TV X=particles.X(i);
            TV Y=X+(*V)(i)*scale_velocities;
            if(draw_arrows) OPENGL_SHAPES::Draw_Arrow(X,Y);
            else OpenGL_Line(X,Y);}
        OpenGL_End();
        glPopAttrib();}

    if(mode==GL_SELECT) glPushName(0);
    for(int i=0;i<particles.X.m;i++){
        if(mode==GL_SELECT) glLoadName(i);

        if(colors) OPENGL_COLOR((*colors)(i)).Send_To_GL_Pipeline();
        else default_color.Send_To_GL_Pipeline();

        if(sizes && (*sizes)(i)) OPENGL_SHAPES::Draw_Circle(particles.X(i),(*sizes)(i)*scale_velocities,20,false);
        else{
            OpenGL_Begin(GL_POINTS);
            OpenGL_Vertex(particles.X(i));
            OpenGL_End();}}
    if(mode==GL_SELECT) glPopName();

    glPopAttrib();
    glPopMatrix();
    if(slice && slice->Is_Slice_Mode()) glPopAttrib();
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
template<class T> bool OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
Use_Bounding_Box() const
{
    return draw && valid && particles.X.Size()>0;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
Bounding_Box() const
{
    if(valid && draw) return World_Space_Box(RANGE<TV>::Bounding_Box(particles.X));
    return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Get_Selection_Priority
//#####################################################################
template<class T> int OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
Get_Selection_Priority(ARRAY_VIEW<GLuint> indices)
{
    return 110;
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> bool OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers)
{
    selected_index=indices(0);
    return true;
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
Clear_Selection()
{
    selected_index=-1;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
Print_Selection_Info(std::ostream &output_stream) const
{
    if(selected_index>=0){
        output_stream<<"(total number = "<<particles.Size()<<")"<<std::endl;
        output_stream<<"current index = "<<selected_index<<std::endl;
        particles.Print(output_stream,selected_index);}
}
//#####################################################################
// Function Destroy_Selection_After_Frame_Change
//#####################################################################
template<class T> bool OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
Destroy_Selection_After_Frame_Change()
{
    return true;
}
//#####################################################################
// Function Selection_Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
Selection_Bounding_Box() const
{
    return World_Space_Box(RANGE<TV>(particles.X(selected_index)));
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
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
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
Show_Colored_Wireframe()
{
    wireframe_only=!wireframe_only;
}
//#####################################################################
// Function Toggle_Draw_Velocities
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
Toggle_Draw_Velocities()
{
    draw_velocities=!draw_velocities;
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
Increase_Vector_Size()
{
    scale_velocities*=(T)1.1;
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
Decrease_Vector_Size()
{
    scale_velocities/=(T)1.1;
}
//#####################################################################
// Function Toggle_Arrowhead
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
Toggle_Arrowhead()
{
    draw_arrows=!draw_arrows;
}
//#####################################################################
// Function Command_Prompt_Response
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
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
            //Store_Point_Colors(false);
            //Set_Point_Colors(indices,OPENGL_COLOR::Yellow());
        }}
}
//#####################################################################
// Function Set_Slice
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
Set_Slice(OPENGL_SLICE *slice_input)
{
    slice=slice_input;
    Slice_Has_Changed();
}
//#####################################################################
// Function Command_Prompt
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>::
Command_Prompt()
{
    OPENGL_WORLD<T>::Singleton()->Prompt_User("Command: ",{[this](){Command_Prompt_Response();},""});
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_COMPONENT_DEBUG_PARTICLES_2D<float>;
template class OPENGL_COMPONENT_DEBUG_PARTICLES_2D<double>;
}

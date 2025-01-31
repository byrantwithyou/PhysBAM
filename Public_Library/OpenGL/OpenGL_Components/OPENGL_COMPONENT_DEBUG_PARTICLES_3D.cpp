//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_DEBUG_PARTICLES_3D.h>
#include <sstream>
#include <stdexcept>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
OPENGL_COMPONENT_DEBUG_PARTICLES_3D(const VIEWER_DIR& viewer_dir,const std::string& filename)
    :OPENGL_COMPONENT<T>(viewer_dir,"Particles 3D"),particles(*new GEOMETRY_PARTICLES<TV>),
    debug_objects(*new ARRAY<DEBUG_OBJECT<TV> >),
    debug_text(*new ARRAY<DEBUG_TEXT<TV> >),default_color(OPENGL_COLOR::White()),
    velocity_color(OPENGL_COLOR(1,(T).078,(T).576)),filename(filename)
{
    viewer_callbacks.Set("show_colored_wireframe",{[this](){Show_Colored_Wireframe();},"Show colored wireframe"});
    viewer_callbacks.Set("toggle_draw_velocities",{[this](){Toggle_Draw_Velocities();},"Toggle draw velocities"});
    viewer_callbacks.Set("increase_vector_size",{[this](){Increase_Vector_Size();},"Increase vector size"});
    viewer_callbacks.Set("decrease_vector_size",{[this](){Decrease_Vector_Size();},"Decrease vector size"});
    viewer_callbacks.Set("toggle_arrowhead",{[this](){Toggle_Arrowhead();},"Toggle arrow heads"});
    viewer_callbacks.Set("command_prompt",{[this](){Command_Prompt();},"Command prompt"});

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
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Set_Frame()
{
    
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
    if(!valid || !draw) return;

    GLint mode=0;
    glGetIntegerv(GL_RENDER_MODE,&mode);

    if(mode==GL_SELECT) glPushName(select_none);
    if(slice && slice->Is_Slice_Mode()){
        glPushAttrib(GL_ENABLE_BIT);
        slice->Enable_Clip_Planes();}

    if(mode==GL_SELECT) glLoadName(select_object);
    if(mode==GL_SELECT) glPushName(0);
    for(int i=0;i<debug_objects.m;i++)
        if(debug_objects(i).draw_vertices){
            if(mode==GL_SELECT) glLoadName(i);
            for(int j=0;j<debug_objects(i).type;j++)
                OPENGL_SHAPES::Draw_Dot(debug_objects(i).X(j),OPENGL_COLOR(1,0,1),5);}
    if(mode==GL_SELECT) glPopName();
    
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

    if(mode==GL_SELECT) glPushName(0);
    for(int i=0;i<debug_objects.m;i++)
        if(debug_objects(i).type==DEBUG_OBJECT<TV>::triangle){
            if(mode==GL_SELECT) glLoadName(i);
            OpenGL_Begin(GL_TRIANGLES);
            OPENGL_MATERIAL::Plastic(OPENGL_COLOR(debug_objects(i).color)).Send_To_GL_Pipeline(GL_FRONT);
            OPENGL_MATERIAL::Plastic(OPENGL_COLOR(debug_objects(i).bgcolor)).Send_To_GL_Pipeline(GL_BACK);
            OpenGL_Triangle(debug_objects(i).X(0),debug_objects(i).X(1),debug_objects(i).X(2));
            OpenGL_End();}
    if(mode==GL_SELECT) glPopName();

    if(wireframe_only) glPopAttrib();

    if(mode==GL_SELECT) glPushName(0);
    for(int i=0;i<debug_objects.m;i++)
        if(debug_objects(i).type==DEBUG_OBJECT<TV>::segment){
            if(mode==GL_SELECT) glLoadName(i);
            int lw=(i==selected_object)?4:2;
            if(debug_objects(i).separation && debug_objects(i).bgcolor!=debug_objects(i).color){
                TV t=debug_objects(i).X(1)-debug_objects(i).X(0),n=t.Unit_Orthogonal_Vector()*debug_objects(i).separation;
                OPENGL_SHAPES::Draw_Segment(debug_objects(i).X(0)+n,debug_objects(i).X(1)+n,OPENGL_COLOR(debug_objects(i).color),lw);
                OPENGL_SHAPES::Draw_Segment(debug_objects(i).X(0)-n,debug_objects(i).X(1)-n,OPENGL_COLOR(debug_objects(i).bgcolor),lw);}
            else OPENGL_SHAPES::Draw_Segment(debug_objects(i).X(0),debug_objects(i).X(1),OPENGL_COLOR(debug_objects(i).color),lw);}
    if(mode==GL_SELECT) glPopName();

    if(TV::m==2) glPopAttrib();

    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_LIGHTING);
    if(mode==GL_SELECT) glLoadName(select_text);
    if(mode==GL_SELECT) glPushName(0);
    for(int i=0;i<debug_text.m;i++){
        if(mode==GL_SELECT) glLoadName(i);
        OPENGL_COLOR(debug_text(i).color).Send_To_GL_Pipeline();
        OpenGL_String(debug_text(i).X,debug_text(i).text,GLUT_BITMAP_HELVETICA_12);}
    if(mode==GL_SELECT) glPopName();
    glPopAttrib();

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

    ARRAY_VIEW<VECTOR<T,3> >* colors=particles.template Get_Array<VECTOR<T,3> >("color");
    ARRAY_VIEW<T>* sizes=particles.template Get_Array<T>("display_size");
    ARRAY_VIEW<TV>* V=particles.template Get_Array<TV>("V");

    if(mode==GL_SELECT) glLoadName(select_particle);
    if(draw_velocities && V && mode!=GL_SELECT){
        glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
        glDisable(GL_LIGHTING);
        velocity_color.Send_To_GL_Pipeline();
        OpenGL_Begin(GL_LINES);
        if(mode==GL_SELECT) glPushName(0);
        for(int i=0;i<particles.X.m;i++){
            if(mode==GL_SELECT) glLoadName(i);
            TV X=particles.X(i);
            TV Y=X+(*V)(i)*scale_velocities;
            if(draw_arrows) OPENGL_SHAPES::Draw_Arrow(X,Y);
            else OpenGL_Line(X,Y);}
        if(mode==GL_SELECT) glPopName();
        OpenGL_End();
        glPopAttrib();}

    if(mode==GL_SELECT) glPushName(0);
    for(int i=0;i<particles.X.m;i++){
        if(mode==GL_SELECT) glLoadName(i);

        if(colors) OPENGL_COLOR((*colors)(i)).Send_To_GL_Pipeline();
        else default_color.Send_To_GL_Pipeline();

        if(sizes && (*sizes)(i)) OPENGL_SHAPES::Draw_Circle(particles.X(i),(*sizes)(i)*scale_velocities,20,false);
        else{
            if(i==selected_particle) glPointSize(8);
            OpenGL_Begin(GL_POINTS);
            OpenGL_Vertex(particles.X(i));
            OpenGL_End();
            if(i==selected_particle) glPointSize(5);}}
    if(mode==GL_SELECT) glPopName();

    glPopAttrib();
    glPopMatrix();
    if(slice && slice->Is_Slice_Mode()) glPopAttrib();
    if(mode==GL_SELECT) glPopName();
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
template<class T> bool OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Use_Bounding_Box() const
{
    return draw && valid && particles.X.Size()>0;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Bounding_Box() const
{
    if(valid && draw) return World_Space_Box(RANGE<TV>::Bounding_Box(particles.X));
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Get_Selection_Priority
//#####################################################################
template<class T> int OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Get_Selection_Priority(ARRAY_VIEW<GLuint> indices)
{
    if(indices(0)==select_particle) return 110;
    if(indices(0)==select_object) return 109;
    if(indices(0)==select_text) return 108;
    return 0;
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> bool OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers)
{
    PHYSBAM_ASSERT(indices.m==2 && indices(0)!=select_none);
    selected_particle=selected_object=selected_text=-1;
    if(indices(0)==select_particle) selected_particle=indices(1);
    if(indices(0)==select_object) selected_object=indices(1);
    if(indices(0)==select_text) selected_text=indices(1);
    return true;
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Clear_Selection()
{
    selected_particle=selected_object=selected_text=-1;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Print_Selection_Info(std::ostream &output_stream) const
{
    if(selected_particle>=0){
        output_stream<<"debug particle"<<std::endl;
        output_stream<<"(total number = "<<particles.Size()<<")"<<std::endl;
        output_stream<<"current index = "<<selected_particle<<std::endl;
        particles.Print(output_stream,selected_particle);}
    if(selected_object>=0){
        output_stream<<"debug object"<<std::endl;
        output_stream<<"(total number = "<<debug_objects.Size()<<")"<<std::endl;
        output_stream<<"current index = "<<selected_object<<std::endl;
        const DEBUG_OBJECT<TV>& obj=debug_objects(selected_object);
        output_stream<<"type = "<<(obj.type==2?"segment":"triangle")<<std::endl;
        output_stream<<"vertices = "<<obj.X.Array_View(0,(int)obj.type)<<std::endl;}
    if(selected_text>=0){
        output_stream<<"debug text"<<std::endl;
        output_stream<<"(total number = "<<debug_text.Size()<<")"<<std::endl;
        output_stream<<"current index = "<<selected_text<<std::endl;
        const DEBUG_TEXT<TV>& text=debug_text(selected_text);
        output_stream<<"string = "<<text.text<<std::endl;
        output_stream<<"position = "<<text.X<<std::endl;}
}
//#####################################################################
// Function Destroy_Selection_After_Frame_Change
//#####################################################################
template<class T> bool OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Destroy_Selection_After_Frame_Change()
{
    return true;
}
//#####################################################################
// Function Selection_Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Selection_Bounding_Box() const
{
    RANGE<TV> box;
    if(selected_particle>=0)
        box.Enlarge_To_Include_Point(particles.X(selected_particle));
    if(selected_object>=0){
        const DEBUG_OBJECT<TV>& obj=debug_objects(selected_object);
        box.Enlarge_To_Include_Points(obj.X.Array_View(0,(int)obj.type));}
    if(selected_text>=0)
        box.Enlarge_To_Include_Point(debug_text(selected_text).X);

    return World_Space_Box(box);
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Reinitialize()
{
    if(!draw) return;
    valid=true;

    std::string frame_filename=viewer_dir.current_directory+"/"+filename;
        
    try{
        FILE_ISTREAM input;
        Safe_Open_Input(input,frame_filename);
        Read_Binary(input,particles,debug_objects,debug_text);}
    catch(FILESYSTEM_ERROR&){valid=false;}
}
//#####################################################################
// Function Show_Colored_Wireframe
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Show_Colored_Wireframe()
{
    wireframe_only=!wireframe_only;
}
//#####################################################################
// Function Toggle_Draw_Velocities
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Toggle_Draw_Velocities()
{
    draw_velocities=!draw_velocities;
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Increase_Vector_Size()
{
    scale_velocities*=(T)1.1;
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Decrease_Vector_Size()
{
    scale_velocities/=(T)1.1;
}
//#####################################################################
// Function Toggle_Arrowhead
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Toggle_Arrowhead()
{
    draw_arrows=!draw_arrows;
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
            //Store_Point_Colors(false);
            //Set_Point_Colors(indices,OPENGL_COLOR::Yellow());
        }}
}
//#####################################################################
// Function Set_Slice
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>::
Set_Slice(OPENGL_SLICE *slice_input)
{
    slice=slice_input;
    Slice_Has_Changed();
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
namespace PhysBAM{
template class OPENGL_COMPONENT_DEBUG_PARTICLES_3D<float>;
template class OPENGL_COMPONENT_DEBUG_PARTICLES_3D<double>;
}

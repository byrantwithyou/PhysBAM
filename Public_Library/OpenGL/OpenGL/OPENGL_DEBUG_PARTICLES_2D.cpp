//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <OpenGL/OpenGL/OPENGL_DEBUG_PARTICLES_2D.h>
#include <OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_DEBUG_PARTICLES_2D<T>::
OPENGL_DEBUG_PARTICLES_2D(STREAM_TYPE stream_type,GEOMETRY_PARTICLES<TV>& particle_input,ARRAY<DEBUG_OBJECT<TV> >& debug_objects_input,const OPENGL_COLOR& color_input)
    :OPENGL_OBJECT<T>(stream_type),particles(particle_input),debug_objects(debug_objects_input),default_color(color_input),velocity_color(OPENGL_COLOR(1,(T).078,(T).576)),
    draw_velocities(false),draw_arrows(true),scale_velocities((T).025)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_DEBUG_PARTICLES_2D<T>::
~OPENGL_DEBUG_PARTICLES_2D()
{
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
template<class T> bool OPENGL_DEBUG_PARTICLES_2D<T>::
Use_Bounding_Box() const
{
    return particles.X.Size()>0;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_DEBUG_PARTICLES_2D<T>::
Bounding_Box() const
{
    return World_Space_Box(RANGE<TV>::Bounding_Box(particles.X));
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_DEBUG_PARTICLES_2D<T>::
Display() const
{
    for(int i=0;i<debug_objects.m;i++)
        if(debug_objects(i).draw_vertices)
            for(int i=0;i<debug_objects(i).type;i++)
                OPENGL_SHAPES::Draw_Dot(debug_objects(i).X(i),OPENGL_COLOR(1,0,1),5);

    glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,1);
    glDisable(GL_CULL_FACE);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    for(int i=0;i<debug_objects.m;i++)
        if(debug_objects(i).type==DEBUG_OBJECT<TV>::triangle){
            OpenGL_Begin(GL_TRIANGLES);
            OPENGL_MATERIAL::Plastic(OPENGL_COLOR(debug_objects(i).color)).Send_To_GL_Pipeline(GL_FRONT);
            OPENGL_MATERIAL::Plastic(OPENGL_COLOR(debug_objects(i).bgcolor)).Send_To_GL_Pipeline(GL_BACK);
            OpenGL_Triangle(debug_objects(i).X(0),debug_objects(i).X(1),debug_objects(i).X(2));
            OpenGL_End();}

    for(int i=0;i<debug_objects.m;i++)
        if(debug_objects(i).type==DEBUG_OBJECT<TV>::segment){
            if(debug_objects(i).separation && debug_objects(i).bgcolor!=debug_objects(i).color){
                TV t=debug_objects(i).X(1)-debug_objects(i).X(0),n=t.Unit_Orthogonal_Vector()*debug_objects(i).separation;
                OPENGL_SHAPES::Draw_Segment(debug_objects(i).X(0)+n,debug_objects(i).X(1)+n,OPENGL_COLOR(debug_objects(i).color),2);
                OPENGL_SHAPES::Draw_Segment(debug_objects(i).X(0)-n,debug_objects(i).X(1)-n,OPENGL_COLOR(debug_objects(i).bgcolor),2);}
            else OPENGL_SHAPES::Draw_Segment(debug_objects(i).X(0),debug_objects(i).X(1),OPENGL_COLOR(debug_objects(i).color),2);}

    glPopAttrib();


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
}
//#####################################################################
// Function Select_Point
//#####################################################################
template<class T> void OPENGL_DEBUG_PARTICLES_2D<T>::
Select_Point(int index)
{
    //Set_Point_Color(index,OPENGL_COLOR::Yellow());
}
//#####################################################################
// Function Select_Points
//#####################################################################
template<class T> void OPENGL_DEBUG_PARTICLES_2D<T>::
Select_Points(const ARRAY<int> &indices)
{
    //Set_Point_Colors(indices,OPENGL_COLOR::Yellow());
}
//#####################################################################
// Function Clear_Selection
//#####################################################################
template<class T> void OPENGL_DEBUG_PARTICLES_2D<T>::
Clear_Selection()
{
    //Store_Point_Colors(false);
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_DEBUG_PARTICLES_2D<T>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    LOG::cout<<"Get_Selection "<<buffer_size<<std::endl;
    if(buffer_size==1){
        OPENGL_SELECTION_DEBUG_PARTICLES_2D<T> *selection=new OPENGL_SELECTION_DEBUG_PARTICLES_2D<T>(this);
        selection->index=buffer[0];
        return selection;}
    else return 0;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_DEBUG_PARTICLES_2D<T>::
Highlight_Selection(OPENGL_SELECTION<T>* selection)
{
    if(selection->type!=OPENGL_SELECTION<T>::DEBUG_PARTICLES_2D) return;
    OPENGL_SELECTION_DEBUG_PARTICLES_2D<T> *real_selection=(OPENGL_SELECTION_DEBUG_PARTICLES_2D<T>*)selection;
    Select_Point(real_selection->index);
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_DEBUG_PARTICLES_2D<T>::
Clear_Highlight()
{
    Clear_Selection();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_SELECTION_DEBUG_PARTICLES_2D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    return object->Bounding_Box();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_DEBUG_PARTICLES_2D<T>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION<T>* selection) const
{
    if(selection->type!=OPENGL_SELECTION<T>::DEBUG_PARTICLES_2D) return;
    output_stream<<"Particle "<<dynamic_cast<OPENGL_SELECTION_DEBUG_PARTICLES_2D<T>&>(*selection).index<<std::endl;
}
namespace PhysBAM{
template class OPENGL_DEBUG_PARTICLES_2D<float>;
template class OPENGL_DEBUG_PARTICLES_2D<double>;
}


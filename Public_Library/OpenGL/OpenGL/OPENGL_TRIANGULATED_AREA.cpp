//#####################################################################
// Copyright 2004, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Matrices/MATRIX_3X3.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <OpenGL/OpenGL/OPENGL_TRIANGULATED_AREA.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_TRIANGULATED_AREA<T>::
OPENGL_TRIANGULATED_AREA(STREAM_TYPE stream_type,TRIANGULATED_AREA<T>& triangulated_area_input,
    const bool draw_vertices_input,const OPENGL_COLOR& vertex_color_input,const OPENGL_COLOR& segment_color_input,
    const OPENGL_COLOR& triangle_color_input,const OPENGL_COLOR& triangle_inverted_color_input,
    ARRAY<OPENGL_COLOR>* color_map_input)
    :OPENGL_OBJECT<T>(stream_type),triangulated_area(triangulated_area_input),vertex_color(vertex_color_input),
    segment_color(segment_color_input),triangle_color(triangle_color_input),
    triangle_inverted_color(triangle_inverted_color_input),velocity_color(OPENGL_COLOR::Yellow()),selected_vertex(-1),
    selected_segment(-1),selected_triangle(-1),color_map(color_map_input),draw_vertices(draw_vertices_input),
    draw_velocities(false),velocity_scale(0.025)
{}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA<T>::
Display() const
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    GLint mode=0;
    glGetIntegerv(GL_RENDER_MODE,&mode);

    if(mode == GL_SELECT){
        glPushName(0);
        Draw_Vertices_For_Selection();
        if(triangulated_area.mesh.segment_mesh){
            glLoadName(1);
            Draw_Segments_For_Selection();}
        glLoadName(2);
        Draw_Triangles_For_Selection();
        glPopName();}
    else
    {
        glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_POLYGON_BIT);
        glDisable(GL_LIGHTING);
        glEnable(GL_POLYGON_OFFSET_FILL);glEnable(GL_POLYGON_OFFSET_LINE);
        glPolygonOffset(0,2);
        Draw_Triangles();
        glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
        glPolygonOffset(0,1);
        segment_color.Send_To_GL_Pipeline();
        Draw_Triangles(false);
        Draw_Vertices();
        glPopAttrib();}

    if(mode != GL_SELECT){
        if(selected_vertex>=0 || selected_segment>=0 || selected_triangle>=0){
            glPushAttrib(GL_ENABLE_BIT);
            glDisable(GL_DEPTH_TEST);
            if(selected_vertex>=0){
                OPENGL_SELECTION::Draw_Highlighted_Vertex(triangulated_area.particles.X(selected_vertex),selected_vertex);}
            else if(selected_segment>=0){
                int node1,node2;triangulated_area.mesh.segment_mesh->elements(selected_segment).Get(node1,node2);
                OPENGL_SELECTION::Draw_Highlighted_Segment(triangulated_area.particles.X(node1),triangulated_area.particles.X(node2),selected_segment);}
            else if(selected_triangle>=0){
                int node1,node2,node3;triangulated_area.mesh.elements(selected_triangle).Get(node1,node2,node3);
                OPENGL_SELECTION::Draw_Highlighted_Triangle_Boundary(triangulated_area.particles.X(node1),triangulated_area.particles.X(node2),triangulated_area.particles.X(node3),selected_triangle);}
            glPopAttrib();
        }
    }
    if(draw_velocities && triangulated_area.particles.store_velocity){
        glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT);
        glDisable(GL_LIGHTING);
        velocity_color.Send_To_GL_Pipeline();
        OpenGL_Begin(GL_LINES);
        for(int t=0;t<triangulated_area.particles.Size();t++)
            OPENGL_SHAPES::Draw_Arrow(triangulated_area.particles.X(t),triangulated_area.particles.X(t)+velocity_scale*triangulated_area.particles.V(t));
        OpenGL_End();
        glPopAttrib();}

    glPopMatrix();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_TRIANGULATED_AREA<T>::
Bounding_Box() const
{
    RANGE<VECTOR<T,3> > box;
    for(int i=0;i<triangulated_area.particles.Size();i++)
        box.Enlarge_To_Include_Point(World_Space_Point(triangulated_area.particles.X(i)));
    return box;
}
//#####################################################################
// Function Get_Selection_Priority
//#####################################################################
template<class T> int OPENGL_TRIANGULATED_AREA<T>::
Get_Selection_Priority(ARRAY_VIEW<GLuint> indices)
{
    if(!indices.m) return -1;
    PHYSBAM_ASSERT(indices.m==2);
    const static int priority[]={77,76,75};
    return priority[indices(0)];
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> bool OPENGL_TRIANGULATED_AREA<T>::
Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers)
{
    if(indices(0)==0) selected_vertex=indices(1);
    else if(indices(0)==1) selected_segment=indices(1);
    else if(indices(0)==2) selected_triangle=indices(1);
    else return false;
    return true;
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA<T>::
Clear_Selection()
{
    selected_vertex=-1;
    selected_segment=-1;
    selected_triangle=-1;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA<T>::
Print_Selection_Info(std::ostream &output_stream) const
{
    Print_Selection_Info(output_stream,0);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA<T>::
Print_Selection_Info(std::ostream &output_stream,MATRIX<T,3>* transform) const
{
    if(selected_vertex>=0){
        output_stream<<"Vertex "<<selected_vertex<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(triangulated_area.particles.X(selected_vertex))<<std::endl;}
        triangulated_area.particles.Print(output_stream,selected_vertex);}
    else if(selected_segment>=0){
        PHYSBAM_ASSERT(triangulated_area.mesh.segment_mesh);
        int node1,node2;triangulated_area.mesh.segment_mesh->elements(selected_segment).Get(node1,node2);
        output_stream<<"Segment "<<selected_segment<<" ("<<node1<<", "<<node2<<")"<<std::endl;
        output_stream<<std::endl;
        output_stream<<"Vertex "<<node1<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(triangulated_area.particles.X(node1))<<std::endl;}
        triangulated_area.particles.Print(output_stream,node1);
        output_stream<<std::endl;
        output_stream<<"Vertex "<<node2<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(triangulated_area.particles.X(node2))<<std::endl;}
        triangulated_area.particles.Print(output_stream,node2);}
    else if(selected_triangle>=0){
        int node1,node2,node3;triangulated_area.mesh.elements(selected_triangle).Get(node1,node2,node3);
        output_stream<<"Triangle "<<selected_triangle<<" ("<<node1<<", "<<node2<<", "<<node3<<")"<<std::endl;
        output_stream<<std::endl;
        output_stream<<"Vertex "<<node1<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(triangulated_area.particles.X(node1))<<std::endl;}
        triangulated_area.particles.Print(output_stream,node1);
        output_stream<<std::endl;
        output_stream<<"Vertex "<<node2<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(triangulated_area.particles.X(node2))<<std::endl;}
        triangulated_area.particles.Print(output_stream,node2);
        output_stream<<std::endl;
        output_stream<<"Vertex "<<node3<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(triangulated_area.particles.X(node3))<<std::endl;}
        triangulated_area.particles.Print(output_stream,node3);}
}
//#####################################################################
// Function Draw_Vertices
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA<T>::
Draw_Vertices() const
{
    if(draw_vertices){
        vertex_color.Send_To_GL_Pipeline();
        OpenGL_Begin(GL_POINTS);
        if(!triangulated_area.mesh.neighbor_nodes) triangulated_area.mesh.Initialize_Neighbor_Nodes();
        for(int i=0;i<triangulated_area.particles.Size();i++)
            if((*triangulated_area.mesh.neighbor_nodes)(i).m)
                OpenGL_Vertex(triangulated_area.particles.X(i));
        OpenGL_End();}
}
//#####################################################################
// Function Draw_Vertices_For_Selection
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA<T>::
Draw_Vertices_For_Selection() const
{
    OPENGL_SELECTION::Draw_Vertices_For_Selection(triangulated_area.mesh,triangulated_area.particles);
}
//#####################################################################
// Function Draw_Segments
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA<T>::
Draw_Segments() const
{
    PHYSBAM_ASSERT(triangulated_area.mesh.segment_mesh);
    segment_color.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_LINES);
    for(int i=0;i<triangulated_area.mesh.segment_mesh->elements.m;i++){
        int node1,node2;triangulated_area.mesh.segment_mesh->elements(i).Get(node1,node2);
        OpenGL_Line(triangulated_area.particles.X(node1),triangulated_area.particles.X(node2));}
    OpenGL_End();
}
//#####################################################################
// Function Draw_Segments_For_Selection
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA<T>::
Draw_Segments_For_Selection() const
{
    PHYSBAM_ASSERT(triangulated_area.mesh.segment_mesh);
    glPushAttrib(GL_LINE_BIT);
    glLineWidth(OPENGL_PREFERENCES::selection_line_width);
    glPushName(0);
    for(int i=0;i<triangulated_area.mesh.segment_mesh->elements.m;i++){
        int node1,node2;
        triangulated_area.mesh.segment_mesh->elements(i).Get(node1,node2);
        glLoadName(i);
        OpenGL_Begin(GL_LINES);
        OpenGL_Line(triangulated_area.particles.X(node1),triangulated_area.particles.X(node2));
        OpenGL_End();}
    glPopName();
    glPopAttrib();
}
//#####################################################################
// Function Draw_Triangles
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA<T>::
Draw_Triangles(const bool use_color_map) const
{
    OpenGL_Begin(GL_TRIANGLES);
    if(use_color_map) triangle_color.Send_To_GL_Pipeline();
    for(int i=0;i<triangulated_area.mesh.elements.m;i++){
        if(color_map && use_color_map) (*color_map)(i).Send_To_GL_Pipeline();
        int node1,node2,node3;triangulated_area.mesh.elements(i).Get(node1,node2,node3);
        OpenGL_Triangle(triangulated_area.particles.X(node1),triangulated_area.particles.X(node2),triangulated_area.particles.X(node3));}
    OpenGL_End();

    OpenGL_Begin(GL_TRIANGLES);
    if(use_color_map) triangle_inverted_color.Send_To_GL_Pipeline();
    for(int i=0;i<triangulated_area.mesh.elements.m;i++){
        if(color_map && use_color_map) (*color_map)(i).Send_To_GL_Pipeline();
        int node1,node2,node3;triangulated_area.mesh.elements(i).Get(node1,node2,node3);
        OpenGL_Triangle(triangulated_area.particles.X(node2),triangulated_area.particles.X(node1),triangulated_area.particles.X(node3));}
    OpenGL_End();
}
//#####################################################################
// Function Draw_Triangles_For_Selection
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA<T>::
Draw_Triangles_For_Selection() const
{
    glPushName(0);
    for(int i=0;i<triangulated_area.mesh.elements.m;i++){
        int node1,node2,node3;triangulated_area.mesh.elements(i).Get(node1,node2,node3);
        glLoadName(i);
        OpenGL_Begin(GL_TRIANGLES);
        OpenGL_Triangle(triangulated_area.particles.X(node1),triangulated_area.particles.X(node2),triangulated_area.particles.X(node3));
        OpenGL_Triangle(triangulated_area.particles.X(node2),triangulated_area.particles.X(node1),triangulated_area.particles.X(node3));
        OpenGL_End();}
    glPopName();
}
//#####################################################################
// Selection object functions
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_TRIANGULATED_AREA<T>::
Selection_Bounding_Box() const
{
    if(selected_vertex>=0) return World_Space_Box(RANGE<TV>(triangulated_area.particles.X(selected_vertex)));
    if(selected_segment>=0) return World_Space_Box(RANGE<TV>::Bounding_Box(triangulated_area.particles.X.Subset(triangulated_area.mesh.segment_mesh->elements(selected_segment))));
    if(selected_triangle>=0) return World_Space_Box(RANGE<TV>::Bounding_Box(triangulated_area.particles.X.Subset(triangulated_area.mesh.elements(selected_triangle))));
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_TRIANGULATED_AREA<float>;
template class OPENGL_TRIANGULATED_AREA<double>;
}

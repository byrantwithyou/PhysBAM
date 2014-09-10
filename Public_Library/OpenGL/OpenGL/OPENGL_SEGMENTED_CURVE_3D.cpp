//#####################################################################
// Copyright 2002-2008, Kevin Der, Ronald Fedkiw, Eran Guendelman, Sergey Koltakov, Michael Lentine, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Data_Structures/QUEUE.h>
#include <Tools/Math_Tools/sign.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_3D.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
using namespace PhysBAM;
//#####################################################################
// Function Camera_Oriented_Normal
//#####################################################################
template<class TV>
TV Camera_Oriented_Normal(const TV& camera_direction,const VECTOR<TV,2>& Xs)
{
    TV tangent=(Xs[0]-Xs[1]).Normalized();
    TV normal=camera_direction.Projected_Orthogonal_To_Unit_Direction(tangent).Normalized();
    if(TV::Dot_Product(normal,camera_direction)<0) normal=-normal;
    return normal;
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_SEGMENTED_CURVE_3D<T>::
Display() const
{   
    if(use_solid_color){glPushAttrib(GL_LIGHTING_BIT);glDisable(GL_LIGHTING);color.Send_To_GL_Pipeline();}
    else OPENGL_MATERIAL::Plastic(color).Send_To_GL_Pipeline();

    int edge,node1,node2,len;
    ARRAY<int> selected_edges=Get_Selected_Edges();
    len=selected_edges.m/2+2;
    HASHTABLE<int,bool> seen(len);
    if(hide_unselected){
        for(int i=0;i<selected_edges.m;i++){edge=selected_edges(i);parent_curve->curve.mesh.elements(edge).Get(node1,node2);seen.Set(node1,true);seen.Set(node2,true);}
        for(int i=0;i<curve.mesh.elements.m;i++){const VECTOR<int,2> edge=curve.mesh.elements(i);
            if(seen.Contains(edge[0])||seen.Contains(edge[1])){seen.Set(edge[0],true);seen.Set(edge[1],true);}}}
    
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    GLint mode=0;
    glGetIntegerv(GL_RENDER_MODE, &mode);
    if(mode==GL_SELECT){
        glPushName(1);
        Draw_Vertices_For_Selection();
        glLoadName(2);
        Draw_Segments_For_Selection();
        glPopName();}
    else
    {
        if(use_solid_color){
            OpenGL_Begin(GL_LINES);
            for(int t=0;t<curve.mesh.elements.m;t++){
                int i=curve.mesh.elements(t)(0),j=curve.mesh.elements(t)(1);
                if(!hide_unselected || (seen.Contains(i) && seen.Contains(j))){
                    OpenGL_Vertex(curve.particles.X(i));OpenGL_Vertex(curve.particles.X(j));}}
            OpenGL_End();}
        else{
            if(smooth_normals){
                Initialize_Vertex_Normals();
                OpenGL_Begin(GL_LINES);
                for(int i=0;i<segment_nodes.m;i++){int p=segment_nodes(i);
                    const ARRAY<int>& incident=(*curve.mesh.incident_elements)(p);
                    for(int j=0;j<incident.m;j++){const VECTOR<int,2>& nodes=curve.mesh.elements(incident(j));
                        if(!hide_unselected || (seen.Contains(nodes[0]) && seen.Contains(nodes[1]))){
                            OpenGL_Normal(vertex_normals.Get(nodes[0]));OpenGL_Vertex(curve.particles.X(nodes[0]));
                            OpenGL_Normal(vertex_normals.Get(nodes[1]));OpenGL_Vertex(curve.particles.X(nodes[1]));}}}
                OpenGL_End();}
            else{
                OpenGL_Begin(GL_LINES);
                TV camera_direction=OPENGL_WORLD<T>::Singleton()->Get_Camera_Position()-OPENGL_WORLD<T>::Singleton()->Get_Target_Position().Normalized();
                for(int s=0;s<curve.mesh.elements.m;s++){
                    const VECTOR<int,2>& nodes=curve.mesh.elements(s);
                    if(!hide_unselected || (seen.Contains(nodes[0]) && seen.Contains(nodes[1]))){
                        VECTOR<TV,2> Xs(curve.particles.X.Subset(nodes));
                        for(int line_vertices=0;line_vertices<2;line_vertices++) OpenGL_Normal(Camera_Oriented_Normal(camera_direction,Xs));
                        OpenGL_Vertex(Xs[0]);OpenGL_Vertex(Xs[1]);}}
                OpenGL_End();}}

        if(current_selection){
            if(current_selection->type==OPENGL_SELECTION<T>::SEGMENTED_CURVE_VERTEX_3D){
                int index=((OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_3D<T> *)current_selection)->index;
                OPENGL_SELECTION<T>::Draw_Highlighted_Vertex(curve.particles.X(index),index);} 
            else if(current_selection->type==OPENGL_SELECTION<T>::SEGMENTED_CURVE_SEGMENT_3D){
                int index=((OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_3D<T> *)current_selection)->index;
                int node1,node2;curve.mesh.elements(index).Get(node1,node2);
                OPENGL_SELECTION<T>::Draw_Highlighted_Segment(curve.particles.X(node1),curve.particles.X(node2),index);} 
            else if(current_selection->type==OPENGL_SELECTION<T>::SEGMENTED_CURVE_3D) {
                int node1,node2;
                ARRAY<VECTOR<TV,2> > lines;
                for(int i=0;i<selected_edges.m;i++){
                    curve.mesh.elements(selected_edges(i)).Get(node1,node2);
                    lines.Append(VECTOR<TV,2>(curve.particles.X(node1),curve.particles.X(node2)));}
                OPENGL_SELECTION<T>::Draw_Highlighted_Curve(lines);}}}

    if(draw_vertices){
        vertex_color.Send_To_GL_Pipeline();
        glPointSize(5.0f);
        OpenGL_Begin(GL_POINTS);
        for(int t=0;t<curve.particles.Size();t++) OpenGL_Vertex(curve.particles.X(t));
        OpenGL_End();}

    if(use_solid_color) glPopAttrib();

    glPopMatrix();
}
//#####################################################################
// Function Initialize_Vertex_Normals
//#####################################################################
template<class T> void OPENGL_SEGMENTED_CURVE_3D<T>::
Initialize_Vertex_Normals() const
{
    TV camera_direction=TV(OPENGL_WORLD<T>::Singleton()->Get_Camera_Position()-OPENGL_WORLD<T>::Singleton()->Get_Target_Position()).Normalized();
    bool incident_segments_defined=(curve.mesh.incident_elements!=0);if(!incident_segments_defined) curve.mesh.Initialize_Incident_Elements();
    segment_nodes.Remove_All();curve.mesh.elements.Flattened().Get_Unique(segment_nodes);
    vertex_normals.Remove_All();
    for(int i=0;i<segment_nodes.m;i++){int p=segment_nodes(i);
        const ARRAY<int>& incident=(*curve.mesh.incident_elements)(p);
        TV normal;
        for(int j=0;j<incident.m;j++){const VECTOR<int,2>& nodes=curve.mesh.elements(incident(j));VECTOR<TV,2> Xs(curve.particles.X.Subset(nodes));
            normal+=Camera_Oriented_Normal(camera_direction,Xs);}
        normal.Normalize();
        vertex_normals.Set(p,normal);}
}
//#####################################################################
// Function Turn_Smooth_Shading_On
//#####################################################################
template<class T> void OPENGL_SEGMENTED_CURVE_3D<T>::
Turn_Smooth_Shading_On()
{
    Initialize_Vertex_Normals();
    smooth_normals=true;
}
//#####################################################################
// Function Turn_Smooth_Shading_Off
//#####################################################################
template<class T> void OPENGL_SEGMENTED_CURVE_3D<T>::
Turn_Smooth_Shading_Off()
{
    smooth_normals=false;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_SEGMENTED_CURVE_3D<T>::
Bounding_Box() const
{
    return World_Space_Box(RANGE<TV>::Bounding_Box(curve.particles.X));
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_SEGMENTED_CURVE_3D<T>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    OPENGL_SELECTION<T>* selection = 0;
    if(buffer_size==2)
    {
        if(buffer[0]==1) selection = Get_Vertex_Selection(buffer[1]);
        else if(buffer[0]==2) selection = Get_Segment_Selection(buffer[1]);
    }
    return selection;
}
//#####################################################################
// Function Get_Selected_Edges
//#####################################################################
template<class T> ARRAY<int> OPENGL_SEGMENTED_CURVE_3D<T>::
Get_Selected_Edges() const
{
    ARRAY<int> sel_edges;
    if(!(parent_curve && parent_curve->current_selection)) return sel_edges;
    if(parent_curve->current_selection->type==OPENGL_SELECTION<T>::SEGMENTED_CURVE_SEGMENT_3D){
        int index=((OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_3D<T> *)parent_curve->current_selection)->index;
        sel_edges.Append(index);}
    else if(parent_curve->current_selection->type==OPENGL_SELECTION<T>::SEGMENTED_CURVE_3D) {
        int edge;
        int index=((OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_3D<T> *)parent_curve->current_selection)->index;
        ARRAY<int> adjacent;
        parent_curve->curve.mesh.Initialize_Adjacent_Elements();
        QUEUE<int> edges(parent_curve->curve.mesh.adjacent_elements->m);
        HASHTABLE<int,bool> seen(parent_curve->curve.mesh.adjacent_elements->m);
        edges.Enqueue(index);
        while(!edges.Empty()){
            edge=edges.Dequeue();
            sel_edges.Append(edge);
            adjacent=(*parent_curve->curve.mesh.adjacent_elements)(edge);
            for(int i=0;i<adjacent.m;i++) if(!seen.Contains(adjacent(i))){seen.Set(adjacent(i),true);edges.Enqueue(adjacent(i));}}}
    return sel_edges;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_SEGMENTED_CURVE_3D<T>::
Highlight_Selection(OPENGL_SELECTION<T>* selection)
{
    delete current_selection;current_selection = 0;
    // Make a copy of selection
    if(selection->type==OPENGL_SELECTION<T>::SEGMENTED_CURVE_VERTEX_3D) 
        current_selection=new OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_3D<T>(this,((OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_3D<T>*)selection)->index);
    else if(selection->type==OPENGL_SELECTION<T>::SEGMENTED_CURVE_SEGMENT_3D)
        current_selection=new OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_3D<T>(this,((OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_3D<T>*)selection)->index);
    else if(selection->type==OPENGL_SELECTION<T>::SEGMENTED_CURVE_3D)
        current_selection=new OPENGL_SELECTION_SEGMENTED_CURVE_3D<T>(this,((OPENGL_SELECTION_SEGMENTED_CURVE_3D<T>*)selection)->index);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_SEGMENTED_CURVE_3D<T>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION<T>* selection) const
{
    if(selection->type==OPENGL_SELECTION<T>::SEGMENTED_CURVE_VERTEX_3D){
        int index=((OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_3D<T> *)selection)->index;
        output_stream<<"Vertex "<<index<<std::endl;
        curve.particles.Print(output_stream,index);}
    else if(selection->type==OPENGL_SELECTION<T>::SEGMENTED_CURVE_SEGMENT_3D){
        int index=((OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_3D<T> *)selection)->index;
        int node1,node2;curve.mesh.elements(index).Get(node1,node2);
        output_stream<<"Segment "<<index<<" ("<<node1<<", "<<node2<<")"<<std::endl;
        output_stream<<std::endl;
        output_stream<<"Vertex "<<node1<<std::endl;
        curve.particles.Print(output_stream,node1);
        output_stream<<std::endl;
        output_stream<<"Vertex "<<node2<<std::endl;
        curve.particles.Print(output_stream,node2);}
    else if(selection->type==OPENGL_SELECTION<T>::SEGMENTED_CURVE_3D){
        int node1,node2;
        int index=((OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_3D<T> *)selection)->index;
        ARRAY<int> edges=Get_Selected_Edges(); 
        output_stream<<"Curve "<<index<<std::endl;
        output_stream<<std::endl;
        for(int i=0;i<edges.m;i++) {curve.mesh.elements(edges(i)).Get(node1,node2);output_stream<<"("<<node1<<", "<<node2<<")"<<std::endl;}}
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_SEGMENTED_CURVE_3D<T>::
Clear_Highlight()
{
    delete current_selection;current_selection = 0;
}
//#####################################################################
// Function Get_Vertex_Selection
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_SEGMENTED_CURVE_3D<T>::
Get_Vertex_Selection(int index)
{
    return new OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_3D<T>(this, index);
}
//#####################################################################
// Function Get_Segment_Selection
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_SEGMENTED_CURVE_3D<T>::
Get_Segment_Selection(int index)
{
    return new OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_3D<T>(this, index);
}
//#####################################################################
// Function Get_Curve_Selection
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_SEGMENTED_CURVE_3D<T>::
Get_Curve_Selection(int index)
{
    return new OPENGL_SELECTION_SEGMENTED_CURVE_3D<T>(this, index);
}
//#####################################################################
// Function Draw_Vertices_For_Selection
//#####################################################################
template<class T>
void OPENGL_SEGMENTED_CURVE_3D<T>::
Draw_Vertices_For_Selection() const
{
    OPENGL_SELECTION<T>::Draw_Vertices_For_Selection(curve.mesh,curve.particles);
}
//#####################################################################
// Function Draw_Segments_For_Selection
//#####################################################################
template<class T>
void OPENGL_SEGMENTED_CURVE_3D<T>::
Draw_Segments_For_Selection() const
{
    glPushAttrib(GL_LINE_BIT);
    glLineWidth(OPENGL_PREFERENCES::selection_line_width);
    glPushName(0);
    for(int i=0;i<curve.mesh.elements.m;i++){
        int node1,node2;curve.mesh.elements(i).Get(node1,node2);
        glLoadName(i);
        OpenGL_Begin(GL_LINES);
        OpenGL_Line(curve.particles.X(node1),curve.particles.X(node2));
        OpenGL_End();}
    glPopName();
    glPopAttrib();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T>
RANGE<VECTOR<T,3> > OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_3D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const SEGMENTED_CURVE<VECTOR<T,3> > &curve=((OPENGL_SEGMENTED_CURVE_3D<T> *)object)->curve;
    RANGE<VECTOR<T,3> > box(VECTOR<T,3>(curve.particles.X(index)));
    return object->World_Space_Box(box);
}

template<class T>
RANGE<VECTOR<T,3> > OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_3D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const SEGMENTED_CURVE<VECTOR<T,3> > &curve=((OPENGL_SEGMENTED_CURVE_3D<T> *)object)->curve;
    return object->World_Space_Box(RANGE<VECTOR<T,3> >::Bounding_Box(curve.particles.X.Subset(curve.mesh.elements(index))));
}

template<class T>
RANGE<VECTOR<T,3> > OPENGL_SELECTION_SEGMENTED_CURVE_3D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const OPENGL_SEGMENTED_CURVE_3D<T> *curve=(OPENGL_SEGMENTED_CURVE_3D<T> *)object;
    return curve->Bounding_Box();
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_SEGMENTED_CURVE_3D<float>;
template class OPENGL_SEGMENTED_CURVE_3D<double>;
}

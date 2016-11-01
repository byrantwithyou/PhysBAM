//#####################################################################
// Copyright 2002-2005, Eran Guendelman, Eilene Hao, Neil Molino, Robert Bridson, Igor Neverov, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Matrices/MATRIX_4X4.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Geometry/Topology_Based_Geometry/B_SPLINE_PATCH.h>
#include <OpenGL/OpenGL/OPENGL_B_SPLINE_PATCH.h>
#include <OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor from geometry, two-sided
//#####################################################################
template<class T> OPENGL_B_SPLINE_PATCH<T>::
OPENGL_B_SPLINE_PATCH(STREAM_TYPE stream_type,B_SPLINE_PATCH<VECTOR<T,3>>& patch_input,const OPENGL_MATERIAL& front_material_input,const OPENGL_MATERIAL& back_material_input)
    :OPENGL_OBJECT<T>(stream_type),patch(patch_input),two_sided(true),front_material(front_material_input),
    back_material(back_material_input),use_display_list(false),
    owns_display_list(false),selected_vertex(-1),selected_element(-1),
    current_node(1),highlight_current_node(false),wireframe_only(false), draw_particles(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_B_SPLINE_PATCH<T>::
~OPENGL_B_SPLINE_PATCH()
{
    if(owns_display_list) glDeleteLists(display_list_id,1);
}
//#####################################################################
// Function Highlight_Current_Node
//#####################################################################
template<class T> void OPENGL_B_SPLINE_PATCH<T>::
Highlight_Current_Node() const
{
    int node=current_node%patch.particles.Size();
    OPENGL_SHAPES::Draw_Dot(patch.particles.X(node),OPENGL_COLOR(1,0,1),7);
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_B_SPLINE_PATCH<T>::
Display() const
{
    if(two_sided){
        glDisable(GL_CULL_FACE);
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
        front_material.Send_To_GL_Pipeline(GL_FRONT);
        back_material.Send_To_GL_Pipeline(GL_BACK);}
    else front_material.Send_To_GL_Pipeline();

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    if(draw_particles) for(int i=0;i<patch.particles.Size();i++) OPENGL_SHAPES::Draw_Dot(patch.particles.X(i),OPENGL_COLOR(1,0,1),7);

    GLint mode=0;
    glGetIntegerv(GL_RENDER_MODE,&mode);

    if(wireframe_only){
        glPushAttrib(GL_POLYGON_BIT);
        glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);}

    if(mode == GL_SELECT){
        glPushName(0);
        Draw_Vertices_For_Selection();
        glLoadName(1);
        Draw_Elements_For_Selection();
        glPopName();}
    else if(use_display_list) glCallList(display_list_id);
    else Draw();
    if(wireframe_only) glPopAttrib();

    if(mode != GL_SELECT){
        if(selected_vertex>=0)
            OPENGL_SELECTION::Draw_Highlighted_Vertex(patch.particles.X(selected_vertex),selected_vertex);
        if(selected_element>=0) Draw_Selected_Element(selected_element);
        if(highlight_current_node) Highlight_Current_Node();}

    if(two_sided){
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,0);
        glEnable(GL_CULL_FACE);}

    glPopMatrix();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_B_SPLINE_PATCH<T>::
Bounding_Box() const
{
    return World_Space_Box(RANGE<VECTOR<T,3> >::Bounding_Box(patch.particles.X));
}
//#####################################################################
// Function Get_Selection_Priority
//#####################################################################
template<class T> int OPENGL_B_SPLINE_PATCH<T>::
Get_Selection_Priority(ARRAY_VIEW<GLuint> indices)
{
    PHYSBAM_ASSERT(indices.m==2);
    const static int priority[]={74,73};
    return priority[indices(0)];
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> bool OPENGL_B_SPLINE_PATCH<T>::
Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers)
{
    if(indices(0) == 0) selected_vertex=indices(1);
    else if(indices(0) == 1) selected_element=indices(1);
    else return false;
    return true;
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_B_SPLINE_PATCH<T>::
Clear_Selection()
{
    selected_vertex=-1;
    selected_element=-1;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_B_SPLINE_PATCH<T>::
Print_Selection_Info(std::ostream& output_stream) const
{
    Print_Selection_Info(output_stream,0);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_B_SPLINE_PATCH<T>::
Print_Selection_Info(std::ostream& output_stream,MATRIX<T,4>* transform) const
{
    if(selected_vertex>=0){
        output_stream<<"Vertex "<<selected_vertex<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(patch.particles.X(selected_vertex))<<std::endl;}
        patch.particles.Print(output_stream,selected_vertex);}
    else if(selected_element>=0){
        output_stream<<"Patch ("<<selected_element<<")"<<std::endl;
        output_stream<<std::endl;
        VECTOR<int,16> a=patch.Control_Points_For_Element(selected_element);
        for(int i=0;i<a.m;i++){
            output_stream<<"Control Point "<<a(i)<<std::endl;
            if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(patch.particles.X(a(i)))<<std::endl;}
            patch.particles.Print(output_stream,a(i));}}
}
//#####################################################################
// Function Set_Front_Material
//#####################################################################
template<class T> void OPENGL_B_SPLINE_PATCH<T>::
Set_Front_Material(const OPENGL_MATERIAL& material_input)
{
    front_material=material_input;
}
//#####################################################################
// Function Set_Back_Material
//#####################################################################
template<class T> void OPENGL_B_SPLINE_PATCH<T>::
Set_Back_Material(const OPENGL_MATERIAL& material_input)
{
    back_material=material_input;
}
//#####################################################################
// Function
//#####################################################################
template<class T> int OPENGL_B_SPLINE_PATCH<T>::
Create_Display_List()
{
    PHYSBAM_ASSERT(!owns_display_list);
    owns_display_list=true;
    use_display_list=true;
    display_list_id=glGenLists(1);
    Reinitialize_Display_List();
    return display_list_id;
}
//#####################################################################
// Function
//#####################################################################
template<class T> void OPENGL_B_SPLINE_PATCH<T>::
Use_Display_List(int input_display_list_id)
{
    PHYSBAM_ASSERT(!owns_display_list);
    use_display_list=true;
    display_list_id=input_display_list_id;
}
//#####################################################################
// Function Reinitialize_Display_List
//#####################################################################
template<class T> void OPENGL_B_SPLINE_PATCH<T>::
Reinitialize_Display_List()
{
    glNewList(display_list_id,GL_COMPILE);
    Draw();
    glEndList();
}
//#####################################################################
// Function Draw
//#####################################################################
template<class T> void OPENGL_B_SPLINE_PATCH<T>::
Draw() const
{
    for(int i=0;i<patch.knots_s.m-1;i++){
        T s0=patch.knots_s(i);
        T s1=patch.knots_s(i+1);
        for(int ii=0;ii<substeps;ii++){
            OpenGL_Begin(GL_TRIANGLE_STRIP);
            for(int j=0;j<patch.knots_t.m-1;j++){
                T t0=patch.knots_t(j);
                T t1=patch.knots_t(j+1);
                for(int jj=0;jj<=substeps;jj++){
                    for(int k=0;k<2;k++){
                        VECTOR<VECTOR<T,3>,2> tangents;
                        T s = s0 + (ii+k)*(s1-s0)/substeps;
                        T t = t0 + jj*(t1-t0)/substeps;
                        VECTOR<T,3> x=patch.Evaluate(s,t,&tangents);
                        VECTOR<T,3> normal=tangents(0).Cross(tangents(1)).Normalized();
                        OpenGL_Normal(normal);OpenGL_Vertex(x);}}}
                OpenGL_End();}}
}

//#####################################################################
// Function Draw_Selected_Element
//#####################################################################
template<class T> void OPENGL_B_SPLINE_PATCH<T>::
Draw_Selected_Element(int index) const
{
    VECTOR<int,16> a=patch.Control_Points_For_Element(index);
    for(int i=0;i<a.m;i++)
            OPENGL_SELECTION::Draw_Highlighted_Vertex(patch.particles.X(a(i)),a(i));
    RANGE<VECTOR<T,2>> r=patch.Range_For_Element(index);
    glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    glLineWidth(OPENGL_PREFERENCES::highlighted_line_width);
    const OPENGL_COLOR& color=OPENGL_PREFERENCES::selection_highlight_color;
    color.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_LINE_LOOP);
    T s0=r.min_corner(0);
    T t0=r.min_corner(1);
    T s=s0;
    T t=t0;
    VECTOR<T,2> dr=r.Edge_Lengths();
    for(int ii=0;ii<=substeps;ii++){
        s = s0 + ii*dr(0)/substeps;
        VECTOR<T,3> x=patch.Evaluate(s,t);
        OpenGL_Vertex(x);}
    for(int jj=0;jj<=substeps;jj++){
        t = t0 + jj*dr(1)/substeps;
        VECTOR<T,3> x=patch.Evaluate(s,t);
        OpenGL_Vertex(x);}
    for(int ii=substeps;ii>=0;ii--){
        s = s0 + ii*dr(0)/substeps;
        VECTOR<T,3> x=patch.Evaluate(s,t);
        OpenGL_Vertex(x);}
    for(int jj=substeps;jj>=0;jj--){
        t = t0 + jj*dr(1)/substeps;
        VECTOR<T,3> x=patch.Evaluate(s,t);
        OpenGL_Vertex(x);}
    OpenGL_End();
    glPopAttrib();
}
//#####################################################################
// Function Draw_Vertices_For_Selection
//#####################################################################
template<class T> void OPENGL_B_SPLINE_PATCH<T>::
Draw_Vertices_For_Selection() const
{
    glPushAttrib(GL_POINT_BIT);
    glPointSize(OPENGL_PREFERENCES::selection_point_size);
    glPushName(0);
    for(int i=0;i<patch.control_points.domain.max_corner(0)-1;i++){
        for(int j=0;j<patch.control_points.domain.max_corner(1)-1;j++){
            int p=patch.control_points(i,j);
            glLoadName(p);
            OpenGL_Begin(GL_POINTS);
            OpenGL_Vertex(patch.particles.X(p));
            OpenGL_End();}}
    glPopName();
    glPopAttrib();
}
//#####################################################################
// Function Draw_Patches_For_Selection
//#####################################################################
template<class T> void OPENGL_B_SPLINE_PATCH<T>::
Draw_Elements_For_Selection() const
{
    glPushName(0);
    for(int e=0;e<patch.m;e++){
        RANGE<VECTOR<T,2>> r=patch.Range_For_Element(e);
        VECTOR<T,2> dr=r.Edge_Lengths();
        glLoadName(e);
        OpenGL_Begin(GL_TRIANGLE_STRIP);
            for(int ii=0;ii<substeps;ii++){
                for(int jj=0;jj<=substeps;jj++){
                    for(int k=0;k<2;k++){
                        T s=r.min_corner(0)+(ii+k)*dr(0)/substeps;
                        T t=r.min_corner(1)+jj*dr(1)/substeps;
                        VECTOR<T,3> x=patch.Evaluate(s,t);
                        OpenGL_Vertex(x);}}}
            OpenGL_End();}
    glPopName();
}
//#####################################################################
//
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_B_SPLINE_PATCH<T>::
Selection_Bounding_Box() const
{
    if(selected_vertex>=0) return World_Space_Box(RANGE<TV>(patch.particles.X(selected_vertex)));
    if(selected_element>=0) return World_Space_Box(RANGE<TV>::Bounding_Box(patch.particles.X.Subset(patch.Control_Points_For_Element(selected_element))));
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_B_SPLINE_PATCH<float>;
template class OPENGL_B_SPLINE_PATCH<double>;
}

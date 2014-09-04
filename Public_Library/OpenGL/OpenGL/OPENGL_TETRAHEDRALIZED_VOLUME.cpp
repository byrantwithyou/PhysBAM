//#####################################################################
// Copyright 2005-2009, Zhaosheng Bao, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Sergey Koltakov, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Matrices/MATRIX_4X4.h>
#include <Geometry/Basic_Geometry/PLANE.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Topology/TRIANGLE_MESH.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <OpenGL/OpenGL/OPENGL_TETRAHEDRALIZED_VOLUME.h>
using namespace PhysBAM;
namespace{
static VECTOR<int,4> Spring_Nodes(unsigned char pair_id,const VECTOR<int,4>& n)
{
    switch(pair_id){
        case 0: return VECTOR<int,4>(n[0],n[1],n[3],n[2]); // point face
        case 1: return VECTOR<int,4>(n[1],n[0],n[2],n[3]); // point face
        case 2: return VECTOR<int,4>(n[2],n[0],n[3],n[1]); // point face
        case 3: return VECTOR<int,4>(n[3],n[0],n[1],n[2]); // point face
        case 4: return VECTOR<int,4>(n[0],n[1],n[2],n[3]); // edge edge
        case 5: return VECTOR<int,4>(n[1],n[2],n[0],n[3]); // edge edge
        case 6: return VECTOR<int,4>(n[0],n[2],n[3],n[1]); // edge edge
        default: PHYSBAM_FATAL_ERROR();}
}
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_TETRAHEDRALIZED_VOLUME<T>::
OPENGL_TETRAHEDRALIZED_VOLUME(const OPENGL_MATERIAL& material_input,const OPENGL_MATERIAL& inverted_material_input)
    :material(material_input),inverted_material(inverted_material_input),use_inverted_material(true),mesh(0),particles(0),current_tetrahedron(0),current_node(0),
    current_boundary_triangle(0),boundary_only(true),draw_subsets(false),cutaway_mode(0),cutaway_fraction((T).5),color_map(0),smooth_normals(false),
    vertex_normals(0),current_selection(0)
{
    Initialize();
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_TETRAHEDRALIZED_VOLUME<T>::
OPENGL_TETRAHEDRALIZED_VOLUME(TETRAHEDRON_MESH* mesh_input,GEOMETRY_PARTICLES<VECTOR<T,3> >* particles_input,const OPENGL_MATERIAL& material_input,
    const OPENGL_MATERIAL& inverted_material_input,bool initialize,ARRAY<OPENGL_COLOR>* color_map_input)
    :material(material_input),inverted_material(inverted_material_input),use_inverted_material(true),mesh(mesh_input),particles(particles_input),current_tetrahedron(0),current_node(0),
    current_boundary_triangle(0),boundary_only(true),draw_subsets(false),cutaway_mode(0),cutaway_fraction((T).5),color_map(color_map_input),smooth_normals(false),
    vertex_normals(0),current_selection(0)
{
    if(initialize)Initialize();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_TETRAHEDRALIZED_VOLUME<T>::
~OPENGL_TETRAHEDRALIZED_VOLUME()
{
    delete vertex_normals;
}
//#####################################################################
// Function Display_In_Color
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Display() const
{
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    GLint mode=0;
#ifndef USE_OPENGLES
    glGetIntegerv(GL_RENDER_MODE,&mode);
#endif
    glDisable(GL_CULL_FACE);
#ifndef USE_OPENGLES
    if(mode==GL_SELECT){
        glPushName(1);
        Draw_Vertices_For_Selection();
        glLoadName(2);
        Draw_Tetrahedra_For_Selection();
        glPopName();}
    else if(cutaway_mode){
#else
    if(cutaway_mode){
#endif
        if(boundary_only) Draw_Boundary_Triangles(cutaway_mesh);
        else Draw_Wireframe_Mesh(cutaway_mesh);}
    else{
        if(spectrum.m) Draw_In_Color_From_Spectrum();
        if(color_map) Draw_In_Color_From_Color_Map();
        if(boundary_only) Draw_Boundary_Triangles(*mesh);
        else Draw_Wireframe_Mesh(*mesh);
        if(draw_subsets){Draw_Subset();Draw_Subset_Particles();Draw_Subset_Triangles();}
        //Draw_Current_Tetrahedron();
        //Highlight_Boundary_Nodes_Of_Current_Tetrahedron();
        //Highlight_Boundary_Normal_Vectors_Of_Current_Tetrahedron();
        //Highlight_Current_Node();
        //Highlight_Nodes_Of_Minimum_Valence();
        //Highlight_Current_Boundary_Triangle();

        if(current_selection){
            if(current_selection->type==OPENGL_SELECTION<T>::TETRAHEDRALIZED_VOLUME_VERTEX){
                int index=((OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_VERTEX<T>*)current_selection)->index;
                OPENGL_SELECTION<T>::Draw_Highlighted_Vertex(particles->X(index),index);}
            else if(current_selection->type==OPENGL_SELECTION<T>::TETRAHEDRALIZED_VOLUME_TETRAHEDRON){
                int index=((OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_TETRAHEDRON<T>*)current_selection)->index;
                const VECTOR<int,4>& element_nodes=mesh->elements(index);
                ARRAY_VIEW<const TV> X(particles->X);
                OPENGL_SELECTION<T>::Draw_Highlighted_Tetrahedron_Boundary(X(element_nodes[0]),X(element_nodes[1]),X(element_nodes[2]),X(element_nodes[3]),index);
                T distance;TV min_normal,weights;
                int spring=Find_Shortest_Spring(element_nodes,distance,min_normal,weights);
                VECTOR<int,4> spring_nodes;
                spring_nodes=Spring_Nodes(spring,element_nodes);

                glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
                glDisable(GL_LIGHTING);
                glLineWidth(OPENGL_PREFERENCES::highlighted_line_width*2);
                OPENGL_COLOR colors[]={OPENGL_COLOR::Red(),OPENGL_COLOR::Blue(),OPENGL_COLOR::Green(),OPENGL_COLOR::Yellow(),OPENGL_COLOR::Cyan(),OPENGL_COLOR::Magenta(),
                                       OPENGL_COLOR::White()};
                if(spring>=0){
                    colors[spring].Send_To_GL_Pipeline();
                    ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
                    if(spring<4) OpenGL_Line(X(spring_nodes[0]),TRIANGLE_3D<T>(X.Subset(spring_nodes.Remove_Index(0))).Point_From_Barycentric_Coordinates(weights),vertices);
                    else if(spring<7) OpenGL_Line((1-weights.x)*X(spring_nodes[0])+weights.x*X(spring_nodes[1]),(1-weights.y)*X(spring_nodes[2])+weights.y*X(spring_nodes[3]),vertices);
                    OpenGL_Draw_Arrays(GL_LINES,3,vertices);
                    glPopAttrib();}}}}
    glPopMatrix();

    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,0);
    glEnable(GL_CULL_FACE);
}
//#####################################################################
// Function Find_Shortest_Spring
//#####################################################################
template<class T> int OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Find_Shortest_Spring(const VECTOR<int,4>& element_nodes,T& minimum_signed_distance,TV& minimum_normal,TV& weights) const
{
    ARRAY_VIEW<const TV> X(particles->X);
    // Find Shortest Altitude
    int maximum_cross_squared_index=-1;T maximum_cross_squared=(T)-FLT_MAX;TV maximum_cross;
    for(unsigned char h=0;h<4;h++){
        VECTOR<int,4> spring_nodes=Spring_Nodes(h,element_nodes);
        TV u_cross_v=TV::Cross_Product(X(spring_nodes[2])-X(spring_nodes[1]),X(spring_nodes[3])-X(spring_nodes[1]));
        T u_cross_v_squared=u_cross_v.Magnitude_Squared();
        if(u_cross_v_squared>maximum_cross_squared){maximum_cross_squared_index=h;maximum_cross_squared=u_cross_v_squared;maximum_cross=u_cross_v;}}
    for(unsigned char h=4;h<7;h++){
        VECTOR<int,4> spring_nodes=Spring_Nodes(h,element_nodes);
        TV u_cross_v=TV::Cross_Product(X(spring_nodes[1])-X(spring_nodes[0]),X(spring_nodes[3])-X(spring_nodes[2]));
        T u_cross_v_squared=u_cross_v.Magnitude_Squared();
        if(u_cross_v_squared>maximum_cross_squared){maximum_cross_squared_index=h;maximum_cross_squared=u_cross_v_squared;maximum_cross=u_cross_v;}}
    TV u,v;VECTOR<int,4> spring_nodes=Spring_Nodes(maximum_cross_squared_index,element_nodes);
    if(maximum_cross_squared_index<4){u=X(spring_nodes[2])-X(spring_nodes[1]);v=X(spring_nodes[3])-X(spring_nodes[1]);}
    else{u=X(spring_nodes[1])-X(spring_nodes[0]);v=X(spring_nodes[3])-X(spring_nodes[2]);}
    T u_length_squared=u.Magnitude_Squared(),v_length_squared=v.Magnitude_Squared();
    if(abs(maximum_cross_squared)<sqr(sin((T)pi/(T)180))*u_length_squared*v_length_squared){
        VECTOR<int,2> edge_index;T max_distance_squared=0;
        for(int i=0;i<3;i++) for(int j=i+1;j<4;j++){
            T distance_squared=(X(element_nodes[i])-X(element_nodes[j])).Magnitude_Squared();
            if(distance_squared>max_distance_squared){edge_index=VECTOR<int,2>(i,j);max_distance_squared=distance_squared;}}
            VECTOR<int,2> other_nodes=element_nodes.Remove_Index(edge_index[1]).Remove_Index(edge_index[0]);
            if((X(element_nodes[edge_index[0]])-X(other_nodes[1])).Magnitude_Squared()<(X(element_nodes[edge_index[0]])-X(other_nodes[0])).Magnitude_Squared())
                other_nodes=other_nodes.Reversed();
            VECTOR<int,2> edge1(element_nodes[edge_index[0]],other_nodes[1]),edge2(other_nodes[0],element_nodes[edge_index[1]]);
            bool found=false;VECTOR<int,4> spring_nodes;
            for(maximum_cross_squared_index=4;maximum_cross_squared_index<7;maximum_cross_squared_index++){
                spring_nodes=Spring_Nodes(maximum_cross_squared_index,element_nodes);
                if(spring_nodes.Slice<0,1>()==edge1 || spring_nodes.Slice<2,3>()==edge1 || spring_nodes.Slice<0,1>().Reversed()==edge1 || spring_nodes.Slice<2,3>().Reversed()==edge1){
                    found=true;break;}}
            PHYSBAM_ASSERT(found);
            minimum_normal=(X(edge1[0])-X(edge2[1])).Orthogonal_Vector().Normalized();
            minimum_signed_distance=0;
            TV midpoint=(T).5*(X(edge1[1])+X(edge2[0]));
            weights=VECTOR<T,3>(SEGMENT_3D<T>(X(spring_nodes[0]),X(spring_nodes[1])).Interpolation_Fraction(midpoint),
                SEGMENT_3D<T>(X(spring_nodes[2]),X(spring_nodes[3])).Interpolation_Fraction(midpoint),0);}
    else{
        minimum_normal=maximum_cross.Normalized();
        minimum_signed_distance=TV::Dot_Product(minimum_normal,X(spring_nodes[0])-X(spring_nodes[2]));
        if(maximum_cross_squared_index<4){
            weights=TRIANGLE_3D<T>::Clamped_Barycentric_Coordinates(X(spring_nodes[0]),X(spring_nodes[1]),X(spring_nodes[2]),X(spring_nodes[3]));}
        else if(maximum_cross_squared_index<7){
            VECTOR<T,2> dummy_weights;
            SEGMENT_3D<T>(X(spring_nodes[0]),X(spring_nodes[1])).Shortest_Vector_Between_Segments(SEGMENT_3D<T>(X(spring_nodes[2]),X(spring_nodes[3])),dummy_weights);
            weights=dummy_weights.Append(0);}}
    return maximum_cross_squared_index;

}
//#####################################################################
// Function Highlight_Nodes_Of_Minimum_Valence
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Highlight_Nodes_Of_Minimum_Valence() const
{
    bool neighbor_nodes_defined=mesh->neighbor_nodes!=0;if(!neighbor_nodes_defined) mesh->Initialize_Neighbor_Nodes();
    for(int p=0;p<particles->Size();p++) if((*mesh->neighbor_nodes)(p).m==minimum_valence)
        OPENGL_SHAPES::Draw_Dot(particles->X(p),OPENGL_COLOR(1,1,0),4);
    if(!neighbor_nodes_defined){delete mesh->neighbor_nodes;mesh->neighbor_nodes=0;}
}
//#####################################################################
// Function Highlight_Current_Node
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Highlight_Current_Node() const
{
    OPENGL_SHAPES::Draw_Dot(particles->X(current_node),OPENGL_COLOR(1,0,1),7);
}
//#####################################################################
// Function Highlight_Boundary_Normal_Vectors of Current_Tetrahedron
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Highlight_Boundary_Normal_Vectors_Of_Current_Tetrahedron() const
{
    int t=current_tetrahedron,i,j,k,l;mesh->elements(t).Get(i,j,k,l);
    VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k),xl=particles->X(l);
    OPENGL_SHAPES::Draw_Dot(((T)1/3)*(xj+xl+xk),OPENGL_COLOR((float).7,0,0),5);
    OPENGL_SHAPES::Draw_Vector(((T)1/3)*(xj+xl+xk),PLANE<T>::Normal(xj,xl,xk),OPENGL_COLOR(1,0,1),2);
    OPENGL_SHAPES::Draw_Vector(((T)1/3)*(xi+xk+xl),PLANE<T>::Normal(xi,xk,xl),OPENGL_COLOR(1,0,1),1);
    OPENGL_SHAPES::Draw_Vector(((T)1/3)*(xi+xl+xj),PLANE<T>::Normal(xi,xl,xj),OPENGL_COLOR(1,0,1),1);
    OPENGL_SHAPES::Draw_Vector(((T)1/3)*(xi+xj+xk),PLANE<T>::Normal(xi,xj,xk),OPENGL_COLOR(1,0,1),1);
}
//#####################################################################
// Function Highlight_Boundary_Nodes_Of_Current_Tetrahedron
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Highlight_Boundary_Nodes_Of_Current_Tetrahedron() const
{
    int t=current_tetrahedron,i,j,k,l;mesh->elements(t).Get(i,j,k,l);
    VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k),xl=particles->X(l);
    if((*mesh->node_on_boundary)(i)) OPENGL_SHAPES::Draw_Dot(xi,OPENGL_COLOR(1,0,0));
    if((*mesh->node_on_boundary)(j)) OPENGL_SHAPES::Draw_Dot(xj,OPENGL_COLOR(1,0,0));
    if((*mesh->node_on_boundary)(k)) OPENGL_SHAPES::Draw_Dot(xk,OPENGL_COLOR(1,0,0));
    if((*mesh->node_on_boundary)(l)) OPENGL_SHAPES::Draw_Dot(xl,OPENGL_COLOR(1,0,0));
}
//#####################################################################
// Function Draw_Wireframe_Mesh
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_Wireframe_Mesh(const TETRAHEDRON_MESH& tetrahedron_mesh) const
{
    glDisable(GL_LIGHTING);
    material.diffuse.Send_To_GL_Pipeline();
    ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
    for(int t=0;t<tetrahedron_mesh.elements.m;t++){
        int i,j,k,l;tetrahedron_mesh.elements(t).Get(i,j,k,l);
        OpenGL_Line(particles->X(i),particles->X(j),vertices);
        OpenGL_Line(particles->X(j),particles->X(k),vertices);
        OpenGL_Line(particles->X(k),particles->X(l),vertices);
        OpenGL_Line(particles->X(l),particles->X(i),vertices);
        OpenGL_Line(particles->X(i),particles->X(k),vertices);
        OpenGL_Line(particles->X(j),particles->X(l),vertices);}
    OpenGL_Draw_Arrays(GL_LINES,3,vertices);
    glEnable(GL_LIGHTING);
}
//#####################################################################
// Function Draw_Current_Tetrahedron
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_Current_Tetrahedron() const
{
    int t=current_tetrahedron,i,j,k,l;mesh->elements(t).Get(i,j,k,l);
    VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k),xl=particles->X(l);
    TETRAHEDRON<T> tet(xi,xj,xk,xl);
    ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
    ARRAY<GLfloat> normals;
    for(int i=0;i<4;i++){
        for(int plane_vertices=0;plane_vertices<3;plane_vertices++)
            OpenGL_Normal(tet.triangle(i).Normal(),normals);
        OpenGL_Triangle(tet.triangle(i).X,vertices);}
    OpenGL_Draw_Arrays_With_Normals(GL_TRIANGLES,3,vertices,normals);
}
//#####################################################################
// Function Set_Color_From_Aspect_Ratio
//#####################################################################
template<class T>
void Set_Color_From_Aspect_Ratio(const TRIANGLE_3D<T>& triangle)
{
    if(abs(triangle.Aspect_Ratio()-sqrt((T)2))<1e-3) OPENGL_MATERIAL::Plastic(OPENGL_COLOR(1.f,.1f,.3f)).Send_To_GL_Pipeline();
    else OPENGL_MATERIAL::Plastic(OPENGL_COLOR(.5f,1.f,.3f)).Send_To_GL_Pipeline();
}
//#####################################################################
// Function Draw_Subset
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_Subset() const
{
    OPENGL_MATERIAL::Plastic(OPENGL_COLOR(VECTOR<double,3>((T)1,(T).7,(T).7))).Send_To_GL_Pipeline();
    for(int tet_index=0;tet_index<subset.m;tet_index++){
        int t=subset(tet_index);int i,j,k,l;mesh->elements(t).Get(i,j,k,l);
        VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k),xl=particles->X(l);
        TETRAHEDRON<T> tet(xi,xj,xk,xl);
        ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;ARRAY<GLfloat> normals;
        //Set_Color_From_Aspect_Ratio(tet.triangle1);
        for(int i=0;i<4;i++){
            for(int plane_vertices=0;plane_vertices<3;plane_vertices++)
                OpenGL_Normal(tet.triangle(i).Normal(),normals);
            OpenGL_Triangle(tet.triangle(i).X,vertices);}
        OpenGL_Draw_Arrays_With_Normals(GL_TRIANGLES,3,vertices,normals);}
    material.Send_To_GL_Pipeline();
}
//#####################################################################
// Function Draw_Subset_Triangles
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_Subset_Triangles() const
{
    OPENGL_MATERIAL::Plastic(OPENGL_COLOR(VECTOR<double,3>((T).7,(T).7,(T)1))).Send_To_GL_Pipeline();
    ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;ARRAY<GLfloat> normals;
    for(int t=0;t<subset_triangles.m;t++){
        int i,j,k,tri;tri=subset_triangles(t);
        mesh->boundary_mesh->elements(tri).Get(i,j,k);
        VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k);
        for(int plane_vertices=0;plane_vertices<3;plane_vertices++) OpenGL_Normal(PLANE<T>::Normal(xi,xj,xk),normals);
        OpenGL_Triangle(xi,xj,xk,vertices);}
    OpenGL_Draw_Arrays_With_Normals(GL_TRIANGLES,3,vertices,normals);
    material.Send_To_GL_Pipeline();
}
//#####################################################################
// Function Draw_Subset_Particles
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_Subset_Particles() const
{
    for(int p=0;p<subset_particles.m;p++) OPENGL_SHAPES::Draw_Dot(particles->X(subset_particles(p)),OPENGL_COLOR(1,1,0));
}
//#####################################################################
// Function Set_Color_From_Spectrum
//#####################################################################
template<class T>
void Set_Color_From_Spectrum(int tetrahedron,const ARRAY<T>& spectrum,T spectrum_max,T spectrum_min)
{
    T red=spectrum(tetrahedron)-spectrum_min/(spectrum_max-spectrum_min);
    PHYSBAM_ASSERT(0<=red && red<=1);
    OPENGL_MATERIAL::Plastic(OPENGL_COLOR((float)red,.1f,.1f)).Send_To_GL_Pipeline();
}
//#####################################################################
// Function Draw_In_Color_From_Spectrum
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_In_Color_From_Spectrum() const
{
    if(spectrum.m!=mesh->elements.m) PHYSBAM_FATAL_ERROR("Spectrum has incorrect size");
    T spectrum_max=spectrum.Max(),spectrum_min=spectrum.Min();
    for(int t=0;t<mesh->elements.m;t++){
        int i,j,k,l;mesh->elements(t).Get(i,j,k,l);
        VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k),xl=particles->X(l);
        TETRAHEDRON<T> tet(xi,xj,xk,xl);
        ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;ARRAY<GLfloat> normals;
        for(int i=0;i<4;i++){
            Set_Color_From_Spectrum(t,spectrum,spectrum_max,spectrum_min);
            for(int plane_vertices=0;plane_vertices<3;plane_vertices++)
                OpenGL_Normal(tet.triangle(i).Normal(),normals);
            OpenGL_Triangle(tet.triangle(i).X,vertices);}
        OpenGL_Draw_Arrays_With_Normals(GL_TRIANGLES,3,vertices,normals);}
    material.Send_To_GL_Pipeline();
}
//#####################################################################
// Function Draw_In_Color_From_Spectrum
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_In_Color_From_Color_Map() const
{
    glPushAttrib(GL_ENABLE_BIT|GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    for(int t=0;t<mesh->elements.m;t++){
        (*color_map)(t).Send_To_GL_Pipeline();
        int i,j,k,l;mesh->elements(t).Get(i,j,k,l);
        VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k),xl=particles->X(l);
        TETRAHEDRON<T> tet(xi,xj,xk,xl);
        ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;ARRAY<GLfloat> normals;
        for(int i=0;i<4;i++){
            for(int plane_vertices=0;plane_vertices<3;plane_vertices++)
                OpenGL_Normal(tet.triangle(i).Normal(),normals);
            OpenGL_Triangle(tet.triangle(i).X,vertices);}
        OpenGL_Draw_Arrays_With_Normals(GL_TRIANGLES,3,vertices,normals);}
    glPopAttrib();
}
//#####################################################################
// Function Draw_Boundary_Triangles
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_Boundary_Triangles(const TETRAHEDRON_MESH& tetrahedron_mesh) const
{
    if(use_inverted_material){
        glDisable(GL_CULL_FACE);
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
        material.Send_To_GL_Pipeline(GL_FRONT);
        inverted_material.Send_To_GL_Pipeline(GL_BACK);}
    else material.Send_To_GL_Pipeline();

    ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;ARRAY<GLfloat> normals;
    if(smooth_normals)
        for(int t=0;t<tetrahedron_mesh.boundary_mesh->elements.m;t++){
            int i,j,k;tetrahedron_mesh.boundary_mesh->elements(t).Get(i,j,k);
            OpenGL_Normal((*vertex_normals)(i),normals);OpenGL_Vertex(particles->X(i),vertices);
            OpenGL_Normal((*vertex_normals)(j),normals);OpenGL_Vertex(particles->X(j),vertices);
            OpenGL_Normal((*vertex_normals)(k),normals);OpenGL_Vertex(particles->X(k),vertices);}
    else
        for(int t=0;t<tetrahedron_mesh.boundary_mesh->elements.m;t++){
            int i,j,k;tetrahedron_mesh.boundary_mesh->elements(t).Get(i,j,k);
            VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k);
            for(int plane_vertices=0;plane_vertices<3;plane_vertices++) OpenGL_Normal(PLANE<T>::Normal(xi,xj,xk),normals);
            OpenGL_Triangle(xi,xj,xk,vertices);}
    OpenGL_Draw_Arrays_With_Normals(GL_TRIANGLES,3,vertices,normals);

    if(use_inverted_material){
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,0);
        glEnable(GL_CULL_FACE);}
}
//#####################################################################
// Function Highlight_Current_Boundary_Triangle
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Highlight_Current_Boundary_Triangle() const
{
    glDisable(GL_LIGHTING);
    glColor3f(1,1,1);
    ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;ARRAY<GLfloat> normals;
    int i,j,k;mesh->boundary_mesh->elements(current_boundary_triangle).Get(i,j,k);
    VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k);
    for(int plane_vertices=0;plane_vertices<3;plane_vertices++) OpenGL_Normal(PLANE<T>::Normal(xi,xj,xk),normals);
    OpenGL_Triangle(xi,xj,xk,vertices);
    OpenGL_Draw_Arrays_With_Normals(GL_TRIANGLES,3,vertices,normals);
    OPENGL_SHAPES::Draw_Dot((T)1./3*(xi+xj+xk),OPENGL_COLOR(1,1,1),8);
    glEnable(GL_LIGHTING);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Bounding_Box() const
{
    return World_Space_Box(RANGE<VECTOR<T,3> >::Bounding_Box(particles->X));
}
//#####################################################################
// Function Turn_Smooth_Shading_On
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Turn_Smooth_Shading_On()
{
    if(!vertex_normals)Initialize_Vertex_Normals();smooth_normals=true;
}
//#####################################################################
// Function Turn_Smooth_Shading_Off
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Turn_Smooth_Shading_Off()
{
    smooth_normals=false;
}
//#####################################################################
// Function Initialize_Vertex_Normals
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Initialize_Vertex_Normals()
{
    delete vertex_normals;vertex_normals=new ARRAY<VECTOR<T,3> >(particles->Size());
    for(int t=0;t<mesh->boundary_mesh->elements.m;t++){
        int i,j,k;mesh->boundary_mesh->elements(t).Get(i,j,k);
        VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k);
        VECTOR<T,3> normal=PLANE<T>::Normal(xi,xj,xk);
        (*vertex_normals)(i)+=normal;(*vertex_normals)(j)+=normal;(*vertex_normals)(k)+=normal;}
    for(int p=0;p<particles->Size();p++)(*vertex_normals)(p).Normalize();
}
//#####################################################################
// Function Update_Cutaway_Plane
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Update_Cutaway_Plane()
{
    TETRAHEDRALIZED_VOLUME<T> tetrahedralized_volume(*mesh,*particles);
    tetrahedralized_volume.Update_Bounding_Box();
    RANGE<VECTOR<T,3> > box=*tetrahedralized_volume.bounding_box;
    ARRAY<bool> inside(particles->Size());T threshold;
    switch(cutaway_mode){
        case 1:threshold=box.min_corner.x+cutaway_fraction*(box.max_corner.x-box.min_corner.x);for(int p=0;p<particles->Size();p++)inside(p)=particles->X(p).x<threshold;break;
        case 2:threshold=box.max_corner.x+cutaway_fraction*(box.min_corner.x-box.max_corner.x);for(int p=0;p<particles->Size();p++)inside(p)=particles->X(p).x>threshold;break;
        case 3:threshold=box.min_corner.y+cutaway_fraction*(box.max_corner.y-box.min_corner.y);for(int p=0;p<particles->Size();p++)inside(p)=particles->X(p).y<threshold;break;
        case 4:threshold=box.max_corner.y+cutaway_fraction*(box.min_corner.y-box.max_corner.y);for(int p=0;p<particles->Size();p++)inside(p)=particles->X(p).y>threshold;break;
        case 5:threshold=box.min_corner.z+cutaway_fraction*(box.max_corner.z-box.min_corner.z);for(int p=0;p<particles->Size();p++)inside(p)=particles->X(p).z<threshold;break;
        case 6:threshold=box.max_corner.z+cutaway_fraction*(box.min_corner.z-box.max_corner.z);for(int p=0;p<particles->Size();p++)inside(p)=particles->X(p).z>threshold;break;}
    cutaway_mesh.number_nodes=mesh->number_nodes;
    cutaway_mesh.elements.Remove_All();
    for(int t=0;t<mesh->elements.m;t++){
        int i,j,k,l;mesh->elements(t).Get(i,j,k,l);
        if(inside(i) || inside(j) || inside(k) || inside(l)) cutaway_mesh.elements.Append(VECTOR<int,4>(i,j,k,l));}
    cutaway_mesh.Initialize_Boundary_Mesh();
}
//#####################################################################
// Function Display_Subset
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Display_Subset()
{
    material.Send_To_GL_Pipeline();
    glColor3f((float).7,(float).7,(float).7);
    ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;ARRAY<GLfloat> normals;
    for(int tet_index=0;tet_index<subset.m;tet_index++){
        int t=subset(tet_index);
        int i,j,k,l;mesh->elements(t).Get(i,j,k,l);
        VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k),xl=particles->X(l);
        TETRAHEDRON<T> tet(xi,xj,xk,xl);
        for(int i=0;i<4;i++){
            for(int plane_vertices=0;plane_vertices<3;plane_vertices++)
                OpenGL_Normal(tet.triangle(i).Normal(),normals);
            OpenGL_Triangle(tet.triangle(i).X,vertices);}}
    OpenGL_Draw_Arrays_With_Normals(GL_TRIANGLES,3,vertices,normals);
    vertices.Resize(0);
    for(int t=0;t<mesh->boundary_mesh->elements.m;t++){
        int i,j,k;mesh->boundary_mesh->elements(t).Get(i,j,k);
        VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k);
        OpenGL_Line(xi,xj,vertices);
        OpenGL_Line(xj,xk,vertices);
        OpenGL_Line(xi,xk,vertices);}
    OpenGL_Draw_Arrays(GL_LINES,3,vertices);
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Get_Selection(GLuint* buffer,int buffer_size)
{
    OPENGL_SELECTION<T>* selection=0;
    if(buffer_size==2){
        if(buffer[0]==1)selection=Get_Vertex_Selection(buffer[1]);
        else if(buffer[0]==2)selection=Get_Tetrahedron_Selection(buffer[1]);}
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Highlight_Selection(OPENGL_SELECTION<T>* selection)
{
    delete current_selection;current_selection=0;
    // Make a copy of selection
    if(selection->type==OPENGL_SELECTION<T>::TETRAHEDRALIZED_VOLUME_VERTEX)
        current_selection=new OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_VERTEX<T>(this,((OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_VERTEX<T>*)selection)->index);
    else if(selection->type==OPENGL_SELECTION<T>::TETRAHEDRALIZED_VOLUME_TETRAHEDRON)
        current_selection=new OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_TETRAHEDRON<T>(this,((OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_TETRAHEDRON<T>*)selection)->index);
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Clear_Highlight()
{
    delete current_selection;current_selection=0;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION<T>* selection) const
{
    Print_Selection_Info(output_stream,selection,0);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION<T>* selection,MATRIX<T,4>* transform) const
{
    if(selection->type==OPENGL_SELECTION<T>::TETRAHEDRALIZED_VOLUME_VERTEX){
        int index=((OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_VERTEX<T>*)selection)->index;
        output_stream<<"Vertex "<<index<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(particles->X(index))<<std::endl;}
        particles->Print(output_stream,index);}
    else if(selection->type==OPENGL_SELECTION<T>::TETRAHEDRALIZED_VOLUME_TETRAHEDRON){
        int index=((OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_TETRAHEDRON<T>*)selection)->index;
        const VECTOR<int,4>& nodes=mesh->elements(index);
        output_stream<<"Tetrahedron "<<index<<" ("<<nodes[0]<<", "<<nodes[1]<<", "<<nodes[2]<<", "<<nodes[3]<<")"<<std::endl;
        output_stream<<"Signed Volume = "<<TETRAHEDRON<T>(particles->X.Subset(nodes)).Signed_Volume()<<std::endl;
        output_stream<<std::endl;
        output_stream<<"Vertex "<<nodes[0]<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(particles->X(nodes[0]))<<std::endl;}
        particles->Print(output_stream,nodes[0]);
        output_stream<<std::endl;
        output_stream<<"Vertex "<<nodes[1]<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(particles->X(nodes[1]))<<std::endl;}
        particles->Print(output_stream,nodes[1]);
        output_stream<<std::endl;
        output_stream<<"Vertex "<<nodes[2]<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(particles->X(nodes[2]))<<std::endl;}
        particles->Print(output_stream,nodes[2]);
        output_stream<<std::endl;
        output_stream<<"Vertex "<<nodes[3]<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(particles->X(nodes[3]))<<std::endl;}
        particles->Print(output_stream,nodes[3]);

        T distance;TV min_normal,weights;
        int spring=Find_Shortest_Spring(nodes,distance,min_normal,weights);
        output_stream<<"Shortest Altitude is ";
        if(spring<4) output_stream<<"point-face ";
        else if(spring<7) output_stream<<"edge-edge ";
        output_stream<<spring<<" nodes="<<Spring_Nodes(spring,nodes)<<std::endl;
        output_stream<<"   distance="<<distance<<" normal="<<min_normal<<" weights="<<weights<<std::endl;
    }

}
//#####################################################################
// Function Get_Vertex_Selection
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Get_Vertex_Selection(int index)
{
    return new OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_VERTEX<T>(this,index);
}
//#####################################################################
// Function Get_Tetrahedron_Selection
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Get_Tetrahedron_Selection(int index)
{
    return new OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_TETRAHEDRON<T>(this,index);
}
//#####################################################################
// Function Draw_Vertices_For_Selection
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_Vertices_For_Selection() const
{
    OPENGL_SELECTION<T>::Draw_Vertices_For_Selection(*mesh,*particles);
}
//#####################################################################
// Function Draw_Tetrahedra_For_Selection
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_Tetrahedra_For_Selection() const
{
#ifndef USE_OPENGLES
    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_CULL_FACE);
    glPushName(0);
    for(int t=0;t<mesh->elements.m;t++){
        glLoadName(t);
        ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
        int i,j,k,l;mesh->elements(t).Get(i,j,k,l);
        OpenGL_Triangle(particles->X(i),particles->X(j),particles->X(k),vertices);
        OpenGL_Triangle(particles->X(i),particles->X(k),particles->X(l),vertices);
        OpenGL_Triangle(particles->X(i),particles->X(l),particles->X(j),vertices);
        OpenGL_Triangle(particles->X(l),particles->X(k),particles->X(j),vertices);
        OpenGL_Draw_Arrays(GL_TRIANGLES,3,vertices);}
    glPopName();
    glPopAttrib();
#endif
}
//#####################################################################
// Selection object functions
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_VERTEX<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const OPENGL_TETRAHEDRALIZED_VOLUME<T>& volume=dynamic_cast<OPENGL_TETRAHEDRALIZED_VOLUME<T>&>(*object);
    return object->World_Space_Box(RANGE<VECTOR<T,3> >(volume.particles->X(index)));
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_TETRAHEDRON<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const OPENGL_TETRAHEDRALIZED_VOLUME<T>& volume=dynamic_cast<OPENGL_TETRAHEDRALIZED_VOLUME<T>&>(*object);
    return object->World_Space_Box(RANGE<VECTOR<T,3> >::Bounding_Box(volume.particles->X.Subset(volume.mesh->elements(index))));
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_TETRAHEDRALIZED_VOLUME<float>;
template class OPENGL_TETRAHEDRALIZED_VOLUME<double>;
}

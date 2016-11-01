//#####################################################################
// Copyright 2005-2009, Zhaosheng Bao, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Sergey Koltakov, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Matrices/MATRIX_4X4.h>
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
OPENGL_TETRAHEDRALIZED_VOLUME(STREAM_TYPE stream_type,TETRAHEDRON_MESH* mesh_input,GEOMETRY_PARTICLES<VECTOR<T,3> >* particles_input,const OPENGL_MATERIAL& material_input,
    const OPENGL_MATERIAL& inverted_material_input,bool initialize,ARRAY<OPENGL_COLOR>* color_map_input)
    :OPENGL_OBJECT<T>(stream_type),material(material_input),inverted_material(inverted_material_input),
    use_inverted_material(true),mesh(mesh_input),particles(particles_input),current_tetrahedron(0),current_node(0),
    current_boundary_triangle(0),boundary_only(true),draw_subsets(false),cutaway_mode(0),cutaway_fraction((T).5),
    color_map(color_map_input),smooth_normals(false),vertex_normals(0),selected_vertex(-1),selected_tet(-1)
{
    if(initialize){
        if(!mesh->boundary_mesh)
            mesh->Initialize_Boundary_Mesh(); // Neighboring nodes is no longer initialized here to conserve memory.
        if(!mesh->node_on_boundary)
            mesh->Initialize_Node_On_Boundary();
        minimum_valence=mesh->Minimum_Valence();
        mesh->Initialize_Boundary_Nodes();}
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
    glGetIntegerv(GL_RENDER_MODE,&mode);
    glDisable(GL_CULL_FACE);
    if(mode==GL_SELECT){
        glPushName(0);
        Draw_Vertices_For_Selection();
        glLoadName(1);
        Draw_Tetrahedra_For_Selection();
        glPopName();}
    else if(cutaway_mode){
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

        if(selected_vertex>=0)
            OPENGL_SELECTION::Draw_Highlighted_Vertex(particles->X(selected_vertex),selected_vertex);
        else if(selected_tet>=0){
            const VECTOR<int,4>& element_nodes=mesh->elements(selected_tet);
            ARRAY_VIEW<const TV> X(particles->X);
            VECTOR<TV,4> Y(X.Subset(element_nodes));
            OPENGL_SELECTION::Draw_Highlighted_Tetrahedron_Boundary(Y(0),Y(1),Y(2),Y(3),selected_tet);
            T distance;TV min_normal,weights;
            int spring=Find_Shortest_Spring(element_nodes,distance,min_normal,weights);
            VECTOR<int,4> spring_nodes=Spring_Nodes(spring,element_nodes);

            glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
            glDisable(GL_LIGHTING);
            glLineWidth(OPENGL_PREFERENCES::highlighted_line_width*2);
            OPENGL_COLOR colors[]={OPENGL_COLOR::Red(),OPENGL_COLOR::Blue(),OPENGL_COLOR::Green(),OPENGL_COLOR::Yellow(),OPENGL_COLOR::Cyan(),OPENGL_COLOR::Magenta(),
                                   OPENGL_COLOR::White()};
            if(spring>=0){
                colors[spring].Send_To_GL_Pipeline();
                OpenGL_Begin(GL_LINES);
                if(spring<4) OpenGL_Line(X(spring_nodes[0]),TRIANGLE_3D<T>(X.Subset(spring_nodes.Remove_Index(0))).Point_From_Barycentric_Coordinates(weights));
                else if(spring<7) OpenGL_Line((1-weights.x)*X(spring_nodes[0])+weights.x*X(spring_nodes[1]),(1-weights.y)*X(spring_nodes[2])+weights.y*X(spring_nodes[3]));
                OpenGL_End();
                glPopAttrib();}}}
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
    OpenGL_Begin(GL_LINES);
    for(int t=0;t<tetrahedron_mesh.elements.m;t++){
        int i,j,k,l;tetrahedron_mesh.elements(t).Get(i,j,k,l);
        OpenGL_Line(particles->X(i),particles->X(j));
        OpenGL_Line(particles->X(j),particles->X(k));
        OpenGL_Line(particles->X(k),particles->X(l));
        OpenGL_Line(particles->X(l),particles->X(i));
        OpenGL_Line(particles->X(i),particles->X(k));
        OpenGL_Line(particles->X(j),particles->X(l));}
    OpenGL_End();
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
    OpenGL_Begin(GL_TRIANGLES);
    for(int i=0;i<4;i++){
        for(int plane_vertices=0;plane_vertices<3;plane_vertices++)
            OpenGL_Normal(tet.triangle(i).Normal());
        OpenGL_Triangle(tet.triangle(i).X);}
    OpenGL_End();
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
        OpenGL_Begin(GL_TRIANGLES);
        //Set_Color_From_Aspect_Ratio(tet.triangle1);
        for(int i=0;i<4;i++){
            for(int plane_vertices=0;plane_vertices<3;plane_vertices++)
                OpenGL_Normal(tet.triangle(i).Normal());
            OpenGL_Triangle(tet.triangle(i).X);}
        OpenGL_End();}
    material.Send_To_GL_Pipeline();
}
//#####################################################################
// Function Draw_Subset_Triangles
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_Subset_Triangles() const
{
    OPENGL_MATERIAL::Plastic(OPENGL_COLOR(VECTOR<double,3>((T).7,(T).7,(T)1))).Send_To_GL_Pipeline();
    OpenGL_Begin(GL_TRIANGLES);
    for(int t=0;t<subset_triangles.m;t++){
        int i,j,k,tri;tri=subset_triangles(t);
        mesh->boundary_mesh->elements(tri).Get(i,j,k);
        VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k);
        for(int plane_vertices=0;plane_vertices<3;plane_vertices++) OpenGL_Normal(PLANE<T>::Normal(xi,xj,xk));
        OpenGL_Triangle(xi,xj,xk);}
    OpenGL_End();
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
        OpenGL_Begin(GL_TRIANGLES);
        for(int i=0;i<4;i++){
            Set_Color_From_Spectrum(t,spectrum,spectrum_max,spectrum_min);
            for(int plane_vertices=0;plane_vertices<3;plane_vertices++)
                OpenGL_Normal(tet.triangle(i).Normal());
            OpenGL_Triangle(tet.triangle(i).X);}
        OpenGL_End();}
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
        OpenGL_Begin(GL_TRIANGLES);
        for(int i=0;i<4;i++){
            for(int plane_vertices=0;plane_vertices<3;plane_vertices++)
                OpenGL_Normal(tet.triangle(i).Normal());
            OpenGL_Triangle(tet.triangle(i).X);}
        OpenGL_End();}
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

    OpenGL_Begin(GL_TRIANGLES);
    if(smooth_normals)
        for(int t=0;t<tetrahedron_mesh.boundary_mesh->elements.m;t++){
            int i,j,k;tetrahedron_mesh.boundary_mesh->elements(t).Get(i,j,k);
            OpenGL_Normal((*vertex_normals)(i));OpenGL_Vertex(particles->X(i));
            OpenGL_Normal((*vertex_normals)(j));OpenGL_Vertex(particles->X(j));
            OpenGL_Normal((*vertex_normals)(k));OpenGL_Vertex(particles->X(k));}
    else
        for(int t=0;t<tetrahedron_mesh.boundary_mesh->elements.m;t++){
            int i,j,k;tetrahedron_mesh.boundary_mesh->elements(t).Get(i,j,k);
            VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k);
            for(int plane_vertices=0;plane_vertices<3;plane_vertices++) OpenGL_Normal(PLANE<T>::Normal(xi,xj,xk));
            OpenGL_Triangle(xi,xj,xk);}
    OpenGL_End();

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
    OpenGL_Begin(GL_TRIANGLES);
    int i,j,k;mesh->boundary_mesh->elements(current_boundary_triangle).Get(i,j,k);
    VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k);
    for(int plane_vertices=0;plane_vertices<3;plane_vertices++) OpenGL_Normal(PLANE<T>::Normal(xi,xj,xk));
    OpenGL_Triangle(xi,xj,xk);
    OpenGL_End();
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
    OpenGL_Begin(GL_TRIANGLES);
    for(int tet_index=0;tet_index<subset.m;tet_index++){
        int t=subset(tet_index);
        int i,j,k,l;mesh->elements(t).Get(i,j,k,l);
        VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k),xl=particles->X(l);
        TETRAHEDRON<T> tet(xi,xj,xk,xl);
        for(int i=0;i<4;i++){
            for(int plane_vertices=0;plane_vertices<3;plane_vertices++)
                OpenGL_Normal(tet.triangle(i).Normal());
            OpenGL_Triangle(tet.triangle(i).X);}}
    OpenGL_End();
    OpenGL_Begin(GL_LINES);
    for(int t=0;t<mesh->boundary_mesh->elements.m;t++){
        int i,j,k;mesh->boundary_mesh->elements(t).Get(i,j,k);
        VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k);
        OpenGL_Line(xi,xj);
        OpenGL_Line(xj,xk);
        OpenGL_Line(xi,xk);}
    OpenGL_End();
}
//#####################################################################
// Function Get_Selection_Priority
//#####################################################################
template<class T> int OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Get_Selection_Priority(ARRAY_VIEW<GLuint> indices)
{
    if(!indices.m) return -1;
    PHYSBAM_ASSERT(indices.m==2);
    const static int priority[]={85,84};
    return priority[indices(0)];
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> bool OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers)
{
    if(indices(0)==0) selected_vertex=indices(1);
    else if(indices(0)==1) selected_tet=indices(1);
    else return false;
    return true;
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Clear_Selection()
{
    selected_vertex=-1;
    selected_tet=-1;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Print_Selection_Info(std::ostream &output_stream) const
{
    Print_Selection_Info(output_stream,0);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Print_Selection_Info(std::ostream &output_stream,MATRIX<T,4>* transform) const
{
    if(selected_vertex>=0){
        output_stream<<"Vertex "<<selected_vertex<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(particles->X(selected_vertex))<<std::endl;}
        particles->Print(output_stream,selected_vertex);}
    else if(selected_tet>=0){
        const VECTOR<int,4>& nodes=mesh->elements(selected_tet);
        output_stream<<"Tetrahedron "<<selected_tet<<" ("<<nodes[0]<<", "<<nodes[1]<<", "<<nodes[2]<<", "<<nodes[3]<<")"<<std::endl;
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
        output_stream<<"   distance="<<distance<<" normal="<<min_normal<<" weights="<<weights<<std::endl;}
}
//#####################################################################
// Function Draw_Vertices_For_Selection
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_Vertices_For_Selection() const
{
    OPENGL_SELECTION::Draw_Vertices_For_Selection(*mesh,*particles);
}
//#####################################################################
// Function Draw_Tetrahedra_For_Selection
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_Tetrahedra_For_Selection() const
{
    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_CULL_FACE);
    glPushName(0);
    for(int t=0;t<mesh->elements.m;t++){
        glLoadName(t);
        OpenGL_Begin(GL_TRIANGLES);
        int i,j,k,l;mesh->elements(t).Get(i,j,k,l);
        OpenGL_Triangle(particles->X(i),particles->X(j),particles->X(k));
        OpenGL_Triangle(particles->X(i),particles->X(k),particles->X(l));
        OpenGL_Triangle(particles->X(i),particles->X(l),particles->X(j));
        OpenGL_Triangle(particles->X(l),particles->X(k),particles->X(j));
        OpenGL_End();}
    glPopName();
    glPopAttrib();
}
//#####################################################################
// Function Selection_Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Selection_Bounding_Box() const
{
    if(selected_vertex>=0) return World_Space_Box(RANGE<TV>(particles->X(selected_vertex)));
    if(selected_tet>=0) return World_Space_Box(RANGE<TV>::Bounding_Box(particles->X.Subset(mesh->elements(selected_tet))));
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_TETRAHEDRALIZED_VOLUME<float>;
template class OPENGL_TETRAHEDRALIZED_VOLUME<double>;
}

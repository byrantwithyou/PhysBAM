//#####################################################################
// Copyright 2003-2008, Eran Guendelman, Sergey Koltakov, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
using namespace PhysBAM;
//#####################################################################
// Function Draw_Highlighted_Vertex
//#####################################################################
template<class TV> void PhysBAM::OPENGL_SELECTION::
Draw_Highlighted_Vertex(const TV& position,int id,const OPENGL_COLOR& color)
{
    glPushAttrib(GL_POINT_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    glPointSize(OPENGL_PREFERENCES::highlighted_point_size);
    color.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_POINTS);
    OpenGL_Vertex(position);
    OpenGL_End();
    if(id>=0) OpenGL_String(position,LOG::sprintf("   %d",id));
    glPopAttrib();
}
//#####################################################################
// Function Draw_Highlighted_Segment
//#####################################################################
template<class TV> void PhysBAM::OPENGL_SELECTION::
Draw_Highlighted_Segment(const TV& x0,const TV& x1,int id,const OPENGL_COLOR& color)
{
    glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    glLineWidth(OPENGL_PREFERENCES::highlighted_line_width);
    color.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_LINES);
    OpenGL_Line(x0,x1);
    OpenGL_End();
    if(id>=0) OpenGL_String((typename TV::SCALAR).5*(x0+x1),LOG::sprintf("   %d",id));
    glPopAttrib();
}
//#####################################################################
// Function Draw_Highlighted_Curve
//#####################################################################
template<class TV> void PhysBAM::OPENGL_SELECTION::
Draw_Highlighted_Curve(const ARRAY<VECTOR<TV,2> >& X,int id,const OPENGL_COLOR& color)
{
    TV total=TV();
    glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    glLineWidth(OPENGL_PREFERENCES::highlighted_line_width);
    color.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_LINES);
    for(int i=0;i<X.m;i++) {OpenGL_Line(X(i).x,X(i).y);total+=X(i).x;}
    OpenGL_End();
    if(id>=0) OpenGL_String((typename TV::SCALAR)1./X.m*(total),LOG::sprintf("   %d",id));
    glPopAttrib();
}
//#####################################################################
// Function Draw_Highlighted_Triangle_Boundary
//#####################################################################
template<class TV> void PhysBAM::OPENGL_SELECTION::
Draw_Highlighted_Triangle_Boundary(const TV& x0,const TV& x1,const TV& x2,int id,const OPENGL_COLOR& color)
{
    glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    glLineWidth(OPENGL_PREFERENCES::highlighted_line_width);
    color.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_LINE_LOOP);
    OpenGL_Vertex(x0);
    OpenGL_Vertex(x1);
    OpenGL_Vertex(x2);
    OpenGL_End();
    if(id>=0) OpenGL_String(((typename TV::SCALAR)1/3)*(x0+x1+x2),LOG::sprintf("%d",id));
    glPopAttrib();
}
//#####################################################################
// Function Draw_Highlighted_Tetrahedron_Boundary
//#####################################################################
template<class TV> void PhysBAM::OPENGL_SELECTION::
Draw_Highlighted_Tetrahedron_Boundary(const TV& x0,const TV& x1,const TV& x2,const TV& x3,int id,const OPENGL_COLOR& color)
{
    typedef typename TV::SCALAR T;
    glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    glLineWidth(OPENGL_PREFERENCES::highlighted_line_width);
    color.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_LINE_STRIP);
    OpenGL_Vertex(x0);
    OpenGL_Vertex(x1);
    OpenGL_Vertex(x2);
    OpenGL_Vertex(x0);
    OpenGL_Vertex(x3);
    OpenGL_Vertex(x1);
    OpenGL_End();
    OpenGL_Begin(GL_LINES);
    OpenGL_Vertex(x3);
    OpenGL_Vertex(x2);
    OpenGL_End();
    if(id>=0) OpenGL_String((T).5*(x0+x1+x2+x3),LOG::sprintf("%d",id));
    glPopAttrib();
}
//#####################################################################
// Function Draw_Highlighted_Quad
//#####################################################################
template<class TV> void PhysBAM::OPENGL_SELECTION::
Draw_Highlighted_Quad(const TV& x00,const TV& x11,const OPENGL_COLOR& color)
{
    glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    glLineWidth(OPENGL_PREFERENCES::highlighted_line_width);
    color.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_LINE_LOOP);
    OpenGL_Quad_2D(x00,x11);
    OpenGL_End();
    glPopAttrib();
}
//#####################################################################
// Function Draw_Highlighted_Quad
//#####################################################################
template<class TV> void PhysBAM::OPENGL_SELECTION::
Draw_Highlighted_Quad(const TV& node1,const TV& node2,const TV& node3,const TV& node4,const OPENGL_COLOR& color)
{
    glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    glLineWidth(OPENGL_PREFERENCES::highlighted_line_width);
    color.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_LINE_LOOP);
    OpenGL_Vertex(node1);OpenGL_Vertex(node2);OpenGL_Vertex(node4);OpenGL_Vertex(node3);
    OpenGL_End();
    glPopAttrib();
}
//#####################################################################
// Function Draw_Highlighted_Box
//#####################################################################
template<class TV> void PhysBAM::OPENGL_SELECTION::
Draw_Highlighted_Box(const TV& x000,const TV& x111,const OPENGL_COLOR& color)
{
    glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    glLineWidth(OPENGL_PREFERENCES::highlighted_line_width);
    color.Send_To_GL_Pipeline();
    TV x001(x000.x,x000.y,x111.z),x010(x000.x,x111.y,x000.z),x011(x000.x,x111.y,x111.z),
        x100(x111.x,x000.y,x000.z),x101(x111.x,x000.y,x111.z),x110(x111.x,x111.y,x000.z);
    OpenGL_Begin(GL_LINE_LOOP);
    OpenGL_Vertex(x000);
    OpenGL_Vertex(x010);
    OpenGL_Vertex(x011);
    OpenGL_Vertex(x001);
    OpenGL_End();
    OpenGL_Begin(GL_LINE_LOOP);
    OpenGL_Vertex(x100);
    OpenGL_Vertex(x110);
    OpenGL_Vertex(x111);
    OpenGL_Vertex(x101);
    OpenGL_End();
    OpenGL_Begin(GL_LINES);
    OpenGL_Vertex(x000);
    OpenGL_Vertex(x100);
    OpenGL_Vertex(x010);
    OpenGL_Vertex(x110);
    OpenGL_Vertex(x011);
    OpenGL_Vertex(x111);
    OpenGL_Vertex(x001);
    OpenGL_Vertex(x101);
    OpenGL_End();
    glPopAttrib();
}
//#####################################################################
// Function Draw_Vertices_For_Selection
//#####################################################################
template<class TV,int d> void PhysBAM::OPENGL_SELECTION::
Draw_Vertices_For_Selection(const SIMPLEX_MESH<d>& mesh,const GEOMETRY_PARTICLES<TV>& particles)
{
    glPushAttrib(GL_POINT_BIT);
    glPointSize(OPENGL_PREFERENCES::selection_point_size);
    glPushName(0);
    ARRAY<int> particles_in_mesh;
    Get_Unique(particles_in_mesh,mesh.elements.Flattened());
    ARRAY<typename OPENGL_POLICY<typename TV::SCALAR>::T_GL >vertices;
    for(int i=0;i<particles_in_mesh.m;i++){const int p=particles_in_mesh(i);
        glLoadName(p);
        OpenGL_Begin(GL_POINTS);
        OpenGL_Vertex(particles.X(p));
        OpenGL_End();
    }
    glPopName();
    glPopAttrib();
}
namespace PhysBAM{
template void OPENGL_SELECTION::Draw_Highlighted_Box<VECTOR<double,3> >(VECTOR<double,3> const&,VECTOR<double,3> const&,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Highlighted_Box<VECTOR<float,3> >(VECTOR<float,3> const&,VECTOR<float,3> const&,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Highlighted_Curve<VECTOR<double,3> >(ARRAY<VECTOR<VECTOR<double,3>,2>,int> const&,int,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Highlighted_Curve<VECTOR<float,3> >(ARRAY<VECTOR<VECTOR<float,3>,2>,int> const&,int,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Highlighted_Quad<VECTOR<double,2> >(VECTOR<double,2> const&,VECTOR<double,2> const&,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Highlighted_Quad<VECTOR<float,2> >(VECTOR<float,2> const&,VECTOR<float,2> const&,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Highlighted_Segment<VECTOR<double,2> >(VECTOR<double,2> const&,VECTOR<double,2> const&,int,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Highlighted_Segment<VECTOR<double,3> >(VECTOR<double,3> const&,VECTOR<double,3> const&,int,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Highlighted_Segment<VECTOR<float,2> >(VECTOR<float,2> const&,VECTOR<float,2> const&,int,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Highlighted_Segment<VECTOR<float,3> >(VECTOR<float,3> const&,VECTOR<float,3> const&,int,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Highlighted_Tetrahedron_Boundary<VECTOR<double,3> >(VECTOR<double,3> const&,VECTOR<double,3> const&,VECTOR<double,3> const&,VECTOR<double,3> const&,int,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Highlighted_Tetrahedron_Boundary<VECTOR<float,3> >(VECTOR<float,3> const&,VECTOR<float,3> const&,VECTOR<float,3> const&,VECTOR<float,3> const&,int,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Highlighted_Triangle_Boundary<VECTOR<double,2> >(VECTOR<double,2> const&,VECTOR<double,2> const&,VECTOR<double,2> const&,int,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Highlighted_Triangle_Boundary<VECTOR<double,3> >(VECTOR<double,3> const&,VECTOR<double,3> const&,VECTOR<double,3> const&,int,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Highlighted_Triangle_Boundary<VECTOR<float,2> >(VECTOR<float,2> const&,VECTOR<float,2> const&,VECTOR<float,2> const&,int,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Highlighted_Triangle_Boundary<VECTOR<float,3> >(VECTOR<float,3> const&,VECTOR<float,3> const&,VECTOR<float,3> const&,int,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Highlighted_Vertex<VECTOR<double,1> >(VECTOR<double,1> const&,int,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Highlighted_Vertex<VECTOR<double,2> >(VECTOR<double,2> const&,int,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Highlighted_Vertex<VECTOR<double,3> >(VECTOR<double,3> const&,int,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Highlighted_Vertex<VECTOR<float,1> >(VECTOR<float,1> const&,int,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Highlighted_Vertex<VECTOR<float,2> >(VECTOR<float,2> const&,int,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Highlighted_Vertex<VECTOR<float,3> >(VECTOR<float,3> const&,int,OPENGL_COLOR const&);
template void OPENGL_SELECTION::Draw_Vertices_For_Selection<VECTOR<double,2>,1>(SIMPLEX_MESH<1> const&,GEOMETRY_PARTICLES<VECTOR<double,2> > const&);
template void OPENGL_SELECTION::Draw_Vertices_For_Selection<VECTOR<double,2>,2>(SIMPLEX_MESH<2> const&,GEOMETRY_PARTICLES<VECTOR<double,2> > const&);
template void OPENGL_SELECTION::Draw_Vertices_For_Selection<VECTOR<double,3>,1>(SIMPLEX_MESH<1> const&,GEOMETRY_PARTICLES<VECTOR<double,3> > const&);
template void OPENGL_SELECTION::Draw_Vertices_For_Selection<VECTOR<double,3>,2>(SIMPLEX_MESH<2> const&,GEOMETRY_PARTICLES<VECTOR<double,3> > const&);
template void OPENGL_SELECTION::Draw_Vertices_For_Selection<VECTOR<double,3>,3>(SIMPLEX_MESH<3> const&,GEOMETRY_PARTICLES<VECTOR<double,3> > const&);
template void OPENGL_SELECTION::Draw_Vertices_For_Selection<VECTOR<float,2>,1>(SIMPLEX_MESH<1> const&,GEOMETRY_PARTICLES<VECTOR<float,2> > const&);
template void OPENGL_SELECTION::Draw_Vertices_For_Selection<VECTOR<float,2>,2>(SIMPLEX_MESH<2> const&,GEOMETRY_PARTICLES<VECTOR<float,2> > const&);
template void OPENGL_SELECTION::Draw_Vertices_For_Selection<VECTOR<float,3>,1>(SIMPLEX_MESH<1> const&,GEOMETRY_PARTICLES<VECTOR<float,3> > const&);
template void OPENGL_SELECTION::Draw_Vertices_For_Selection<VECTOR<float,3>,2>(SIMPLEX_MESH<2> const&,GEOMETRY_PARTICLES<VECTOR<float,3> > const&);
template void OPENGL_SELECTION::Draw_Vertices_For_Selection<VECTOR<float,3>,3>(SIMPLEX_MESH<3> const&,GEOMETRY_PARTICLES<VECTOR<float,3> > const&);
}

//#####################################################################
// Copyright 2002-2005, Eran Guendelman, Robert Bridson, Sergey Koltakov, Neil Molino, Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace OPENGL_PRIMITIVES
//#####################################################################
#ifndef __OPENGL_PRIMITIVES__
#define __OPENGL_PRIMITIVES__
#include <cstdio>
#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

#ifndef __APPLE__
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#else
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#endif
#include <Tools/Math_Tools/constants.h>
#include <Tools/Matrices/FRAME.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Vectors/VECTOR_2D.h>
#include <Geometry/Basic_Geometry/PLANE.h>
#include <OpenGL/OpenGL/OPENGL_EPS_OUTPUT.h>
#include <OpenGL/OpenGL/OPENGL_POLICY.h>
namespace PhysBAM{

#define WANT_OPENGL_EPS_OUTPUT // Uncomment this if you want support for dumping eps.
extern OPENGL_EPS_OUTPUT<float>* opengl_eps_output;
#ifdef WANT_OPENGL_EPS_OUTPUT
#define IF_OPENGL_EPS_OUTPUT(x) if(opengl_eps_output){x;}
#else
#define IF_OPENGL_EPS_OUTPUT(x)
#endif

inline void glMultMatrix(float* matrix){glMultMatrixf(matrix);}
inline void glMultMatrix(double* matrix){glMultMatrixd(matrix);}

template<int d1>
inline void OpenGL_Draw_Spline(GLenum mode,int resolution,const VECTOR<VECTOR<GLfloat,3>,d1>& control_points)
{
    IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Draw_Bezier(control_points));
    glMap1f(GL_MAP1_VERTEX_3,0.0,1.0,3,d1,&control_points(0)(0));
    glEnable(GL_MAP1_VERTEX_3);
    glMapGrid1d(resolution,0.0,1.0);
    glEvalMesh1(mode,0,resolution);
}
template<int d1>
inline void OpenGL_Draw_Spline(GLenum mode,int resolution,const VECTOR<VECTOR<GLdouble,3>,d1>& control_points)
{
    IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Draw_Bezier(control_points));
    glMap1d(GL_MAP1_VERTEX_3,0.0,1.0,3,d1,&control_points(0)(0));
    glEnable(GL_MAP1_VERTEX_3);
    glMapGrid1d(resolution,0.0,1.0);
    glEvalMesh1(mode,0,resolution);
}
template<class T,int d1,int m>
inline void OpenGL_Draw_Spline(GLenum mode,int resolution,const VECTOR<VECTOR<T,m>,d1>& control_points)
{
    VECTOR<VECTOR<T,3>,d1> pts(control_points);
    OpenGL_Draw_Spline(mode,resolution,pts);
}

inline void OpenGL_Begin(GLenum mode)
{glBegin(mode);IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Begin(mode));}

inline void OpenGL_End()
{glEnd();IF_OPENGL_EPS_OUTPUT(opengl_eps_output->End());}

inline void OpenGL_Eps_Emit(const char* str)
{IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Emit(str));}

inline void OpenGL_Rotate(const ROTATION<VECTOR<float,3> >& r)
{float angle;VECTOR<float,3> axis;r.Get_Angle_Axis(angle,axis);glRotatef(angle*180/(float)pi,axis.x,axis.y,axis.z);}

inline void OpenGL_Rotate(const ROTATION<VECTOR<double,3> >& r)
{double angle;VECTOR<double,3> axis;r.Get_Angle_Axis(angle,axis);glRotated(angle*180/pi,axis.x,axis.y,axis.z);}

inline void OpenGL_Translate(const VECTOR<float,2>& v)
{glTranslatef(v.x,v.y,(float)0);}

inline void OpenGL_Translate(const VECTOR<double,2>& v)
{glTranslated(v.x,v.y,(double)0);}

inline void OpenGL_Translate(const VECTOR<float,3>& v)
{glTranslatef(v.x,v.y,v.z);}

inline void OpenGL_Translate(const VECTOR<double,3>& v)
{glTranslated(v.x,v.y,v.z);}

inline void OpenGL_Scale(const VECTOR<float,3>& v)
{glScalef(v.x,v.y,v.z);}

inline void OpenGL_Scale(const VECTOR<double,3>& v)
{glScaled(v.x,v.y,v.z);}

template<class T>
inline void OpenGL_Transform(const FRAME<VECTOR<T,3> >& frame)
{OpenGL_Translate(frame.t);OpenGL_Rotate(frame.r.Normalized());}

inline void OpenGL_Texture(const VECTOR<float,2>& texture)
{glTexCoord2f(texture[0],texture[1]);}

inline void OpenGL_Texture(const VECTOR<double,2>& texture)
{glTexCoord2d(texture[0],texture[1]);}

inline void OpenGL_Color(const GLfloat* color)
{IF_OPENGL_EPS_OUTPUT((VECTOR<float,3>(color[0],color[1],color[2])));glColor4f(color[0],color[1],color[2],color[3]);}

inline void OpenGL_Color(const GLdouble* color)
{IF_OPENGL_EPS_OUTPUT((VECTOR<double,3>(color[0],color[1],color[2])));glColor4d(color[0],color[1],color[2],color[3]);}

inline void OpenGL_Vertex(const VECTOR<float,3>& v)
{glVertex3f(v.x,v.y,v.z);IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(v));}

inline void OpenGL_Vertex(const VECTOR<float,2>& v)
{glVertex2f(v.x,v.y);IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(v));}

inline void OpenGL_Vertex(const VECTOR<float,1>& v)
{glVertex2f(v.x,(float)0);IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(VECTOR<float,2>(v.x,0)));}

inline void OpenGL_Normal(const VECTOR<float,3>& n)
{glNormal3f(n.x,n.y,n.z);}

inline void OpenGL_RasterPos(const VECTOR<float,1>& v)
{glRasterPos2f(v.x,(float)0);}

inline void OpenGL_RasterPos(const VECTOR<float,2>& v)
{glRasterPos2f(v.x,v.y);}

inline void OpenGL_RasterPos(const VECTOR<float,3>& v)
{glRasterPos3f(v.x,v.y,v.z);}

inline void OpenGL_Vertex(const VECTOR<double,3>& v)
{glVertex3d(v.x,v.y,v.z);IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(v));}

inline void OpenGL_Vertex(const VECTOR<double,2>& v)
{glVertex2d(v.x,v.y);IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(v));}

inline void OpenGL_Vertex(const VECTOR<double,1>& v)
{glVertex2d(v.x,(double)0);IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(VECTOR<double,2>(v.x,0)));}

template<class T,int d>
inline void OpenGL_Line(const VECTOR<T,d>& a,const VECTOR<T,d>& b)
{OpenGL_Vertex(a);OpenGL_Vertex(b);}

template<class T,int d>
inline void OpenGL_Triangle(const VECTOR<T,d>& a, const VECTOR<T,d>& b, const VECTOR<T,d>& c)
{OpenGL_Vertex(a);OpenGL_Vertex(b);OpenGL_Vertex(c);}

template<class T,int d>
inline void OpenGL_Triangle(const VECTOR<VECTOR<T,d>,3>& a)
{OpenGL_Vertex(a(0));OpenGL_Vertex(a(1));OpenGL_Vertex(a(2));}

inline void OpenGL_Normal(const VECTOR<double,3>& n)
{glNormal3d(n.x,n.y,n.z);}

inline void OpenGL_RasterPos(const VECTOR<double,1>& v)
{glRasterPos2f(v.x,(double)0);}

inline void OpenGL_RasterPos(const VECTOR<double,2>& v)
{glRasterPos2d(v.x,v.y);}

inline void OpenGL_RasterPos(const VECTOR<double,3>& v)
{glRasterPos3d(v.x,v.y,v.z);}

inline void OpenGL_LookFrom(const FRAME<VECTOR<double,3> >& frame)
{OpenGL_Rotate(frame.r.Inverse().Normalized());OpenGL_Translate(-frame.t);}

template<class TV>
inline void OpenGL_String(const TV& position,const std::string& str,void* font=GLUT_BITMAP_HELVETICA_12)
{OpenGL_RasterPos(position);
for(unsigned int j=0;j<str.length();j++) glutBitmapCharacter(font,str[j]);}

template<class T>
inline void OpenGL_Quad_2D(const VECTOR<T,2> &bottom_left,const VECTOR<T,2> &top_right)
{OpenGL_Vertex(bottom_left);OpenGL_Vertex(VECTOR<T,2>(bottom_left.x,top_right.y));
 OpenGL_Vertex(top_right);OpenGL_Vertex(VECTOR<T,2>(top_right.x,bottom_left.y));
 IF_OPENGL_EPS_OUTPUT(opengl_eps_output->Vertex(VECTOR<T,2>(bottom_left.x,bottom_left.y));opengl_eps_output->Vertex(VECTOR<T,2>(bottom_left.x,top_right.y));
     opengl_eps_output->Vertex(VECTOR<T,2>(top_right.x,top_right.y));opengl_eps_output->Vertex(VECTOR<T,2>(top_right.x,bottom_left.y)));}

template<class T>
inline void OpenGL_Triangle_Strip_2D(const VECTOR<T,2> &bottom_left,const VECTOR<T,2> &top_right)
{OpenGL_Vertex(bottom_left);OpenGL_Vertex(VECTOR<T,2>(bottom_left.x,top_right.y));
 OpenGL_Vertex(VECTOR<T,2>(top_right.x,bottom_left.y));OpenGL_Vertex(top_right);}

template<class T>
inline void OpenGL_Quad(const VECTOR<T,3> &bottom_left,const VECTOR<T,3> &right,const VECTOR<T,3>& up)
{OpenGL_Vertex(bottom_left);OpenGL_Vertex(bottom_left+up);OpenGL_Vertex(bottom_left+up+right);OpenGL_Vertex(bottom_left+right);}

template<class T>
inline void OpenGL_Clip_Plane(GLenum id,const PLANE<T> &plane)
{GLdouble equation[4]={plane.normal.x,plane.normal.y,plane.normal.z,-VECTOR<T,3>::Dot_Product(plane.normal,plane.x0)};
glClipPlane(id,equation);}
}
#endif

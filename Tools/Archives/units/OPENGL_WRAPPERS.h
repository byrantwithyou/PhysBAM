//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header OPENGL_WRAPPERS
//#####################################################################
#ifndef __OPENGL_WRAPPERS__
#define __OPENGL_WRAPPERS__

#include "QUANTITY.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
namespace PhysBAM{
namespace UNITS{

inline void
glRotatef(const QUANTITY<float>& angle,const QUANTITY<float>& x,const QUANTITY<float>& y,const QUANTITY<float>& z)
{::glRotatef(angle.value,x.value,y.value,z.value);}

inline void
glRotated(const QUANTITY<double>& angle,const QUANTITY<double>& x,const QUANTITY<double>& y,const QUANTITY<double>& z)
{::glRotated(angle.value,x.value,y.value,z.value);}

inline void
glTranslatef(const QUANTITY<float>& x,const QUANTITY<float>& y,const QUANTITY<float>& z)
{::glTranslatef(x.value,y.value,z.value);}

inline void
glTranslated(const QUANTITY<double>& x,const QUANTITY<double>& y,const QUANTITY<double>& z)
{::glTranslated(x.value,y.value,z.value);}

inline void
glScalef(const QUANTITY<float>& x,const QUANTITY<float>& y,const QUANTITY<float>& z)
{::glScalef(x.value,y.value,z.value);}

inline void
glScaled(const QUANTITY<double>& x,const QUANTITY<double>& y,const QUANTITY<double>& z)
{::glScaled(x.value,y.value,z.value);}

inline void
glVertex2f(const QUANTITY<float>& x,const QUANTITY<float>& y)
{::glVertex2f(x.value,y.value);}

inline void
glVertex2d(const QUANTITY<double>& x,const QUANTITY<double>& y)
{::glVertex2d(x.value,y.value);}

inline void
glVertex3f(const QUANTITY<float>& x,const QUANTITY<float>& y,const QUANTITY<float>& z)
{::glVertex3f(x.value,y.value,z.value);}

inline void
glVertex3d(const QUANTITY<double>& x,const QUANTITY<double>& y,const QUANTITY<double>& z)
{::glVertex3d(x.value,y.value,z.value);}

inline void
glTexCoord2f(const QUANTITY<float>& x,const QUANTITY<float>& y)
{::glTexCoord2f(x.value,y.value);}

inline void
glTexCoord2d(const QUANTITY<double>& x,const QUANTITY<double>& y)
{::glTexCoord2d(x.value,y.value);}

inline void
glRasterPos2f(const QUANTITY<float>& x,const QUANTITY<float>& y)
{::glRasterPos2f(x.value,y.value);}

inline void
glRasterPos2d(const QUANTITY<double>& x,const QUANTITY<double>& y)
{::glRasterPos2d(x.value,y.value);}

inline void
glRasterPos3f(const QUANTITY<float>& x,const QUANTITY<float>& y,const QUANTITY<float>& z)
{::glRasterPos3f(x.value,y.value,z.value);}

inline void
glRasterPos3d(const QUANTITY<double>& x,const QUANTITY<double>& y,const QUANTITY<double>& z)
{::glRasterPos3d(x.value,y.value,z.value);}

inline void
glNormal3f(const QUANTITY<float>& x,const QUANTITY<float>& y,const QUANTITY<float>& z)
{::glNormal3f(x.value,y.value,z.value);}

inline void
glNormal3d(const QUANTITY<double>& x,const QUANTITY<double>& y,const QUANTITY<double>& z)
{::glNormal3d(x.value,y.value,z.value);}

inline void
glColor3f(const QUANTITY<float>& x,const QUANTITY<double>& y,const QUANTITY<double>& z)
{::glColor3f(x.value,y.value,z.value);}

inline void
glColor4fv(const QUANTITY<float> c[4])
{float values[4]={c[0].value,c[1].value,c[2].value,c[3].value};::glColor4fv(values);}

inline void
glMaterialf(GLenum a,GLenum b,const QUANTITY<float>& c)
{::glMaterialf(a,b,c.value);}

inline void
glMaterialfv(GLenum a,GLenum b,const QUANTITY<float> c[4])
{float values[4]={c[0].value,c[1].value,c[2].value,c[3].value};::glMaterialfv(a,b,values);}

inline void
glLightfv(GLenum a,GLenum b,const QUANTITY<float> c[4])
{float values[4]={c[0].value,c[1].value,c[2].value,c[3].value};::glLightfv(a,b,values);}

inline void
glLightModelfv(GLenum a,const QUANTITY<float> c[4])
{float values[4]={c[0].value,c[1].value,c[2].value,c[3].value};::glLightModelfv(a,values);}

inline void
glPointSize(const QUANTITY<float>& x)
{::glPointSize(x.value);}

inline void
glLineWidth(const QUANTITY<float>& x)
{::glLineWidth(x.value);}

inline void
glGetFloatv(GLenum a,QUANTITY<float>* p)
{float values[16],sentinel=FLT_MAX;for(int i=0;i<16;i++) values[i]=sentinel;
::glGetFloatv(a,values);for(int i=0;i<16;i++){if(p[i]==sentinel) break;p[i]=values[i];}}

inline void
glGetDoublev(GLenum a,QUANTITY<double>* p)
{double values[16],sentinel=FLT_MAX;for(int i=0;i<16;i++) values[i]=sentinel;
::glGetDoublev(a,values);for(int i=0;i<16;i++){if(p[i]==sentinel) break;p[i]=values[i];}}

inline void
glMultMatrixf(const QUANTITY<float> p[16])
{float values[16];for(int i=0;i<16;i++) values[i]=p[i].value;::glMultMatrixf(values);}

inline void
gluDisk(GLUquadric* quad,const QUANTITY<double>& inner,const QUANTITY<double>& outer,GLint slices,GLint loops)
{::gluDisk(quad,inner.value,outer.value,slices,loops);}

inline void
gluCylinder(GLUquadric* quad,const QUANTITY<double>& base,const QUANTITY<double>& top,const QUANTITY<double>& height,GLint slices,GLint stacks)
{::gluCylinder(quad,base.value,top.value,height.value,slices,stacks);}

inline void
gluPerspective(const QUANTITY<float>& fovy,const QUANTITY<double>& aspect,const QUANTITY<double>& zNear,const QUANTITY<double>& zFar)
{::gluPerspective(fovy.value,aspect.value,zNear.value,zFar.value);}

inline void
gluPickMatrix(const QUANTITY<double>& x,const QUANTITY<double>& y,const QUANTITY<double>& delX,const QUANTITY<double>& delY,GLint* viewport)
{::gluPickMatrix(x.value,y.value,delX.value,delY.value,viewport);}

inline void
gluProject(const QUANTITY<double>& objX,const QUANTITY<double>& objY,const QUANTITY<double>& objZ,const QUANTITY<double>* model,const QUANTITY<double>* proj,const GLint* view,
    QUANTITY<double>* winX,QUANTITY<double>* winY,QUANTITY<double>* winZ)
{double model_values[16],proj_values[16];for(int i=0;i<16;i++){model_values[i]=model[i].value;proj_values[i]=proj[i].value;}
::gluProject(objX.value,objY.value,objZ.value,model_values,proj_values,view,&winX->value,&winY->value,&winZ->value);}

inline void
gluUnProject(const QUANTITY<double>& winX,const QUANTITY<double>& winY,const QUANTITY<double>& winZ,const QUANTITY<double>* model,const QUANTITY<double>* proj,const GLint* view,
    QUANTITY<double>* objX,QUANTITY<double>* objY,QUANTITY<double>* objZ)
{double model_values[16],proj_values[16];for(int i=0;i<16;i++){model_values[i]=model[i].value;proj_values[i]=proj[i].value;}
::gluUnProject(winX.value,winY.value,winZ.value,model_values,proj_values,view,&objX->value,&objY->value,&objZ->value);}

inline void
glClipPlane(GLenum a,const QUANTITY<double> equation[4])
{double values[4];for(int i=0;i<4;i++) values[i]=equation[4].value;
::glClipPlane(a,values);}

inline void
glutWireCube(const QUANTITY<double>& a)
{::glutWireCube(a.value);}

}
}
#endif

//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_EPS_OUTPUT.h>
#ifdef _WIN32
#include <windows.h>
#endif

#ifndef __APPLE__
#ifndef USE_OPENGLES
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#else
#include <GLES/gl.h>
#endif
#else
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#endif
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_EPS_OUTPUT<T>::
OPENGL_EPS_OUTPUT(const std::string& filename)
    :stream(*FILE_UTILITIES::Safe_Open_Output(filename,false,false))
{
    Head();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_EPS_OUTPUT<T>::
~OPENGL_EPS_OUTPUT()
{
    Tail();
    delete &stream;
}
//#####################################################################
// Function Begin
//#####################################################################
template<class T> void OPENGL_EPS_OUTPUT<T>::
Begin(int mode)
{
    glmode=mode;
}
//#####################################################################
// Function End
//#####################################################################
template<class T> void OPENGL_EPS_OUTPUT<T>::
End()
{
#ifndef USE_OPENGLES
    Transform_Buffer();
    switch(glmode){
        case GL_POINTS: for(int i=0;i<buffer.m;i++) Draw_Point(buffer(i)); break;
        case GL_LINES: for(int i=0;i<buffer.m;i+=2) Draw_Line(buffer(i),buffer(i+1)); break;
        case GL_TRIANGLES: for(int i=0;i<buffer.m;i+=3) Draw_Polygon(i,3); break;
        case GL_QUADS: for(int i=0;i<buffer.m;i+=4) Draw_Polygon(i,4); break;
        case GL_POLYGON: Draw_Polygon(1,buffer.m); break;
        case GL_LINE_LOOP: for(int i=1;i<buffer.m;i++) Draw_Line(buffer.Last(),buffer(1)); break;
        default:LOG::cout<<"Unhandled mode "<<glmode<<std::endl;break;
    }
    buffer.Remove_All();
#endif
}
//#####################################################################
// Function Vertex
//#####################################################################
template<class T> template<class T2> void OPENGL_EPS_OUTPUT<T>::
Vertex(const VECTOR<T2,3>& p)
{
    buffer.Append(VECTOR<T,3>(p));
}
//#####################################################################
// Function Vertex
//#####################################################################
template<class T> template<class T2> void OPENGL_EPS_OUTPUT<T>::
Vertex(const VECTOR<T2,2>& p)
{
    buffer.Append(VECTOR<T,3>(VECTOR<T,2>(p)));
}
//#####################################################################
// Function Vertex
//#####################################################################
template<class T> template<class T2> void OPENGL_EPS_OUTPUT<T>::
Vertex(const VECTOR<T2,1>& p)
{
    buffer.Append(VECTOR<T,3>(p.x,0,0));
}
//#####################################################################
// Function Head
//#####################################################################
template<class T> void OPENGL_EPS_OUTPUT<T>::
Head()
{
    float view[4];
    glGetFloatv(GL_VIEWPORT,view);

    Emit("%!PS-Adobe-3.0 EPSF-3.0\n");
    Emit("%%BoundingBox: ");
    stream<<view[0]<<" "<<view[1]<<" "<<(view[2]+view[0])<<" "<<(view[3]+view[1])<<"\n";
    Emit("0 setlinewidth\n");
    Emit("/pointradius 1 def\n");
}
//#####################################################################
// Function Tail
//#####################################################################
template<class T> void OPENGL_EPS_OUTPUT<T>::
Tail()
{
}
//#####################################################################
// Function Emit
//#####################################################################
template<class T> void OPENGL_EPS_OUTPUT<T>::
Emit(const char* s)
{
    stream<<s;
}
//#####################################################################
// Function Emit
//#####################################################################
template<class T> void OPENGL_EPS_OUTPUT<T>::
Emit(const TV& p)
{
    stream<<p.x<<" "<<p.y<<" ";
}
//#####################################################################
// Function Draw_Point
//#####################################################################
template<class T> void OPENGL_EPS_OUTPUT<T>::
Draw_Point(const TV& p)
{
    Emit(p);
    Emit("pointradius 0 360 arc fill stroke\n");
}
//#####################################################################
// Function Draw_Line
//#####################################################################
template<class T> void OPENGL_EPS_OUTPUT<T>::
Draw_Line(const TV& a,const TV& b)
{
    Emit(a);
    Emit("moveto ");
    Emit(b);
    Emit("lineto stroke\n");
}
//#####################################################################
// Function Draw_Polygon
//#####################################################################
template<class T> void OPENGL_EPS_OUTPUT<T>::
Draw_Polygon(int i,int n)
{
#ifndef USE_OPENGLES
    GLint mode[2]={0};
    glGetIntegerv(GL_POLYGON_MODE,mode);
    const char* str = "closepath stroke\n";
    if(mode[0] == GL_FILL) str = "closepath fill stroke\n";
    Emit(buffer(i));
    Emit("moveto ");
    for(int j=i+1;j<i+n;j++){
        Emit(buffer(j));
        Emit("lineto ");}
    Emit(str);
#endif
}
//#####################################################################
// Function Set_Color
//#####################################################################
template<class T> void OPENGL_EPS_OUTPUT<T>::
Set_Color(const TV& color)
{
    stream<<color.x<<" "<<color.y<<" "<<color.z<<" setrgbcolor"<<std::endl;
}
//#####################################################################
// Function Transform_Buffer
//#####################################################################
template<class T> void OPENGL_EPS_OUTPUT<T>::
Transform_Buffer()
{
#ifndef USE_OPENGLES
    MATRIX<T,4> proj,model;
    T view[4];

    if(sizeof(T)==sizeof(float)){
        glGetFloatv(GL_PROJECTION_MATRIX,(GLfloat*)proj.x);
        glGetFloatv(GL_MODELVIEW_MATRIX,(GLfloat*)model.x);
        glGetFloatv(GL_VIEWPORT,(GLfloat*)view);}
    else{
        glGetDoublev(GL_PROJECTION_MATRIX,(GLdouble*)proj.x);
        glGetDoublev(GL_MODELVIEW_MATRIX,(GLdouble*)model.x);
        glGetDoublev(GL_VIEWPORT,(GLdouble*)view);}

    for(int i=0;i<buffer.m;i++){
        VECTOR<T,4> obj=buffer(i).Append(1);
        VECTOR<T,4> eye=model*obj;
        VECTOR<T,4> clip=proj*eye;
        TV device=clip.Remove_Index(3)/clip(3);
        TV window((device.x+1)*(view[2]/2)+view[0],(device.y+1)*(view[3]/2)+view[1],device.z);
        buffer(i)=window;}
#endif
}
//#####################################################################
// Function Draw_Arrays
//#####################################################################
template<class T> void OPENGL_EPS_OUTPUT<T>::
Draw_Arrays(int mode,int dimension,int length,const void* vertices)
{
    Begin(mode);
    float* p=(float*)vertices;
    if(dimension==1) for(int i=0;i<length;i++) buffer.Append(VECTOR<T,3>(p[i],0,0));
    else if(dimension==2) for(int i=0;i<2*length;i+=2) buffer.Append(VECTOR<T,3>(p[i],p[i+1],0));
    else if(dimension==3) for(int i=0;i<3*length;i+=3) buffer.Append(VECTOR<T,3>(p[i],p[i+1],p[i+2]));
    End();
}
template class OPENGL_EPS_OUTPUT<float>;
template void OPENGL_EPS_OUTPUT<float>::Vertex<float>(VECTOR<float,2> const&);
template void OPENGL_EPS_OUTPUT<float>::Vertex<float>(VECTOR<float,3> const&);
template void OPENGL_EPS_OUTPUT<float>::Vertex<double>(VECTOR<double,2> const&);
template void OPENGL_EPS_OUTPUT<float>::Vertex<double>(VECTOR<double,3> const&);
template void OPENGL_EPS_OUTPUT<float>::Vertex<float>(VECTOR<float,1> const&);
template void OPENGL_EPS_OUTPUT<float>::Vertex<double>(VECTOR<double,1> const&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_EPS_OUTPUT<double>;
template void OPENGL_EPS_OUTPUT<double>::Vertex<float>(VECTOR<float,2> const&);
template void OPENGL_EPS_OUTPUT<double>::Vertex<float>(VECTOR<float,3> const&);
template void OPENGL_EPS_OUTPUT<double>::Vertex<double>(VECTOR<double,2> const&);
template void OPENGL_EPS_OUTPUT<double>::Vertex<double>(VECTOR<double,3> const&);
template void OPENGL_EPS_OUTPUT<double>::Vertex<double>(VECTOR<double,1> const&);
template void OPENGL_EPS_OUTPUT<double>::Vertex<float>(VECTOR<float,1> const&);
#endif

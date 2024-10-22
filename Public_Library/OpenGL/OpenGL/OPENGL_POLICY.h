//#####################################################################
// Copyright 2011, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_POLICY 
//#####################################################################
#ifndef __OPENGL_POLICY__
#define __OPENGL_POLICY__

#ifndef __APPLE__
#include <GL/gl.h>
#else
#include <OpenGL/gl.h>
#endif

namespace PhysBAM{

template<class T> struct OPENGL_POLICY;

//#####################################################################
// float
//#####################################################################
template<> struct OPENGL_POLICY<float>
{
    typedef GLfloat T_GL;
};
//#####################################################################
// double
//#####################################################################
template<> struct OPENGL_POLICY<double>
{
    typedef GLdouble T_GL;
};
//#####################################################################
}
#endif

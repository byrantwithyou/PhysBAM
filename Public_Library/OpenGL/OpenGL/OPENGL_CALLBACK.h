//#####################################################################
// Copyright 2002-2005, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_CALLBACK
//##################################################################### 
#ifndef __OPENGL_CALLBACK__
#define __OPENGL_CALLBACK__
#include <functional>
#include <iostream>
#include <string>
namespace PhysBAM{

struct OPENGL_CALLBACK
{
    std::function<void()> func;
    const char* help;
};
}
#endif


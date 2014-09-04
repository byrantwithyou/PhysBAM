//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __OPENGL_WINDOW_PBUFFER__
#define __OPENGL_WINDOW_PBUFFER__

#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <OpenGL/OpenGL/OPENGL_WINDOW.h>
namespace PhysBAM{
class OPENGL_PBUFFER;

template<class T>
class OPENGL_WINDOW_PBUFFER:public OPENGL_WINDOW<T>
{
    using OPENGL_WINDOW<T>::opengl_world;
    int width,height;
    OPENGL_PBUFFER* pbuffer;

//#####################################################################
public:
    OPENGL_WINDOW_PBUFFER(OPENGL_WORLD<T>& world_input,const std::string& window_title_input,const int width_input,const int height_input);
    virtual ~OPENGL_WINDOW_PBUFFER();
    void Setup_Idle(const bool use) PHYSBAM_OVERRIDE;
    void Setup_Timer(const float wait_milliseconds) PHYSBAM_OVERRIDE;
    void Redisplay() PHYSBAM_OVERRIDE;
    void Main_Loop() PHYSBAM_OVERRIDE;
    void Request_Resize(const int width,const int height) PHYSBAM_OVERRIDE;
    void Request_Move(const int x,const int y) PHYSBAM_OVERRIDE;
    int Width() const PHYSBAM_OVERRIDE;
    int Height() const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

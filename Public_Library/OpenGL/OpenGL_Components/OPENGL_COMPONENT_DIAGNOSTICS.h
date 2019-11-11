//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __OPENGL_COMPONENT_DIAGNOSTICS__
#define __OPENGL_COMPONENT_DIAGNOSTICS__

#include <Core/Arrays/ARRAY.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{

template<class T>
class OPENGL_COMPONENT_DIAGNOSTICS:public OPENGL_COMPONENT<T>
{
    std::string filename;
    bool valid;
    ARRAY<std::string> lines;

//#####################################################################
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_COMPONENT<T>::component_name;using OPENGL_COMPONENT<T>::viewer_dir;
    
    OPENGL_COMPONENT_DIAGNOSTICS(const VIEWER_DIR& viewer_dir,const std::string& filename);
    virtual ~OPENGL_COMPONENT_DIAGNOSTICS();
private:
    void Reinitialize();
    void Print_Selection_Info(std::ostream& ostream) const override;    
    void Set_Frame() override;
    //virtual RANGE<TV> Bounding_Box() const override;
    bool Use_Bounding_Box() const override;
    void Display() const override;
//#####################################################################
};
}
#endif

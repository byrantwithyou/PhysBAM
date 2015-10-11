//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __OPENGL_COMPONENT_DIAGNOSTICS__
#define __OPENGL_COMPONENT_DIAGNOSTICS__

#include <Tools/Arrays/ARRAY.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{

template<class T>
class OPENGL_COMPONENT_DIAGNOSTICS:public OPENGL_COMPONENT<T>
{
    std::string filename;
    int frame_loaded;
    bool valid;
    ARRAY<std::string> lines;

//#####################################################################
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;using OPENGL_COMPONENT<T>::component_name;
    using OPENGL_COMPONENT<T>::is_animation;
    OPENGL_COMPONENT_DIAGNOSTICS(STREAM_TYPE stream_type,const std::string& filename);
    virtual ~OPENGL_COMPONENT_DIAGNOSTICS();
private:
    void Reinitialize();
    void Print_Selection_Info(std::ostream& ostream,OPENGL_SELECTION<T>* selection) const override;    
    bool Valid_Frame(int frame) const override;
    void Set_Frame(int frame) override;
    //virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;
    bool Use_Bounding_Box() const override;
    void Display() const override;
//#####################################################################
};
}
#endif

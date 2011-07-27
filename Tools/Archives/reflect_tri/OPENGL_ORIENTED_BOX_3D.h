//#####################################################################
// Class OPENGL_ORIENTED_BOX_3D
//##################################################################### 
#ifndef __OPENGL_ORIENTED_BOX_3D__
#define __OPENGL_ORIENTED_BOX_3D__

#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>

using namespace PhysBAM;

template<class T>
class OPENGL_ORIENTED_BOX_3D : public OPENGL_OBJECT
{
public: 
    ORIENTED_BOX_3D<T> *box;

    OPENGL_ORIENTED_BOX_3D(ORIENTED_BOX_3D<T> *box_input):box(box_input)
    {}
    ~OPENGL_ORIENTED_BOX_3D()
    {}

    virtual void Display(const int in_color=1) const;
    virtual BOX_3D<float> Bounding_Box() const;
};

#endif

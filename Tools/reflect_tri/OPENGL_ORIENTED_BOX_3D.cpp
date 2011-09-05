#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include "OPENGL_ORIENTED_BOX_3D.h"

template class OPENGL_ORIENTED_BOX_3D<float>;
template class OPENGL_ORIENTED_BOX_3D<double>;

template<class T> void OPENGL_ORIENTED_BOX_3D<T>::
Display(const int in_color) const
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();

    glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT);
    glDisable(GL_LIGHTING);
    glLineWidth(3.0f);

    VECTOR_3D<T> corner2=box->corner+box->edge1;
    glBegin(GL_LINE_LOOP);
    OpenGL_Vertex(box->corner); OpenGL_Vertex(box->corner+box->edge2); 
    OpenGL_Vertex(box->corner+box->edge2+box->edge3); OpenGL_Vertex(box->corner+box->edge3);
    glEnd();
    glBegin(GL_LINE_LOOP);
    OpenGL_Vertex(corner2); OpenGL_Vertex(corner2+box->edge2); 
    OpenGL_Vertex(corner2+box->edge2+box->edge3); OpenGL_Vertex(corner2+box->edge3);
    glEnd();
    glBegin(GL_LINES);
    OpenGL_Vertex(box->corner); OpenGL_Vertex(corner2);
    OpenGL_Vertex(box->corner+box->edge2); OpenGL_Vertex(corner2+box->edge2);
    OpenGL_Vertex(box->corner+box->edge2+box->edge3); OpenGL_Vertex(corner2+box->edge2+box->edge3);
    OpenGL_Vertex(box->corner+box->edge3); OpenGL_Vertex(corner2+box->edge3);
    glEnd();

    glPopAttrib();
    glPopMatrix();
}

template<class T> BOX_3D<float> OPENGL_ORIENTED_BOX_3D<T>::
Bounding_Box() const
{
    return box->Axis_Aligned_Bounding_Box();
}
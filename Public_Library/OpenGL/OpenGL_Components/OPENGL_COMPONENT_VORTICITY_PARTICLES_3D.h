//#####################################################################
// Copyright 2004, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_VORTICITY_PARTICLES_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_VORTICITY_PARTICLES_3D__
#define __OPENGL_COMPONENT_VORTICITY_PARTICLES_3D__

#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_3D.h>

namespace PhysBAM{

template<class T>
class OPENGL_COMPONENT_VORTICITY_PARTICLES_3D:public OPENGL_COMPONENT_PARTICLES_3D<T>
{
    typedef VECTOR<T,3> TV;
    using OPENGL_COMPONENT_PARTICLES_3D<T>::have_velocities;
public:
    using OPENGL_COMPONENT_PARTICLES_3D<T>::particles;using OPENGL_COMPONENT_PARTICLES_3D<T>::opengl_vector_field;

    OPENGL_COMPONENT_VORTICITY_PARTICLES_3D(STREAM_TYPE stream_type,const std::string &filename,bool use_ids_input=true);
private:
    void Reinitialize(bool force=false) override;
//#####################################################################
};
}
#endif

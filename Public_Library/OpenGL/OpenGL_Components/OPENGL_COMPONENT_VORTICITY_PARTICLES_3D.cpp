//#####################################################################
// Copyright 2004-2005, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_VORTICITY_PARTICLES_3D
//#####################################################################
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_VORTICITY_PARTICLES_3D.h>
using namespace PhysBAM;

const ATTRIBUTE_ID ATTRIBUTE_ID_VORTICITY(16);
template<class T> OPENGL_COMPONENT_VORTICITY_PARTICLES_3D<T>::
OPENGL_COMPONENT_VORTICITY_PARTICLES_3D(STREAM_TYPE stream_type,const std::string &filename,bool use_ids_input)
    :OPENGL_COMPONENT_PARTICLES_3D<T>(stream_type,filename,"",use_ids_input,false)
{}

template<class T> void OPENGL_COMPONENT_VORTICITY_PARTICLES_3D<T>::
Reinitialize(bool force)
{
    const ARRAY_VIEW<typename TV::SPIN>& vorticity=*particles->template Get_Array<typename TV::SPIN>(ATTRIBUTE_ID_VORTICITY);
    OPENGL_COMPONENT_PARTICLES_3D<T>::Reinitialize(force);
    if(!have_velocities){
        have_velocities=true;
        opengl_vector_field.vector_field.Resize(particles->Size());
        int idx=0;
        for(int i=0;i<particles->Size();i++)
            opengl_vector_field.vector_field(idx++)=vorticity(i);}
}
namespace PhysBAM{
template class OPENGL_COMPONENT_VORTICITY_PARTICLES_3D<float>;
template class OPENGL_COMPONENT_VORTICITY_PARTICLES_3D<double>;
}

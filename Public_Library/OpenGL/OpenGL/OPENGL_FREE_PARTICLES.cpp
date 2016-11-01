//#####################################################################
// Copyright 2007-2008, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <OpenGL/OpenGL/OPENGL_FREE_PARTICLES.h>
using namespace PhysBAM;
//##################################################################### 
// Constructor
//##################################################################### 
template<class TV> OPENGL_FREE_PARTICLES<TV>::
OPENGL_FREE_PARTICLES(STREAM_TYPE stream_type,DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection,INDIRECT_ARRAY<ARRAY_VIEW<TV> >& points,const OPENGL_COLOR& color,const T point_size)
    :BASE(stream_type,points,color,point_size),deformable_body_collection(deformable_body_collection)
{}
//##################################################################### 
// Function Print_Selection_Info
//##################################################################### 
template<class TV> void OPENGL_FREE_PARTICLES<TV>::
Print_Selection_Info(std::ostream &output_stream) const
{
    deformable_body_collection.particles.Print(output_stream,selected_index);
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_FREE_PARTICLES<VECTOR<float,2> >;
template class OPENGL_FREE_PARTICLES<VECTOR<float,3> >;
template class OPENGL_FREE_PARTICLES<VECTOR<double,2> >;
template class OPENGL_FREE_PARTICLES<VECTOR<double,3> >;
}

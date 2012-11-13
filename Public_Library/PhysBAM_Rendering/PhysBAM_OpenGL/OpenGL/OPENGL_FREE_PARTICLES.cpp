//#####################################################################
// Copyright 2007-2008, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_FREE_PARTICLES.h>
using namespace PhysBAM;
//##################################################################### 
// Constructor
//##################################################################### 
template<class TV> OPENGL_FREE_PARTICLES<TV>::
OPENGL_FREE_PARTICLES(DEFORMABLE_GEOMETRY_COLLECTION<TV>& deformable_geometry_collection,INDIRECT_ARRAY<ARRAY_VIEW<TV> >& points,const OPENGL_COLOR& color,const T point_size)
    :BASE(points,color,point_size),deformable_geometry_collection(deformable_geometry_collection)
{}
//##################################################################### 
// Function Print_Selection_Info
//##################################################################### 
template<class TV> void OPENGL_FREE_PARTICLES<TV>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION* selection) const
{
    BASE::Print_Selection_Info(output_stream,selection);
    if(selection->type!=OPENGL_SELECTION::POINTS_3D) return;
    int particle_index=points.indices(dynamic_cast<OPENGL_SELECTION_POINTS_3D<T>*>(selection)->index);
    deformable_geometry_collection.particles.Print(output_stream,particle_index);
}
//#####################################################################
template class OPENGL_FREE_PARTICLES<VECTOR<float,2> >;
template class OPENGL_FREE_PARTICLES<VECTOR<float,3> >;
template class OPENGL_FREE_PARTICLES<VECTOR<double,2> >;
template class OPENGL_FREE_PARTICLES<VECTOR<double,3> >;

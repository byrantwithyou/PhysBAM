//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_UNIFORM_SLICE
//##################################################################### 
#include <OpenGL/OpenGL/OPENGL_UNIFORM_SLICE.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
using namespace PhysBAM;
//#####################################################################
// Function Print_Slice_Info
//#####################################################################
template<class T> void OPENGL_UNIFORM_SLICE<T>::
Print_Slice_Info(std::ostream& output_stream)
{
    const char* axis_name[]={"","x","y","z"};
    if(mode==CELL_SLICE){output_stream<<"Slice: cell="<<index<<", "<<axis_name[axis]<<"="<<grid.Center(TV_INT()+index)[axis]<<std::endl;}
    else if(mode==NODE_SLICE){output_stream<<"Slice: node="<<index<<", "<<axis_name[axis]<<"="<<grid.Node(TV_INT()+index)[axis]<<std::endl;}
}
//#####################################################################
// Function Update_Clip_Planes
//#####################################################################
template<class T> void OPENGL_UNIFORM_SLICE<T>::
Update_Clip_Planes()
{
    if(mode==NO_SLICE) {
        if(clip_plane_id1!=0) world.Remove_Clipping_Plane(clip_plane_id1); 
        if(clip_plane_id2!=0) world.Remove_Clipping_Plane(clip_plane_id2);
        clip_plane_id1=clip_plane_id2=0;
    }
    else {
        T pos=(mode==NODE_SLICE)?grid.Node(TV_INT()+index)[axis]:grid.Center(TV_INT()+index)[axis];
        PLANE<T> plane1(VECTOR<T,3>(0,0,0),VECTOR<T,3>(0,0,0)),plane2(VECTOR<T,3>(0,0,0),VECTOR<T,3>(0,0,0));
        plane1.normal[axis]=1;plane1.x0[axis]=pos-grid.dX[axis]/1.9;
        plane2.normal[axis]=-1;plane2.x0[axis]=pos+grid.dX[axis]/1.9;
        if(clip_plane_id1==0) clip_plane_id1=world.Add_Clipping_Plane(plane1);
        else world.Set_Clipping_Plane(clip_plane_id1,plane1);
        if(clip_plane_id2==0) clip_plane_id2=world.Add_Clipping_Plane(plane2);
        else world.Set_Clipping_Plane(clip_plane_id2,plane2);
    }
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_UNIFORM_SLICE<double>;
template class OPENGL_UNIFORM_SLICE<float>;
}

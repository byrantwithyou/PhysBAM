//#####################################################################
// Copyright 2004, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Utilities/NONCOPYABLE.h>
#include <Ray_Tracing/Rendering/RENDERING_RAY.h>
#include <Ray_Tracing/Rendering/RENDERING_RAY_DEBUG.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> RENDERING_RAY_DEBUG<T>::
RENDERING_RAY_DEBUG()
    :parent(0),hit_object(false)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> RENDERING_RAY_DEBUG<T>::
RENDERING_RAY_DEBUG(const RENDERING_RAY<T>& ray_input)
    :ray(ray_input),parent(0),hit_object(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> RENDERING_RAY_DEBUG<T>::
~RENDERING_RAY_DEBUG()
{
    for(int i=0;i<children.m;i++)delete children(i);
}
//#####################################################################
// Function Add_Child
//#####################################################################
template<class T> void RENDERING_RAY_DEBUG<T>::
Add_Child(RENDERING_RAY<T>& ray_to_add)
{
    RENDERING_RAY_DEBUG* child=new RENDERING_RAY_DEBUG(ray_to_add);
    child->parent=this;
    ray_to_add.debug_ray=child->ray.debug_ray=child;
    children.Append(child);
}
//#####################################################################
// Function Add_Comment
//#####################################################################
template<class T> void RENDERING_RAY_DEBUG<T>::
Add_Comment(const std::string& comment_string)
{
    comments.Append(comment_string);
}
//#####################################################################
// Function Print
//#####################################################################
template<class T> void RENDERING_RAY_DEBUG<T>::
Print(std::ostream& output,const int spaces)
{
    for(int space=0;space<spaces;space++)output<<" ";
    std::string ray_type_str;
    switch(ray.ray_type){
        case RENDERING_RAY<T>::PHOTON_RAY:ray_type_str="PHOTON";break;
        case RENDERING_RAY<T>::PHOTON_GATHER_RAY:ray_type_str="PHOTON_GATHER";break;
        case RENDERING_RAY<T>::DUMMY_RAY:ray_type_str="DUMMY";break;
        case RENDERING_RAY<T>::COLOR_RAY:ray_type_str="COLOR";break;
        case RENDERING_RAY<T>::SHADOW_RAY:ray_type_str="SHADOW";break;
        case RENDERING_RAY<T>::UNKNOWN_RAY:ray_type_str="UNKNOWN";break;}
    output<<"RAY "<<ray_type_str<<" endpoint: "<<ray.ray.endpoint<<" direction: "<<ray.ray.direction<<std::endl;
    for(int i=0;i<children.m;i++)
        children(i)->Print(output,spaces+5);
}
template class RENDERING_RAY_DEBUG<double>;
template class RENDERING_RAY_DEBUG<float>;
}

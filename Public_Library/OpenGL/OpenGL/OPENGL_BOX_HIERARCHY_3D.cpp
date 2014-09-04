//#####################################################################
// Copyright 2004-2005, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_BOX_HIERARCHY_3D
//##################################################################### 
#include <Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
#include <OpenGL/OpenGL/OPENGL_BOX_3D.h>
#include <OpenGL/OpenGL/OPENGL_BOX_HIERARCHY_3D.h>
#include <OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
using namespace PhysBAM;
//#####################################################################
// Function Display
//##################################################################### 
template<class T> void OPENGL_BOX_HIERARCHY_3D<T>::
Display() const
{
    if(hierarchy) Display_Helper(hierarchy->root,1);
}
//#####################################################################
// Function Display_Helper
//##################################################################### 
template<class T> void OPENGL_BOX_HIERARCHY_3D<T>::
Display_Helper(const int cell,const int height) const
{
    if(height>=min_height && height<=max_height){
        OPENGL_BOX_3D<T> opengl_box(hierarchy->box_hierarchy(cell),color);
        opengl_box.Display();}
    if(!hierarchy->Leaf(cell)){
        int left_cell,right_cell;
        hierarchy->children(cell-hierarchy->leaves).Get(left_cell,right_cell);
        Display_Helper(left_cell,height+1);
        Display_Helper(right_cell,height+1);}
}
//#####################################################################
// Function Bounding_Box
//##################################################################### 
template<class T> RANGE<VECTOR<T,3> > OPENGL_BOX_HIERARCHY_3D<T>::
Bounding_Box() const
{
    if(hierarchy) return World_Space_Box(hierarchy->box_hierarchy(hierarchy->root));
    return RANGE<TV>(TV(),TV());
}
//#####################################################################
// Function Decrement_Height
//##################################################################### 
template<class T> void OPENGL_BOX_HIERARCHY_3D<T>::
Decrement_Height()
{
    min_height--;
    max_height--;
}
//#####################################################################
// Function Increment_Height
//##################################################################### 
template<class T> void OPENGL_BOX_HIERARCHY_3D<T>::
Increment_Height()
{
    min_height++;
    max_height++;
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_BOX_HIERARCHY_3D<float>;
template class OPENGL_BOX_HIERARCHY_3D<double>;
}

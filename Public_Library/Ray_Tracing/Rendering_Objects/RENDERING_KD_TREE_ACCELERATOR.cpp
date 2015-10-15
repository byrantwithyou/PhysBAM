//#####################################################################
// Copyright 2004, Ron Fedkiw, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Ray_Tracing/Rendering_Objects/RENDERING_KD_TREE_ACCELERATOR.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> RENDERING_KD_TREE_ACCELERATOR<T>::
RENDERING_KD_TREE_ACCELERATOR()
{
}
//#####################################################################
// Function Add_Object
//#####################################################################
template<class T> void RENDERING_KD_TREE_ACCELERATOR<T>::
Add_Object(RENDERING_OBJECT<T>* object)
{// so this get's into the right bin in render world
    if(object->material_shader)material_shader=object->material_shader; 
    if(object->volumetric_shader)volumetric_shader=object->volumetric_shader;
    objects->Append(object);
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool RENDERING_KD_TREE_ACCELERATOR<T>::
Intersection(RAY<TV>& ray,const int lowest_priority,RENDERING_OBJECT<T>** intersected_object)const
{
    bool hit=false;
    for(int i=0;i<primitives.m;i++){
        if(primitives(i).object->priority>=lowest_priority){
            bool primitive_i_hit=primitives(i).object->Intersection(ray,primitives(i).aggregate_id);
            if(primitive_i_hit){hit=true;*intersected_object=primitives(i).object;}}}
    return hit;
}
//#####################################################################
// Function Preprocess_Efficiency_Structures
//#####################################################################
template<class T> void RENDERING_KD_TREE_ACCELERATOR<T>::
Preprocess_Efficiency_Structures()
{
    for(int i=0;i<objects.m;i++)objects(i)->Get_Aggregate_World_Space_Bounding_Boxes(primitives);
    if(primitives.m>0){
        // construct total bounding box...
        bounding_box.Reset_Bounds(primitives(0).bounding_box);
        for(int i=1;i<objects.m;i++)bounding_box.Enlarge_To_Include_Box(primitives(i).bounding_box);}
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool RENDERING_KD_TREE_ACCELERATOR<T>::
Inside(const TV& location,RENDERING_OBJECT<T>** intersected_object) const
{
    for(int i=0;i<objects.m;i++){
        if(objects(i)->support_transparent_overlapping_objects&&objects(i)->Inside(location)){
            *intersected_object=(RENDERING_OBJECT<T>*)this;return true;}}
    return false;
}
}

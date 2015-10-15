//#####################################################################
// Copyright 2004, Ron Fedkiw, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_KD_TREE_ACCELERATOR
//#####################################################################
#ifndef __RENDERING_KD_TREE_ACCELERATOR__
#define __RENDERING_KD_TREE_ACCELERATOR__

#include <Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T>
class RENDERING_KD_TREE_ACCELERATOR:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::material_shader;using RENDERING_OBJECT<T>::volumetric_shader;

    ARRAY<RENDERING_OBJECT<T>*> objects;
    //KD_TREE_3D<T,RENDERING_OBJECT<T>*> kd_tree;
    ARRAY<RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T> > primitives;
    RANGE<TV> bounding_box;

    RENDERING_KD_TREE_ACCELERATOR();

    void Add_Object(RENDERING_OBJECT<T>* object);
    bool Intersection(RAY<TV>& ray,const int lowest_priority,RENDERING_OBJECT<T>** intersected_object)const override;
    void Preprocess_Efficiency_Structures() override;
    bool Inside(const TV& location,RENDERING_OBJECT<T>** intersected_object) const override;
//#####################################################################
};   
}
#endif


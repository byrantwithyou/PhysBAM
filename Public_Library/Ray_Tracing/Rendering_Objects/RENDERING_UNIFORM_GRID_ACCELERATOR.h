//#####################################################################
// Copyright 2004, Ron Fedkiw, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_UNIFORM_GRID_ACCELERATOR
//#####################################################################
#ifndef __RENDERING_UNIFORM_GRID_ACCELERATOR__
#define __RENDERING_UNIFORM_GRID_ACCELERATOR__

#include <Geometry/Spatial_Acceleration/UNIFORM_BOX_PARTITION.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T>
class RENDERING_UNIFORM_GRID_ACCELERATOR:public RENDERING_OBJECT<T>
{
    using RENDERING_OBJECT<T>::Inside;using RENDERING_OBJECT<T>::Intersection; // silence -Woverloaded-virtual
public:
    using RENDERING_OBJECT<T>::material_shader;using RENDERING_OBJECT<T>::volumetric_shader;

    ARRAY<RENDERING_OBJECT<T>*> objects;
    UNIFORM_BOX_PARTITION<T,RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T>*> uniform_grid;
    ARRAY<RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T> > primitives;
    mutable int operation;
    
    RENDERING_UNIFORM_GRID_ACCELERATOR();
    void Add_Object(RENDERING_OBJECT<T>* object); // so this get's into the right bin in render world
    void Preprocess_Efficiency_Structures(RENDER_WORLD<T>& world) override;
    bool Inside(const VECTOR<T,3>& location,RENDERING_OBJECT<T>** intersected_object) const override;
    TRIANGULATED_SURFACE<T>* Generate_Triangles() const override;
    bool Intersection(RAY<VECTOR<T,3> >& ray,const int lowest_priority,
        const RENDERING_OBJECT<T>** intersected_object) const override;
//#####################################################################
};  
}
#endif


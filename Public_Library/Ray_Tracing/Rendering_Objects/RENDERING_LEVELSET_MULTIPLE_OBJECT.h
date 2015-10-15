//#####################################################################
// Copyright 2005, Jiayi Chong, Jeong-Mo Hong, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_LEVELSET_MULTIPLE_OBJECT
//#####################################################################
#ifndef __RENDERING_LEVELSET_MULTIPLE_OBJECT__
#define __RENDERING_LEVELSET_MULTIPLE_OBJECT__

#include <Incompressible/Level_Sets/LEVELSET_MULTIPLE.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_LEVELSET_MULTIPLE_REGION_OBJECT.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>

namespace PhysBAM{

template<class T_LEVELSET_MULTIPLE>
class RENDERING_LEVELSET_MULTIPLE_OBJECT:public RENDERING_OBJECT<typename T_LEVELSET_MULTIPLE::VECTOR_T::SCALAR>
{
    typedef typename T_LEVELSET_MULTIPLE::VECTOR_T TV;typedef typename TV::SCALAR T;
    using RENDERING_OBJECT<T>::Inside;using RENDERING_OBJECT<T>::Intersection; // silence -Woverloaded-virtual
public:
    using RENDERING_OBJECT<T>::small_number;using RENDERING_OBJECT<T>::inverse_transform;using RENDERING_OBJECT<T>::priority;

    T_LEVELSET_MULTIPLE levelset_multiple;
    ARRAY<RENDERING_LEVELSET_MULTIPLE_REGION_OBJECT<T,T_LEVELSET_MULTIPLE>* > rendering_levelset_multiple_region_objects;
    int number_of_regions;
 
    RENDERING_LEVELSET_MULTIPLE_OBJECT(GRID<TV>& grid_input,ARRAY<ARRAY<T,VECTOR<int,3> > >& phis_input);
    virtual ~RENDERING_LEVELSET_MULTIPLE_OBJECT();
    bool Intersection(RAY<VECTOR<T,3> >& ray,const int lowest_priority,
        const RENDERING_OBJECT<T>** intersected_object) const override;
    int Intersected_Region(RAY<VECTOR<T,3> >& ray) const;
    bool Inside_Region_Only(const VECTOR<T,3>& location,int region_check) const;
    bool Inside(const VECTOR<T,3>& location,RENDERING_OBJECT<T>** intersected_object) const override;
    RANGE<VECTOR<T,3> > Object_Space_Bounding_Box() const  override;
    T Integration_Step(const T phi) const;
    TRIANGULATED_SURFACE<T>* Generate_Triangles() const override;

//#####################################################################
    bool Intersection(RAY<VECTOR<T,3> >& ray,int& region_start,int& region_end,const T thickness=0) const;
//#####################################################################
};   
}
#endif


//#####################################################################
// Copyright 2002-2005, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_PARTICLES
//##################################################################### 
#ifndef __RENDERING_PARTICLES__
#define __RENDERING_PARTICLES__

#include <Tools/Grids_Uniform/GRID.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Intersections/RAY_SPHERE_INTERSECTION.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM {

template<class T>
class RENDERING_PARTICLES:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::small_number;using RENDERING_OBJECT<T>::transform;using RENDERING_OBJECT<T>::Object_Space_Point;
    using RENDERING_OBJECT<T>::Object_Space_Ray;using RENDERING_OBJECT<T>::Inside;using RENDERING_OBJECT<T>::Intersection;
private:
    ARRAY<GEOMETRY_PARTICLES<TV>*,VECTOR<int,3> > particles_array;
    const GRID<TV>& grid;
    T scale;
    ARRAY<ARRAY<int> ,VECTOR<int,3> > particle_to_aggregate_id;
    ARRAY<PAIR<VECTOR<int,3>,int> > aggregate_id_to_particle;
public:

    RENDERING_PARTICLES(const ARRAY<GEOMETRY_PARTICLES<TV>*,VECTOR<int,3> >& p,const GRID<TV>& g,T s);
    virtual ~RENDERING_PARTICLES();

    bool Intersection(RAY<TV>& ray) const override;
    bool Intersection(RAY<TV>& ray,const int aggregate) const  override;
    bool Inside(const TV& location) const  override;
    TV Normal(const TV& location,const int aggregate=1) const  override;
    void Get_Aggregate_World_Space_Bounding_Boxes(ARRAY<RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T> >& primitives)const override;
private:
    void Create_Aggregate_Ids();
//#####################################################################
};
}
#endif

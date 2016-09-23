//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Matrices/FRAME.h>
#include <Geometry/Basic_Geometry/BOUNDED_HORIZONTAL_PLANE.h>
#include <Geometry/Basic_Geometry/BOWL.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Geometry/Basic_Geometry/PLANE.h>
#include <Geometry/Basic_Geometry/RING.h>
#include <Geometry/Basic_Geometry/SMOOTH_GEAR.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Basic_Geometry/TORUS.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_COMBINED.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_COMBINED_EULERIAN.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/MULTIBODY_LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Tessellation/BOUNDED_HORIZONTAL_PLANE_TESSELLATION.h>
#include <Geometry/Tessellation/BOWL_TESSELLATION.h>
#include <Geometry/Tessellation/CYLINDER_TESSELLATION.h>
#include <Geometry/Tessellation/GEAR_TESSELLATION.h>
#include <Geometry/Tessellation/IMPLICIT_OBJECT_TESSELLATION.h>
#include <Geometry/Tessellation/ORIENTED_BOX_TESSELLATION.h>
#include <Geometry/Tessellation/PLANE_TESSELLATION.h>
#include <Geometry/Tessellation/RANGE_TESSELLATION.h>
#include <Geometry/Tessellation/RING_TESSELLATION.h>
#include <Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <Geometry/Tessellation/TORUS_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

namespace PhysBAM{
namespace TESSELLATION{
//#####################################################################
// Function Generate_Triangles_Helper
//#####################################################################
namespace {
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles_Helper(const LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* implicit)
{ // TODO(jontg): There should be a better implementation for this...
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    MARCHING_CUBES<VECTOR<T,3> >::Create_Surface(*surface,implicit->levelset.grid,implicit->levelset.phi);
    return surface;
}

template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles_Helper(const MULTIBODY_LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* implicit)
{
    ARRAY<TRIANGULATED_SURFACE<T>*> triangulated_surfaces(implicit->levelsets->m);
    for(int i=0;i<implicit->levelsets->m;i++) triangulated_surfaces(i)=Generate_Triangles(*(*implicit->levelsets)(i));
    ARRAY<FRAME<VECTOR<T,3> > > identity_frames(implicit->levelsets->m);
    TRIANGULATED_SURFACE<T>* union_object=TRIANGULATED_SURFACE<T>::Union_Mesh_Objects_Relatively(triangulated_surfaces,identity_frames);
    triangulated_surfaces.Delete_Pointers_And_Clean_Memory();
    return union_object;
}

template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles_Helper(const IMPLICIT_OBJECT_COMBINED<VECTOR<T,3> >* implicit)
{
    TRIANGULATED_SURFACE<T>* triangles_1=Generate_Triangles(*implicit->implicit_object1);
    TRIANGULATED_SURFACE<T>* triangles_2=Generate_Triangles(*implicit->implicit_object2);
    assert(triangles_1->particles.Size()==triangles_2->particles.Size());
    for(int i=0;i<triangles_1->particles.Size();i++)
        triangles_1->particles.X(i)=(1-implicit->alpha)*triangles_1->particles.X(i)+implicit->alpha*triangles_2->particles.X(i);
    return triangles_1;
}

template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles_Helper(const IMPLICIT_OBJECT_COMBINED_EULERIAN<VECTOR<T,3> >* implicit)
{PHYSBAM_NOT_IMPLEMENTED();}

template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles_Helper(const IMPLICIT_OBJECT_TRANSFORMED<VECTOR<T,3>,FRAME<VECTOR<T,3> > >* implicit)
{ // TODO(jontg): Better way to templatize this?
    TRIANGULATED_SURFACE<T>* surface=Generate_Triangles(*implicit->object_space_implicit_object);
    for(int i=0;i<surface->particles.Size();i++) surface->particles.X(i)=implicit->World_Space_Point(surface->particles.X(i));
    return surface;
}
}
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> POINT_SIMPLICES_1D<T>* Generate_Triangles(IMPLICIT_OBJECT<VECTOR<T,1> > const& implicit) {PHYSBAM_FATAL_ERROR();}
template<class T> SEGMENTED_CURVE_2D<T>* Generate_Triangles(IMPLICIT_OBJECT<VECTOR<T,2> > const& implicit) {PHYSBAM_FATAL_ERROR();}
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(IMPLICIT_OBJECT<VECTOR<T,3> > const& implicit_input)
{
    typedef VECTOR<T,3> TV;
    if(const LEVELSET_IMPLICIT_OBJECT<TV>* implicit=dynamic_cast<const LEVELSET_IMPLICIT_OBJECT<TV>*>(&implicit_input)) return Generate_Triangles_Helper(implicit);
    else if(const MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>* implicit=dynamic_cast<const MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>*>(&implicit_input)) return Generate_Triangles_Helper(implicit);
    else if(const IMPLICIT_OBJECT_COMBINED<TV>* implicit=dynamic_cast<const IMPLICIT_OBJECT_COMBINED<TV>*>(&implicit_input)) return Generate_Triangles_Helper(implicit); 
    else if(const IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>* implicit=dynamic_cast<const IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>*>(&implicit_input)) return Generate_Triangles_Helper(implicit); 
    else if(const IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* implicit=dynamic_cast<const IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(&implicit_input)) return Generate_Triangles_Helper(implicit); 
    else if(const ANALYTIC_IMPLICIT_OBJECT<PLANE<T> >* implicit=dynamic_cast<const ANALYTIC_IMPLICIT_OBJECT<PLANE<T> >*>(&implicit_input)) return Generate_Triangles(implicit->analytic);
    else if(const ANALYTIC_IMPLICIT_OBJECT<RING<T> >* implicit=dynamic_cast<const ANALYTIC_IMPLICIT_OBJECT<RING<T> >*>(&implicit_input)) return Generate_Triangles(implicit->analytic);
    else if(const ANALYTIC_IMPLICIT_OBJECT<BOWL<T> >* implicit=dynamic_cast<const ANALYTIC_IMPLICIT_OBJECT<BOWL<T> >*>(&implicit_input)) return Generate_Triangles(implicit->analytic);
    else if(const ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >* implicit=dynamic_cast<const ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >*>(&implicit_input)) return Generate_Triangles(implicit->analytic);
    else if(const ANALYTIC_IMPLICIT_OBJECT<SMOOTH_GEAR<TV> >* implicit=dynamic_cast<const ANALYTIC_IMPLICIT_OBJECT<SMOOTH_GEAR<TV> >*>(&implicit_input)) return Generate_Triangles(implicit->analytic);
    else if(const ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<TV> >* implicit=dynamic_cast<const ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<TV> >*>(&implicit_input)) return Generate_Triangles(implicit->analytic);
    else if(const ANALYTIC_IMPLICIT_OBJECT<ORIENTED_BOX<TV> >* implicit=dynamic_cast<const ANALYTIC_IMPLICIT_OBJECT<ORIENTED_BOX<TV> >*>(&implicit_input)) return Generate_Triangles(implicit->analytic);
    else if(const ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >* implicit=dynamic_cast<const ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >*>(&implicit_input)) return Generate_Triangles(implicit->analytic);
    else if(const ANALYTIC_IMPLICIT_OBJECT<TORUS<T> >* implicit=dynamic_cast<const ANALYTIC_IMPLICIT_OBJECT<TORUS<T> >*>(&implicit_input)) return Generate_Triangles(implicit->analytic);
    else if(const ANALYTIC_IMPLICIT_OBJECT<CYLINDER<T> >* implicit=dynamic_cast<const ANALYTIC_IMPLICIT_OBJECT<CYLINDER<T> >*>(&implicit_input)) return Generate_Triangles(implicit->analytic);
    else{
        LOG::cout<<"Trying to generate triangles on an object of type:\n\t"<<typeid(implicit_input).name()<<std::endl;
        PHYSBAM_FATAL_ERROR();}
}
//#####################################################################
#define INSTANTIATION_HELPER(T) \
    template POINT_SIMPLICES_1D<T>* Generate_Triangles(const IMPLICIT_OBJECT<VECTOR<T,1> >&); \
    template SEGMENTED_CURVE_2D<T>* Generate_Triangles(const IMPLICIT_OBJECT<VECTOR<T,2> >&); \
    template TRIANGULATED_SURFACE<T>* Generate_Triangles(const IMPLICIT_OBJECT<VECTOR<T,3> >&); 

INSTANTIATION_HELPER(float);
INSTANTIATION_HELPER(double);
}
}

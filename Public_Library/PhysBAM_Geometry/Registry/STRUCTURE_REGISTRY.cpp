//#####################################################################
// Copyright 2006-2009, Jon Gretarsson, Geoffrey Irving, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOUNDED_HORIZONTAL_PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOWL.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/RING.h>
#include <PhysBAM_Geometry/Basic_Geometry/SMOOTH_GEAR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TORUS.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/MULTIBODY_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TRIANGULATED_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
using namespace PhysBAM;
namespace PhysBAM{
bool Register_Structures(){
    STRUCTURE_REGISTRY<VECTOR<float,1> >::Register<ANALYTIC_IMPLICIT_OBJECT<POINT_SIMPLEX_1D<float> > >();
    STRUCTURE_REGISTRY<VECTOR<float,1> >::Register<ANALYTIC_IMPLICIT_OBJECT<RANGE<VECTOR<float,1> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,1> >::Register<LEVELSET_IMPLICIT_OBJECT<VECTOR<float,1> > >();
    STRUCTURE_REGISTRY<VECTOR<float,1> >::Register<POINT_SIMPLICES_1D<float> >();
    STRUCTURE_REGISTRY<VECTOR<float,1> >::Register<SEGMENTED_CURVE<VECTOR<float,1> > >();
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<VECTOR<float,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<LINE_2D<float> > >();
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<ORIENTED_BOX<VECTOR<float,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<RANGE<VECTOR<float,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<SMOOTH_GEAR<VECTOR<float,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<SPHERE<VECTOR<float,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<EMBEDDED_TRIANGULATED_OBJECT<VECTOR<float,2> > >();
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<FREE_PARTICLES<VECTOR<float,2> > >();
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<IMPLICIT_OBJECT_TRANSFORMED<VECTOR<float,2>,FRAME<VECTOR<float,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<LEVELSET_IMPLICIT_OBJECT<VECTOR<float,2> > >();
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<MULTIBODY_LEVELSET_IMPLICIT_OBJECT<VECTOR<float,2> > >();
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<SEGMENTED_CURVE_2D<float> >();
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<TRIANGULATED_AREA<float> >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<VECTOR<float,3> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<BOWL<float> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<CYLINDER<float> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<ORIENTED_BOX<VECTOR<float,3> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<PLANE<float> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<RANGE<VECTOR<float,3> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<RING<float> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<SMOOTH_GEAR<VECTOR<float,3> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<SPHERE<VECTOR<float,3> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<TORUS<float> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<EMBEDDED_TRIANGULATED_OBJECT<VECTOR<float,3> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<FREE_PARTICLES<VECTOR<float,3> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<IMPLICIT_OBJECT_TRANSFORMED<VECTOR<float,3>,FRAME<VECTOR<float,3> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<LEVELSET_IMPLICIT_OBJECT<VECTOR<float,3> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<MULTIBODY_LEVELSET_IMPLICIT_OBJECT<VECTOR<float,3> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<SEGMENTED_CURVE<VECTOR<float,3> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<TETRAHEDRALIZED_VOLUME<float> >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<TRIANGULATED_SURFACE<float> >();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    STRUCTURE_REGISTRY<VECTOR<double,1> >::Register<ANALYTIC_IMPLICIT_OBJECT<POINT_SIMPLEX_1D<double> > >();
    STRUCTURE_REGISTRY<VECTOR<double,1> >::Register<ANALYTIC_IMPLICIT_OBJECT<RANGE<VECTOR<double,1> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,1> >::Register<LEVELSET_IMPLICIT_OBJECT<VECTOR<double,1> > >();
    STRUCTURE_REGISTRY<VECTOR<double,1> >::Register<POINT_SIMPLICES_1D<double> >();
    STRUCTURE_REGISTRY<VECTOR<double,1> >::Register<SEGMENTED_CURVE<VECTOR<double,1> > >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<VECTOR<double,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<LINE_2D<double> > >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<ORIENTED_BOX<VECTOR<double,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<RANGE<VECTOR<double,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<SMOOTH_GEAR<VECTOR<double,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<SPHERE<VECTOR<double,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<EMBEDDED_TRIANGULATED_OBJECT<VECTOR<double,2> > >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<FREE_PARTICLES<VECTOR<double,2> > >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<IMPLICIT_OBJECT_TRANSFORMED<VECTOR<double,2>,FRAME<VECTOR<double,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<LEVELSET_IMPLICIT_OBJECT<VECTOR<double,2> > >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<MULTIBODY_LEVELSET_IMPLICIT_OBJECT<VECTOR<double,2> > >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<SEGMENTED_CURVE_2D<double> >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<TRIANGULATED_AREA<double> >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<VECTOR<double,3> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<BOWL<double> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<CYLINDER<double> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<ORIENTED_BOX<VECTOR<double,3> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<PLANE<double> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<RANGE<VECTOR<double,3> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<RING<double> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<SMOOTH_GEAR<VECTOR<double,3> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<SPHERE<VECTOR<double,3> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<TORUS<double> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<EMBEDDED_TRIANGULATED_OBJECT<VECTOR<double,3> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<FREE_PARTICLES<VECTOR<double,3> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<IMPLICIT_OBJECT_TRANSFORMED<VECTOR<double,3>,FRAME<VECTOR<double,3> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<LEVELSET_IMPLICIT_OBJECT<VECTOR<double,3> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<MULTIBODY_LEVELSET_IMPLICIT_OBJECT<VECTOR<double,3> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<SEGMENTED_CURVE<VECTOR<double,3> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<TETRAHEDRALIZED_VOLUME<double> >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<TRIANGULATED_SURFACE<double> >();
#endif
    return true;
}
bool registered_structures_asdf=Register_Structures();
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> STRUCTURE_REGISTRY<TV>::
STRUCTURE_REGISTRY()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> STRUCTURE_REGISTRY<TV>::
~STRUCTURE_REGISTRY()
{
    PHYSBAM_ASSERT(registered_structures_asdf);
}
template class STRUCTURE_REGISTRY<VECTOR<float,1> >;
template class STRUCTURE_REGISTRY<VECTOR<float,2> >;
template class STRUCTURE_REGISTRY<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class STRUCTURE_REGISTRY<VECTOR<double,1> >;
template class STRUCTURE_REGISTRY<VECTOR<double,2> >;
template class STRUCTURE_REGISTRY<VECTOR<double,3> >;
#endif

//#####################################################################
// Copyright 2014
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BEZIER_SPLINE
//##################################################################### 
#ifndef __BEZIER_SPLINE__
#define __BEZIER_SPLINE__

#include <Geometry/Topology_Based_Geometry/STRUCTURE.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{
template<class TV> class GEOMETRY_PARTICLES;

template<class TV,int d=3>
class BEZIER_SPLINE:public STRUCTURE<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,d+1> CV;
    bool need_destroy_particles;
public:

    GEOMETRY_PARTICLES<TV>& particles;
    ARRAY<CV> control_points;

    BEZIER_SPLINE();
    BEZIER_SPLINE(GEOMETRY_PARTICLES<TV>& particles);
    virtual ~BEZIER_SPLINE();

//#####################################################################
    static BEZIER_SPLINE* Create();
    static BEZIER_SPLINE* Create(GEOMETRY_PARTICLES<TV>& particles);
    TV Evaluate(int id,T t) const;
    BEZIER_SPLINE<TV,d>* Append_Particles_And_Create_Copy(GEOMETRY_PARTICLES<TV>& new_particles,ARRAY<int>* particle_indices) const;
    void Read(TYPED_ISTREAM& input) PHYSBAM_OVERRIDE;
    void Write(TYPED_OSTREAM& output) const PHYSBAM_OVERRIDE;
//#####################################################################
};
template<class TV> void Smooth_Fit(BEZIER_SPLINE<TV,3>& bs,ARRAY_VIEW<TV> X);
template<class TV> void Smooth_Fit_Loop(BEZIER_SPLINE<TV,3>& bs,ARRAY_VIEW<TV> X);
template<class TV,int d> typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,1>::OBJECT*
Create_Segmented_Curve(const BEZIER_SPLINE<TV,d>& spline,bool same_particles);
}
#endif

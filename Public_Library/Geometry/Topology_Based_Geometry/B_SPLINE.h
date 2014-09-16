//#####################################################################
// Copyright 2014
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class B_SPLINE
//##################################################################### 
#ifndef __B_SPLINE__
#define __B_SPLINE__

#include <Geometry/Topology_Based_Geometry/STRUCTURE.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{
template<class TV> class GEOMETRY_PARTICLES;
template<class TV,int d> class BEZIER_SPLINE;

template<class TV,int d=3>
class B_SPLINE:public STRUCTURE<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,d+1> CV;
    bool need_destroy_particles;
public:

    GEOMETRY_PARTICLES<TV>& particles;
    ARRAY<int> control_points;
    ARRAY<T> knots;

    B_SPLINE();
    B_SPLINE(GEOMETRY_PARTICLES<TV>& particles);
    virtual ~B_SPLINE();

//#####################################################################
    static B_SPLINE* Create();
    static B_SPLINE* Create(GEOMETRY_PARTICLES<TV>& particles);
    TV Evaluate(T t) const;
    B_SPLINE<TV,d>* Append_Particles_And_Create_Copy(GEOMETRY_PARTICLES<TV>& new_particles,ARRAY<int>* particle_indices) const;
    void Read(TYPED_ISTREAM& input) PHYSBAM_OVERRIDE;
    void Write(TYPED_OSTREAM& output) const PHYSBAM_OVERRIDE;
    std::string Name() const PHYSBAM_OVERRIDE;
    static std::string Static_Name();
    std::string Extension() const PHYSBAM_OVERRIDE;
    static std::string Static_Extension();
//#####################################################################
};
template<class TV> void Smooth_Fit(B_SPLINE<TV,3>& bs,ARRAY_VIEW<TV> X);
template<class TV> void Smooth_Fit_Loop(B_SPLINE<TV,3>& bs,ARRAY_VIEW<TV> X);
template<class TV> void Fill_Bezier(BEZIER_SPLINE<TV,3>& bez,const B_SPLINE<TV,3>& bs);
template<class TV,int d> typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,1>::OBJECT*
Create_Segmented_Curve(const B_SPLINE<TV,d>& spline,bool same_particles);
}
#endif

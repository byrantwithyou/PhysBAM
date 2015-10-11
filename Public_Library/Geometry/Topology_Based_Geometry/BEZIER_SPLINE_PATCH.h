//#####################################################################
// Copyright 2014
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BEZIER_SPLINE_PATCH
//##################################################################### 
#ifndef __BEZIER_SPLINE_PATCH__
#define __BEZIER_SPLINE_PATCH__

#include <Geometry/Topology_Based_Geometry/STRUCTURE.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{
template<class TV> class GEOMETRY_PARTICLES;

template<class TV,int d=3>
class BEZIER_SPLINE_PATCH:public STRUCTURE<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,d+1> CV;
    bool need_destroy_particles;
public:

    GEOMETRY_PARTICLES<TV>& particles;
    ARRAY<VECTOR<int,(d+1)*(d+1)> > control_points;

    BEZIER_SPLINE_PATCH();
    BEZIER_SPLINE_PATCH(GEOMETRY_PARTICLES<TV>& particles);
    virtual ~BEZIER_SPLINE_PATCH();

//#####################################################################
    static BEZIER_SPLINE_PATCH* Create();
    static BEZIER_SPLINE_PATCH* Create(GEOMETRY_PARTICLES<TV>& particles);
    TV Evaluate(int id,T s,T t) const;
    BEZIER_SPLINE_PATCH<TV,d>* Append_Particles_And_Create_Copy(GEOMETRY_PARTICLES<TV>& new_particles,ARRAY<int>* particle_indices) const override;
    void Read(TYPED_ISTREAM& input) override;
    void Write(TYPED_OSTREAM& output) const override;
    std::string Name() const override;
    static std::string Static_Name();
    std::string Extension() const override;
    static std::string Static_Extension();
//#####################################################################
};
template<class TV> void Smooth_Fit(BEZIER_SPLINE_PATCH<TV,3>& bs,ARRAY_VIEW<TV> X);
template<class TV> void Smooth_Fit_Loop(BEZIER_SPLINE_PATCH<TV,3>& bs,ARRAY_VIEW<TV> X);
template<class TV,int d> typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,2>::OBJECT*
Create_Triangulated_Object(const BEZIER_SPLINE_PATCH<TV,d>& spline,bool same_particles);
}
#endif

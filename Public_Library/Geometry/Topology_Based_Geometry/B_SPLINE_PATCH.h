//#####################################################################
// Copyright 2014
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class B_SPLINE_PATCH
//##################################################################### 
#ifndef __B_SPLINE_PATCH__
#define __B_SPLINE_PATCH__

#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Matrices/MATRIX_FORWARD.h>
#include <Geometry/Topology_Based_Geometry/STRUCTURE.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{
template<class TV> class GEOMETRY_PARTICLES;
template<class TV,int d> class BEZIER_SPLINE_PATCH;

template<class TV,int d=3>
class B_SPLINE_PATCH:public STRUCTURE<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,2> IV;
    typedef VECTOR<T,(d+1)*(d+1)> WV;
    bool need_destroy_particles;
public:
    GEOMETRY_PARTICLES<TV>& particles;
    ARRAY<int,IV> control_points;
    ARRAY<T> knots_t,knots_s;
    int m;
    bool loop_s,loop_t;

    B_SPLINE_PATCH();
    B_SPLINE_PATCH(GEOMETRY_PARTICLES<TV>& particles);
    virtual ~B_SPLINE_PATCH();

//#####################################################################
    static B_SPLINE_PATCH* Create(); 
    static B_SPLINE_PATCH* Create(GEOMETRY_PARTICLES<TV>& particles);
    TV Evaluate(T s,T t,VECTOR<TV,2>* tangents=0) const;
    void Calculate_Weights(T s,T t,WV& w,WV& dw0,WV& dw1,WV& ddw0,WV& ddw1,WV& ddw2) const;
    VECTOR<int,(d+1)*(d+1)> Control_Points_For_Element(int element) const;
    RANGE<VECTOR<T,2>> Range_For_Element(int element) const;
    B_SPLINE_PATCH<TV,d>* Append_Particles_And_Create_Copy(GEOMETRY_PARTICLES<TV>& new_particles,ARRAY<int>* particle_indices) const;
    void Read(TYPED_ISTREAM& input) PHYSBAM_OVERRIDE;
    void Write(TYPED_OSTREAM& output) const PHYSBAM_OVERRIDE;
    void Wrap(bool new_loop_s,bool new_loop_t);
    std::string Name() const PHYSBAM_OVERRIDE;
    static std::string Static_Name();
    std::string Extension() const PHYSBAM_OVERRIDE;
    static std::string Static_Extension();
//#####################################################################
};
template<class TV> void Smooth_Fit(B_SPLINE_PATCH<TV,3>& bs,const ARRAY<TV,VECTOR<int,2>>& X,bool loop_s=false,bool loop_t=false);
template<class TV> void Fill_Bezier(BEZIER_SPLINE_PATCH<TV,3>& bez,const B_SPLINE_PATCH<TV,3>& bs);
template<class TV,int d> typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,2>::OBJECT*
Create_Triangulated_Surface(const B_SPLINE_PATCH<TV,d>& spline,bool same_particles);
}
#endif

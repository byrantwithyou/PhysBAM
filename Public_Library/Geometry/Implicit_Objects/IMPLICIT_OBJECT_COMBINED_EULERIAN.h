//#####################################################################
// Copyright 2007, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_OBJECT_COMBINED_EULERIAN
//#####################################################################
#ifndef __IMPLICIT_OBJECT_COMBINED_EULERIAN__
#define __IMPLICIT_OBJECT_COMBINED_EULERIAN__

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>

namespace PhysBAM{

template<class TV>
class IMPLICIT_OBJECT_COMBINED_EULERIAN:public IMPLICIT_OBJECT<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::ORIENTED_BOX T_ORIENTED_BOX;
    enum WORKAROUND {d=TV::m};
public:
    typedef IMPLICIT_OBJECT<TV> BASE;
    using BASE::box;

    IMPLICIT_OBJECT<TV> *implicit_object1,*implicit_object2;
    bool owns_implicit_object1,owns_implicit_object2;
    T alpha;T dt; // alpha is the time step fraction, dt is the time step from implicit_object1 to implicit_object2, i.e. the time you want to interpolate at is dt*alpha
    
    IMPLICIT_OBJECT_COMBINED_EULERIAN(IMPLICIT_OBJECT<TV>* implicit_object1_input,bool owns_implicit_object1_input,
        IMPLICIT_OBJECT<TV>* implicit_object2_input, bool owns_implicit_object2_input);

    virtual ~IMPLICIT_OBJECT_COMBINED_EULERIAN();

    void Set_Weights(T alpha_input);
    void Update_Box() override;
    void Update_Minimum_Cell_Size(const int maximum_depth=0) override;
    // TODO: box is not used.  Is it still needed?
    T Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const override; 
    TV Velocity(const TV& location) const override;
    T operator()(const TV& location) const override;
    T Value(const TV& location) const;
    T Extended_Value(const TV& location) const;
    T Signed_Distance(const TV& location) const override; // to make this class compatible with the geometry classes
    T Extended_Phi(const TV& location) const override;
    T Phi_Secondary(const TV& location) const override;
    TV Normal(const TV& location,const int aggregate) const override;
    TV Extended_Normal(const TV& location,const int aggregate=-1) const override;
    void Compute_Normals() override; 
    void Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists=true) override;
    void Rescale(const T scaling_factor) override;
    void Inflate(const T inflation_distance) override;
    bool Lazy_Inside(const TV& location,const T contour_value=0) const override;
    bool Lazy_Inside_And_Value(const TV& location,T& phi,const T contour_value=0) const override;
    bool Lazy_Inside_Extended_Levelset(const TV& location,const T contour_value=0) const override;
    bool Lazy_Inside_Extended_Levelset_And_Value(const TV& location,T& phi_value,const T contour_value=0) const override;
    bool Lazy_Outside(const TV& location,const T contour_value=0) const override;
    bool Lazy_Outside_Extended_Levelset(const TV& location,const T contour_value=0) const override;
    bool Lazy_Outside_Extended_Levelset_And_Value(const TV& location,T& phi_value,const T contour_value=0) const override;
    SYMMETRIC_MATRIX<T,TV::m> Hessian(const TV& X) const override;
    VECTOR<T,d-1> Principal_Curvatures(const TV& X) const override;
    T Integration_Step(const T phi) const override;
    T Minimum_Cell_Size() const override;
    bool Intersection(RAY<TV>& ray,const T thickness) const override;
    virtual void Read(TYPED_ISTREAM& input) override;
    virtual void Write(TYPED_OSTREAM& output) const override;
//#####################################################################
};
}
#endif


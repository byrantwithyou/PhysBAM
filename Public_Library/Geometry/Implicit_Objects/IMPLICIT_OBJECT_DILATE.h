//#####################################################################
// Copyright 2016
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_OBJECT_DILATE
//#####################################################################
#ifndef __IMPLICIT_OBJECT_DILATE__
#define __IMPLICIT_OBJECT_DILATE__

#include <Tools/Math_Tools/RANGE.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Tools/Matrices/MATRIX_POLICY.h>
#include <Tools/Utilities/TYPE_UTILITIES.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Topology_Based_Geometry/STRUCTURE.h>
namespace PhysBAM{

template<class TV>
class IMPLICIT_OBJECT_DILATE:public IMPLICIT_OBJECT<TV>
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
    typedef VECTOR<T,d-1> T_CURVATURES;
public:
    typedef int HAS_TYPED_READ_WRITE;
    typedef TV VECTOR_T;
    
    using IMPLICIT_OBJECT<TV>::box;
    using IMPLICIT_OBJECT<TV>::use_secondary_interpolation;

    IMPLICIT_OBJECT<TV>& io;
    bool owns_io;
    T dilation;

    IMPLICIT_OBJECT_DILATE(IMPLICIT_OBJECT<TV>* o,T dilation)
        :io(*o),owns_io(true),dilation(dilation)
    {
    }
    virtual ~IMPLICIT_OBJECT_DILATE();

    void Use_Secondary_Interpolation(const bool use_secondary_interpolation_input=true)
    {use_secondary_interpolation=use_secondary_interpolation_input;}

//#####################################################################
    void Update_Box() override;
    void Update_Minimum_Cell_Size(const int maximum_depth=0) override;
    T Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const override;
    T operator()(const TV& location) const override;
    T Extended_Phi(const TV& location) const override;
    T Phi_Secondary(const TV& location) const override;
    TV Normal(const TV& location,const int aggregate=-1) const override;
    TV Extended_Normal(const TV& location,const int aggregate=-1) const override;
    void Compute_Normals() override;
    void Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists=true) override;
    void Rescale(const T scaling_factor) override;
    void Translate(const TV& translation) override;
    void Inflate(const T inflation_distance) override;
    bool Lazy_Inside(const TV& location,const T contour_value=0) const override;
    bool Lazy_Inside_And_Value(const TV& location,T& phi_value,const T contour_value=0) const override;
    bool Lazy_Inside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0) const override;
    bool Lazy_Inside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0) const override;
    bool Lazy_Outside(const TV& location,const T contour_value=0) const override;
    bool Lazy_Outside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0) const override;
    bool Lazy_Outside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0) const override;
    TV Velocity(const TV& location) const override;
    // the following only exist in 3d
    SYMMETRIC_MATRIX<T,TV::m> Hessian(const TV& X) const override;
    T_CURVATURES Principal_Curvatures(const TV& X) const override;
    T Integration_Step(const T phi) const override;
    T Minimum_Cell_Size() const override;
    void Read(TYPED_ISTREAM& input) override;
    void Write(TYPED_OSTREAM& output) const override;
//#####################################################################
};
}
#endif

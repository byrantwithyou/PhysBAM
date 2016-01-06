//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Sergey Koltakov, Eran Guendelman, Geoffrey Irving, Neil Molino, Andrew Selle, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_INVERT.h>
#include <Geometry/Intersections/RAY_BOX_INTERSECTION.h>
#include <Geometry/Level_Sets/IMPLICIT_OBJECT_ON_A_RAY.h>
#include <Geometry/Level_Sets/IMPLICIT_OBJECT_ON_A_RAY_SECONDARY_INTERPOLATION.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
namespace PhysBAM{
//#####################################################################
// Destructor
//#####################################################################
template<class TV> IMPLICIT_OBJECT_INVERT<TV>::
~IMPLICIT_OBJECT_INVERT()
{
    if(owns_io) delete &io;
}
//#####################################################################
// Function Update_Box
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_INVERT<TV>::
Update_Box()
{
    box=RANGE<TV>::Full_Box();
}
//#####################################################################
// Function Update_Minimum_Cell_Size
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_INVERT<TV>::
Update_Minimum_Cell_Size(const int maximum_depth)
{
    io.Update_Minimum_Cell_Size(maximum_depth);
}
//#####################################################################
// Function Minimum_Cell_Size_Within_Box
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_INVERT<TV>::
Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const
{
    return io.Minimum_Cell_Size_Within_Box(box);
}
//#####################################################################
// Function operator
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_INVERT<TV>::
operator()(const TV& X) const
{
    return -io(X);
}
//#####################################################################
// Function Extended_Phi
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_INVERT<TV>::
Extended_Phi(const TV& X) const
{
    return -io.Extended_Phi(X);
}
//#####################################################################
// Function Phi_Secondary
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_INVERT<TV>::
Phi_Secondary(const TV& X) const
{
    return -io.Phi_Secondary(X);
}
//#####################################################################
// Function Normal
//#####################################################################
template<class TV> TV IMPLICIT_OBJECT_INVERT<TV>::
Normal(const TV& X,const int aggregate) const
{
    return -io.Normal(X,aggregate);
}
//#####################################################################
// Function Extended_Normal
//#####################################################################
template<class TV> TV IMPLICIT_OBJECT_INVERT<TV>::
Extended_Normal(const TV& X,const int aggregate) const
{
    return -io.Extended_Normal(X,aggregate);
}
//#####################################################################
// Function Compute_Normals
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_INVERT<TV>::
Compute_Normals()
{
    io.Compute_Normals();
}
//#####################################################################
// Function Compute_Cell_Minimum_And_Maximum
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_INVERT<TV>::
Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists)
{
    io.Compute_Cell_Minimum_And_Maximum(recompute_if_exists);
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_INVERT<TV>::
Rescale(const T scaling_factor)
{
    io.Rescale(scaling_factor);
}
//#####################################################################
// Function Translate
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_INVERT<TV>::
Translate(const TV& translation)
{
    io.Translate(translation);
}
//#####################################################################
// Function Inflate
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_INVERT<TV>::
Inflate(const T inflation_distance)
{
    io.Inflate(inflation_distance);
}
//#####################################################################
// Function Inside
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INVERT<TV>::
Inside(const TV& X,const T thickness_over_two) const
{
    return io.Outside(X,thickness_over_two);
}
//#####################################################################
// Function Outside
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INVERT<TV>::
Outside(const TV& X,const T thickness_over_two) const
{
    return io.Inside(X,thickness_over_two);
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INVERT<TV>::
Boundary(const TV& X,const T thickness_over_two) const
{
    return io.Boundary(X,thickness_over_two);
}
//#####################################################################
// Function Lazy_Inside
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INVERT<TV>::
Lazy_Inside(const TV& X,const T contour_value) const
{
    return io.Lazy_Outside(X,-contour_value);
}
//#####################################################################
// Function Lazy_Inside_And_Value
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INVERT<TV>::
Lazy_Inside_And_Value(const TV& X,T& phi_value,const T contour_value) const
{
    bool b=io.Lazy_Inside_And_Value(X,phi_value,-contour_value);
    phi_value=-phi_value;
    return !b;
}
//#####################################################################
// Function Lazy_Inside_Extended_Levelset
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INVERT<TV>::
Lazy_Inside_Extended_Levelset(const TV& X,const T contour_value) const
{
    return io.Lazy_Outside_Extended_Levelset(X,-contour_value);
}
//#####################################################################
// Function Lazy_Inside_Extended_Levelset_And_Value
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INVERT<TV>::
Lazy_Inside_Extended_Levelset_And_Value(const TV& X,T& phi_value,const T contour_value) const
{
    bool b=io.Lazy_Outside_Extended_Levelset_And_Value(X,phi_value,-contour_value);
    phi_value=-phi_value;
    return b;
}
//#####################################################################
// Function Lazy_Outside
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INVERT<TV>::
Lazy_Outside(const TV& X,const T contour_value) const
{
    return io.Lazy_Inside(X,-contour_value);
}
//#####################################################################
// Function Lazy_Outside_Extended_Levelset
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INVERT<TV>::
Lazy_Outside_Extended_Levelset(const TV& X,const T contour_value) const
{
    return io.Lazy_Inside_Extended_Levelset(X,-contour_value);
}
//#####################################################################
// Function Lazy_Outside_Extended_Levelset_And_Value
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INVERT<TV>::
Lazy_Outside_Extended_Levelset_And_Value(const TV& X,T& phi_value,const T contour_value) const
{
    bool b=io.Lazy_Inside_Extended_Levelset_And_Value(X,phi_value,-contour_value);
    phi_value=-phi_value;
    return b;
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INVERT<TV>::
Intersection(RAY<TV>& ray,const T thickness) const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Closest_Point_On_Boundary
//#####################################################################
template<class TV> TV IMPLICIT_OBJECT_INVERT<TV>::
Closest_Point_On_Boundary(const TV& X,const T tolerance,const int max_iterations,T* distance) const
{
    return io.Closest_Point_On_Boundary(X,tolerance,max_iterations,distance);
}
//#####################################################################
// Function Velocity
//#####################################################################
template<class TV> TV IMPLICIT_OBJECT_INVERT<TV>::
Velocity(const TV& X) const
{
    return io.Velocity(X);
}
//#####################################################################
// Function Velocity
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> IMPLICIT_OBJECT_INVERT<TV>::
Hessian(const TV& X) const
{
    return -io.Hessian(X);
}
//#####################################################################
// Function Velocity
//#####################################################################
template<class TV> auto IMPLICIT_OBJECT_INVERT<TV>::
Principal_Curvatures(const TV& X) const -> T_CURVATURES
{
    return -io.Principal_Curvatures(X);
}
//#####################################################################
// Function Integration_Step
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_INVERT<TV>::
Integration_Step(const T phi) const
{
    return io.Integration_Step(phi);
}
//#####################################################################
// Function Minimum_Cell_Size
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_INVERT<TV>::
Minimum_Cell_Size() const
{
    return io.Minimum_Cell_Size();
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_INVERT<TV>::
Read(TYPED_ISTREAM& input)
{
    io.Read(input);
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_INVERT<TV>::
Write(TYPED_OSTREAM& output) const
{
    io.Write(output);
}
template class IMPLICIT_OBJECT_INVERT<VECTOR<float,1> >;
template class IMPLICIT_OBJECT_INVERT<VECTOR<float,2> >;
template class IMPLICIT_OBJECT_INVERT<VECTOR<float,3> >;
template class IMPLICIT_OBJECT_INVERT<VECTOR<double,1> >;
template class IMPLICIT_OBJECT_INVERT<VECTOR<double,2> >;
template class IMPLICIT_OBJECT_INVERT<VECTOR<double,3> >;

}

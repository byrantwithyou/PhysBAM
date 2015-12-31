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
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_INTERSECTION.h>
#include <Geometry/Intersections/RAY_BOX_INTERSECTION.h>
#include <Geometry/Level_Sets/IMPLICIT_OBJECT_ON_A_RAY.h>
#include <Geometry/Level_Sets/IMPLICIT_OBJECT_ON_A_RAY_SECONDARY_INTERPOLATION.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> IMPLICIT_OBJECT_INTERSECTION<TV>::
IMPLICIT_OBJECT_INTERSECTION(IMPLICIT_OBJECT<TV> &a,IMPLICIT_OBJECT<TV> &b)
    :A(a),B(b),owns_A(true),owns_B(true)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> IMPLICIT_OBJECT_INTERSECTION<TV>::
~IMPLICIT_OBJECT_INTERSECTION()
{
    if(owns_A) delete &A;
    if(owns_B) delete &B;
}
//#####################################################################
// Function Update_Box
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_INTERSECTION<TV>::
Update_Box()
{
    A.Update_Box();
    B.Update_Box();
    box=A.box.Intersect(B.box);
}
//#####################################################################
// Function Update_Minimum_Cell_Size
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_INTERSECTION<TV>::
Update_Minimum_Cell_Size(const int maximum_depth)
{
    A.Update_Minimum_Cell_Size(maximum_depth);
    B.Update_Minimum_Cell_Size(maximum_depth);
}
//#####################################################################
// Function Minimum_Cell_Size_Within_Box
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_INTERSECTION<TV>::
Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const
{
    return A.Minimum_Cell_Size_Within_Box(box);
}
//#####################################################################
// Function operator
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_INTERSECTION<TV>::
operator()(const TV& location) const
{
    return max(A(location),B(location));
}
//#####################################################################
// Function Extended_Phi
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_INTERSECTION<TV>::
Extended_Phi(const TV& location) const
{
    return max(A.Extended_Phi(location),B.Extended_Phi(location));
}
//#####################################################################
// Function Phi_Secondary
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_INTERSECTION<TV>::
Phi_Secondary(const TV& location) const
{
    return max(A.Phi_Secondary(location),B.Phi_Secondary(location));
}
//#####################################################################
// Function Normal
//#####################################################################
template<class TV> TV IMPLICIT_OBJECT_INTERSECTION<TV>::
Normal(const TV& location,const int aggregate) const
{
    if(A(location)>B(location)) return A.Normal(location,aggregate);
    return B.Normal(location,aggregate);
}
//#####################################################################
// Function Extended_Normal
//#####################################################################
template<class TV> TV IMPLICIT_OBJECT_INTERSECTION<TV>::
Extended_Normal(const TV& location,const int aggregate) const
{
    if(A.Extended_Phi(location)>B.Extended_Phi(location)) return A.Extended_Normal(location,aggregate);
    return B.Extended_Normal(location,aggregate);
}
//#####################################################################
// Function Compute_Normals
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_INTERSECTION<TV>::
Compute_Normals()
{
    A.Compute_Normals();
    B.Compute_Normals();
}
//#####################################################################
// Function Compute_Cell_Minimum_And_Maximum
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_INTERSECTION<TV>::
Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists)
{
    A.Compute_Cell_Minimum_And_Maximum(recompute_if_exists);
    B.Compute_Cell_Minimum_And_Maximum(recompute_if_exists);
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_INTERSECTION<TV>::
Rescale(const T scaling_factor)
{
    A.Rescale(scaling_factor);
    B.Rescale(scaling_factor);
}
//#####################################################################
// Function Translate
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_INTERSECTION<TV>::
Translate(const TV& translation)
{
    A.Translate(translation);
    B.Translate(translation);
}
//#####################################################################
// Function Inflate
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_INTERSECTION<TV>::
Inflate(const T inflation_distance)
{
    A.Inflate(inflation_distance);
    B.Inflate(inflation_distance);
}
//#####################################################################
// Function Inside
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INTERSECTION<TV>::
Inside(const TV& location,const T thickness_over_two) const
{
    return A.Inside(location,thickness_over_two) && B.Inside(location,thickness_over_two);
}
//#####################################################################
// Function Outside
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INTERSECTION<TV>::
Outside(const TV& location,const T thickness_over_two) const
{
    return A.Outside(location,thickness_over_two) || B.Outside(location,thickness_over_two);
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INTERSECTION<TV>::
Boundary(const TV& location,const T thickness_over_two) const
{
    if(A.Boundary(location,thickness_over_two) && B.Inside(location,-thickness_over_two)) return true;
    if(B.Boundary(location,thickness_over_two) && A.Inside(location,-thickness_over_two)) return true;
    return false;
}
//#####################################################################
// Function Lazy_Inside
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INTERSECTION<TV>::
Lazy_Inside(const TV& location,const T contour_value) const
{
    return A.Lazy_Inside(location,contour_value) && B.Lazy_Inside(location,contour_value);
}
//#####################################################################
// Function Lazy_Inside_And_Value
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INTERSECTION<TV>::
Lazy_Inside_And_Value(const TV& location,T& phi_value,const T contour_value) const
{
    T a,b;
    bool c=A.Lazy_Inside_And_Value(location,a,contour_value);
    bool d=B.Lazy_Inside_And_Value(location,b,contour_value);
    phi_value=max(a,b);
    return c&&d;
}
//#####################################################################
// Function Lazy_Inside_Extended_Levelset
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INTERSECTION<TV>::
Lazy_Inside_Extended_Levelset(const TV& unclamped_X,const T contour_value) const
{
    return A.Lazy_Inside_Extended_Levelset(unclamped_X,contour_value) && B.Lazy_Inside_Extended_Levelset(unclamped_X,contour_value);
}
//#####################################################################
// Function Lazy_Inside_Extended_Levelset_And_Value
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INTERSECTION<TV>::
Lazy_Inside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value) const
{
    T a,b;
    bool c=A.Lazy_Inside_Extended_Levelset_And_Value(unclamped_X,a,contour_value);
    bool d=B.Lazy_Inside_Extended_Levelset_And_Value(unclamped_X,b,contour_value);
    phi_value=max(a,b);
    return c&&d;
}
//#####################################################################
// Function Lazy_Outside
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INTERSECTION<TV>::
Lazy_Outside(const TV& location,const T contour_value) const
{
    return A.Lazy_Outside(location,contour_value) || B.Lazy_Outside(location,contour_value);
}
//#####################################################################
// Function Lazy_Outside_Extended_Levelset
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INTERSECTION<TV>::
Lazy_Outside_Extended_Levelset(const TV& unclamped_X,const T contour_value) const
{
    return A.Lazy_Outside_Extended_Levelset(unclamped_X,contour_value) || B.Lazy_Outside_Extended_Levelset(unclamped_X,contour_value);
}
//#####################################################################
// Function Lazy_Outside_Extended_Levelset_And_Value
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INTERSECTION<TV>::
Lazy_Outside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value) const
{
    T a,b;
    bool c=A.Lazy_Outside_Extended_Levelset_And_Value(unclamped_X,a,contour_value);
    bool d=B.Lazy_Outside_Extended_Levelset_And_Value(unclamped_X,b,contour_value);
    phi_value=max(a,b);
    return c||d;
}
//#####################################################################
// Function Min_Phi
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_INTERSECTION<TV>::
Min_Phi() const
{
    return A.Min_Phi();
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_INTERSECTION<TV>::
Intersection(RAY<TV>& ray,const T thickness) const
{
    RAY<TV> rayB(ray);
    if(!A.Intersection(ray,thickness))
        return false;
    if(!B.Intersection(rayB,thickness)){
        ray=rayB;
        return false;}
    if(ray.t_max<ray.t_max)
        ray=rayB;
    return true;
}
//#####################################################################
// Function Closest_Point_On_Boundary
//#####################################################################
template<class TV> TV IMPLICIT_OBJECT_INTERSECTION<TV>::
Closest_Point_On_Boundary(const TV& location,const T tolerance,const int max_iterations,T* distance) const
{
    if(A(location)>B(location)) return A.Closest_Point_On_Boundary(location,tolerance,max_iterations,distance);
    return B.Closest_Point_On_Boundary(location,tolerance,max_iterations,distance);
}
//#####################################################################
// Function Velocity
//#####################################################################
template<class TV> TV IMPLICIT_OBJECT_INTERSECTION<TV>::
Velocity(const TV& location) const
{
    if(A(location)>B(location)) return A.Velocity(location);
    return B.Velocity(location);
}
//#####################################################################
// Function Velocity
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> IMPLICIT_OBJECT_INTERSECTION<TV>::
Hessian(const TV& X) const
{
    if(A(X)>B(X)) return A.Hessian(X);
    return B.Hessian(X);
}
//#####################################################################
// Function Velocity
//#####################################################################
template<class TV> auto IMPLICIT_OBJECT_INTERSECTION<TV>::
Principal_Curvatures(const TV& X) const -> T_CURVATURES
{
    if(A(X)>B(X)) return A.Principal_Curvatures(X);
    return B.Principal_Curvatures(X);
}
//#####################################################################
// Function Integration_Step
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_INTERSECTION<TV>::
Integration_Step(const T phi) const
{
    return A.Integration_Step(phi);
}
//#####################################################################
// Function Minimum_Cell_Size
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_INTERSECTION<TV>::
Minimum_Cell_Size() const
{
    return A.Minimum_Cell_Size();
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_INTERSECTION<TV>::
Read(TYPED_ISTREAM& input)
{
    A.Read(input);
    B.Read(input);
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_INTERSECTION<TV>::
Write(TYPED_OSTREAM& output) const
{
    A.Write(output);
    B.Write(output);
}
template class IMPLICIT_OBJECT_INTERSECTION<VECTOR<float,1> >;
template class IMPLICIT_OBJECT_INTERSECTION<VECTOR<float,2> >;
template class IMPLICIT_OBJECT_INTERSECTION<VECTOR<float,3> >;
template class IMPLICIT_OBJECT_INTERSECTION<VECTOR<double,1> >;
template class IMPLICIT_OBJECT_INTERSECTION<VECTOR<double,2> >;
template class IMPLICIT_OBJECT_INTERSECTION<VECTOR<double,3> >;

}

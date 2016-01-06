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
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_UNION.h>
#include <Geometry/Intersections/RAY_BOX_INTERSECTION.h>
#include <Geometry/Level_Sets/IMPLICIT_OBJECT_ON_A_RAY.h>
#include <Geometry/Level_Sets/IMPLICIT_OBJECT_ON_A_RAY_SECONDARY_INTERPOLATION.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
namespace PhysBAM{
//#####################################################################
// Destructor
//#####################################################################
template<class TV> IMPLICIT_OBJECT_UNION<TV>::
~IMPLICIT_OBJECT_UNION()
{
    for(int i=0;i<io.m;i++) if(owns_io(i)) delete io(i);
}
//#####################################################################
// Function Update_Box
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_UNION<TV>::
Update_Box()
{
    io(0)->Update_Box();
    box=io(0)->box;
    for(int i=1;i<io.m;i++){
        io(i)->Update_Box();
        box=box.Unite(io(i)->box);}
}
//#####################################################################
// Function Update_Minimum_Cell_Size
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_UNION<TV>::
Update_Minimum_Cell_Size(const int maximum_depth)
{
    for(int i=0;i<io.m;i++) io(i)->Update_Minimum_Cell_Size(maximum_depth);
}
//#####################################################################
// Function Minimum_Cell_Size_Within_Box
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_UNION<TV>::
Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const
{
    return io(0)->Minimum_Cell_Size_Within_Box(box);
}
//#####################################################################
// Function operator
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_UNION<TV>::
operator()(const TV& X) const
{
    T phi=(*io(0))(X);
    for(int i=1;i<io.m;i++)
        phi=min(phi,(*io(i))(X));
    return phi;
}
//#####################################################################
// Function Extended_Phi
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_UNION<TV>::
Extended_Phi(const TV& X) const
{
    T phi=io(0)->Extended_Phi(X);
    for(int i=1;i<io.m;i++)
        phi=min(phi,io(i)->Extended_Phi(X));
    return phi;
}
//#####################################################################
// Function Phi_Secondary
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_UNION<TV>::
Phi_Secondary(const TV& X) const
{
    T phi=io(0)->Phi_Secondary(X);
    for(int i=1;i<io.m;i++)
        phi=min(phi,io(i)->Phi_Secondary(X));
    return phi;
}
//#####################################################################
// Function Active_Levelset
//#####################################################################
template<class TV> int IMPLICIT_OBJECT_UNION<TV>::
Active_Levelset(const TV& X) const
{
    T phi=(*io(0))(X);
    int best_index=0;
    for(int i=1;i<io.m;i++){
        T p=(*io(i))(X);
        if(p<phi){
            phi=p;
            best_index=i;}}
    return best_index;
}
//#####################################################################
// Function Normal
//#####################################################################
template<class TV> TV IMPLICIT_OBJECT_UNION<TV>::
Normal(const TV& X,const int aggregate) const
{
    return io(Active_Levelset(X))->Normal(X,aggregate);
}
//#####################################################################
// Function Extended_Normal
//#####################################################################
template<class TV> TV IMPLICIT_OBJECT_UNION<TV>::
Extended_Normal(const TV& X,const int aggregate) const
{
    T phi=io(0)->Extended_Phi(X);
    int best_index=0;
    for(int i=1;i<io.m;i++){
        T p=io(i)->Extended_Phi(X);
        if(p<phi){
            phi=p;
            best_index=i;}}
    return io(best_index)->Extended_Normal(X,aggregate);
}
//#####################################################################
// Function Compute_Normals
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_UNION<TV>::
Compute_Normals()
{
    for(int i=0;i<io.m;i++) io(i)->Compute_Normals();
}
//#####################################################################
// Function Compute_Cell_Minimum_And_Maximum
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_UNION<TV>::
Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists)
{
    for(int i=0;i<io.m;i++) io(i)->Compute_Cell_Minimum_And_Maximum(recompute_if_exists);
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_UNION<TV>::
Rescale(const T scaling_factor)
{
    for(int i=0;i<io.m;i++) io(i)->Rescale(scaling_factor);
}
//#####################################################################
// Function Translate
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_UNION<TV>::
Translate(const TV& translation)
{
    for(int i=0;i<io.m;i++) io(i)->Translate(translation);
}
//#####################################################################
// Function Inflate
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_UNION<TV>::
Inflate(const T inflation_distance)
{
    for(int i=0;i<io.m;i++) io(i)->Inflate(inflation_distance);
}
//#####################################################################
// Function Inside
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_UNION<TV>::
Inside(const TV& X,const T thickness_over_two) const
{
    for(int i=0;i<io.m;i++)
        if(io(i)->Inside(X,thickness_over_two))
            return true;
    return false;
}
//#####################################################################
// Function Outside
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_UNION<TV>::
Outside(const TV& X,const T thickness_over_two) const
{
    for(int i=0;i<io.m;i++)
        if(!io(i)->Outside(X,thickness_over_two))
            return false;
    return true;
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_UNION<TV>::
Boundary(const TV& X,const T thickness_over_two) const
{
    return !Inside(X,thickness_over_two) && !Outside(X,thickness_over_two);
}
//#####################################################################
// Function Lazy_Inside
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_UNION<TV>::
Lazy_Inside(const TV& X,const T contour_value) const
{
    for(int i=0;i<io.m;i++)
        if(io(i)->Lazy_Inside(X,contour_value))
            return true;
    return false;
}
//#####################################################################
// Function Lazy_Inside_And_Value
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_UNION<TV>::
Lazy_Inside_And_Value(const TV& X,T& phi_value,const T contour_value) const
{
    bool in=io(0)->Lazy_Inside_And_Value(X,phi_value,contour_value);
    for(int i=1;i<io.m;i++){
        T phi=0;
        bool b=io(i)->Lazy_Inside_And_Value(X,phi,contour_value);
        in=in||b;
        phi_value=min(phi_value,phi);}
    return in;
}
//#####################################################################
// Function Lazy_Inside_Extended_Levelset
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_UNION<TV>::
Lazy_Inside_Extended_Levelset(const TV& X,const T contour_value) const
{
    for(int i=0;i<io.m;i++)
        if(io(i)->Lazy_Inside_Extended_Levelset(X,contour_value))
            return true;
    return false;
}
//#####################################################################
// Function Lazy_Inside_Extended_Levelset_And_Value
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_UNION<TV>::
Lazy_Inside_Extended_Levelset_And_Value(const TV& X,T& phi_value,const T contour_value) const
{
    bool in=io(0)->Lazy_Inside_Extended_Levelset_And_Value(X,phi_value,contour_value);
    for(int i=1;i<io.m;i++){
        T phi=0;
        bool b=io(i)->Lazy_Inside_Extended_Levelset_And_Value(X,phi,contour_value);
        in=in||b;
        phi_value=min(phi_value,phi);}
    return in;
}
//#####################################################################
// Function Lazy_Outside
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_UNION<TV>::
Lazy_Outside(const TV& X,const T contour_value) const
{
    for(int i=0;i<io.m;i++)
        if(!io(i)->Lazy_Outside(X,contour_value))
            return false;
    return true;
}
//#####################################################################
// Function Lazy_Outside_Extended_Levelset
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_UNION<TV>::
Lazy_Outside_Extended_Levelset(const TV& X,const T contour_value) const
{
    for(int i=0;i<io.m;i++)
        if(!io(i)->Lazy_Outside_Extended_Levelset(X,contour_value))
            return false;
    return true;
}
//#####################################################################
// Function Lazy_Outside_Extended_Levelset_And_Value
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_UNION<TV>::
Lazy_Outside_Extended_Levelset_And_Value(const TV& X,T& phi_value,const T contour_value) const
{
    bool in=io(0)->Lazy_Inside_Extended_Levelset_And_Value(X,phi_value,contour_value);
    for(int i=1;i<io.m;i++){
        T phi=0;
        bool b=io(i)->Lazy_Inside_Extended_Levelset_And_Value(X,phi,contour_value);
        in=in&&b;
        phi_value=min(phi_value,phi);}
    return in;
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_UNION<TV>::
Intersection(RAY<TV>& ray,const T thickness) const
{
    bool hit=false;
    for(int i=0;i<io.m;i++)
        if(io(i)->Intersection(ray,thickness))
            hit=true;
    return hit;
}
//#####################################################################
// Function Closest_Point_On_Boundary
//#####################################################################
template<class TV> TV IMPLICIT_OBJECT_UNION<TV>::
Closest_Point_On_Boundary(const TV& X,const T tolerance,const int max_iterations,T* distance) const
{
    return io(Active_Levelset(X))->Closest_Point_On_Boundary(X,tolerance,max_iterations,distance);
}
//#####################################################################
// Function Velocity
//#####################################################################
template<class TV> TV IMPLICIT_OBJECT_UNION<TV>::
Velocity(const TV& X) const
{
    return io(Active_Levelset(X))->Velocity(X);
}
//#####################################################################
// Function Velocity
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> IMPLICIT_OBJECT_UNION<TV>::
Hessian(const TV& X) const
{
    return io(Active_Levelset(X))->Hessian(X);
}
//#####################################################################
// Function Velocity
//#####################################################################
template<class TV> auto IMPLICIT_OBJECT_UNION<TV>::
Principal_Curvatures(const TV& X) const -> T_CURVATURES
{
    return io(Active_Levelset(X))->Principal_Curvatures(X);
}
//#####################################################################
// Function Integration_Step
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_UNION<TV>::
Integration_Step(const T phi) const
{
    return io(0)->Integration_Step(phi);
}
//#####################################################################
// Function Minimum_Cell_Size
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_UNION<TV>::
Minimum_Cell_Size() const
{
    return io(0)->Minimum_Cell_Size();
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_UNION<TV>::
Read(TYPED_ISTREAM& input)
{
    for(int i=0;i<io.m;i++) io(i)->Read(input);
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_UNION<TV>::
Write(TYPED_OSTREAM& output) const
{
    for(int i=0;i<io.m;i++) io(i)->Write(output);
}
template class IMPLICIT_OBJECT_UNION<VECTOR<float,1> >;
template class IMPLICIT_OBJECT_UNION<VECTOR<float,2> >;
template class IMPLICIT_OBJECT_UNION<VECTOR<float,3> >;
template class IMPLICIT_OBJECT_UNION<VECTOR<double,1> >;
template class IMPLICIT_OBJECT_UNION<VECTOR<double,2> >;
template class IMPLICIT_OBJECT_UNION<VECTOR<double,3> >;

}

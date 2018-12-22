//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Sergey Koltakov, Eran Guendelman, Geoffrey Irving, Neil Molino, Andrew Selle, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Intersections/RAY_BOX_INTERSECTION.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> IMPLICIT_OBJECT<TV>::
IMPLICIT_OBJECT()
{
    Use_Secondary_Interpolation(false);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> IMPLICIT_OBJECT<TV>::
~IMPLICIT_OBJECT()
{}
//#####################################################################
// Function Iterative_Solver_Tolerance
//#####################################################################
namespace{
template<class T> inline T Iterative_Solver_Tolerance(){STATIC_ASSERT((T)false);}
template<> inline float Iterative_Solver_Tolerance(){return (float).01;}
template<> inline double Iterative_Solver_Tolerance(){return .001;}
}
//#####################################################################
// Function Box
//#####################################################################
template<class TV> RANGE<TV>& IMPLICIT_OBJECT<TV>::
Box()
{
    return box;
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT<TV>::
Intersection(RAY<TV>& ray,const T thickness) const
{
    bool save_semi_infinite=ray.semi_infinite;T save_t_max=ray.t_max;int save_aggregate=ray.aggregate_id;

    T t_start,t_end;
    bool exit_intersection=false;int exit_aggregate=0;T exit_t_max=0;
    int intersect_box=INTERSECTION::Intersects(ray,box,thickness);
    int outside_box=box.Outside(ray.endpoint,thickness);

    if(outside_box && !intersect_box) return false; // missed the box
    else if(outside_box){ // intersected the box from the outside
        TV point=ray.Point(ray.t_max); // intersection point with the box
        point=box.Thickened(-4*thickness).Clamp(point); // moves the point inside the box
        if((*this)(point) <= 0) return true; // level set is on the edge of the box
        else{ // level set is not on the edge of the box
            t_start=ray.t_max; // intersection with the box
            ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate; // box intersection doesn't count
            RAY<TV> new_ray(point,ray.direction);INTERSECTION::Intersects(new_ray,box,thickness);
            if(ray.semi_infinite) t_end=t_start+new_ray.t_max;
            else t_end=min(t_start+new_ray.t_max,ray.t_max);}}
    else if(!intersect_box){t_start=0;t_end=ray.t_max;} // intersected some object inside the box
    else{ // intersects the box from inside
        t_start=0;t_end=ray.t_max;
        exit_intersection=true;exit_t_max=ray.t_max;exit_aggregate=ray.aggregate_id; // save for exiting rays
        ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate;} // box intersection doesn't count

    T tolerance=Iterative_Solver_Tolerance<T>()*thickness;
    auto implicit_surface_on_a_ray=[this,&ray](T x){return (*this)(ray.Point(x));};
    if(!use_secondary_interpolation){
        // set up marching
        // start marching
        T t1=t_start+thickness,phi1=(*this)(ray.Point(t1)),t2=t1+Integration_Step(phi1);
        // march through the line segment
        while(t2 <= t_end){
            T phi2=(*this)(ray.Point(t2));
            if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){
                ray.semi_infinite=false;
                ray.t_max=Bisection_Secant_Root(implicit_surface_on_a_ray,t1,t2,tolerance);
                ray.aggregate_id=-1;
                return true;}
            else{t1=t2;phi1=phi2;t2=t1+Integration_Step(phi1);}}
        // check the last piece of the line segment
        t2=t_end;T phi2=(*this)(ray.Point(t_end));
        if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){
            ray.semi_infinite=false;
            ray.t_max=Bisection_Secant_Root(implicit_surface_on_a_ray,t1,t2,tolerance);
            ray.aggregate_id=-1;
            return true;}
        else if(exit_intersection && phi2 <= 0){ray.semi_infinite=false;ray.t_max=exit_t_max;ray.aggregate_id=exit_aggregate;return true;} // exiting ray
        else return false;}
    else{ // use_secondary_interpolation
        // set up marching
        // start marching
        T t1=t_start+thickness,phi1=(*this)(ray.Point(t1)),t2=t1+Integration_Step(phi1);
        // march through the line segment

        while(t2 <= t_end){
            T phi2=(*this)(ray.Point(t2));
            if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){
                phi1=this->Phi_Secondary(ray.Point(t1));
                phi2=this->Phi_Secondary(ray.Point(t2));
                if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){
                    ray.semi_infinite=false;
                    ray.t_max=Bisection_Secant_Root(implicit_surface_on_a_ray,t1,t2,tolerance);
                    ray.aggregate_id=-1;return true;}
                else{t1=t2;phi1=phi2;t2=t1+Integration_Step(phi1);}}
            else{t1=t2;phi1=phi2;t2=t1+Integration_Step(phi1);}}

        // check the last piece of the line segment
        t2=t_end;T phi2=(*this)(ray.Point(t2));
        if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){
            phi1=this->Phi_Secondary(ray.Point(t1));
            phi2=this->Phi_Secondary(ray.Point(t2));
            if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){
                ray.semi_infinite=false;
                ray.t_max=Bisection_Secant_Root(implicit_surface_on_a_ray,t1,t2,tolerance);
                ray.aggregate_id=-1;return true;}}
        if(exit_intersection && phi2 <= 0){ray.semi_infinite=false;ray.t_max=exit_t_max;ray.aggregate_id=exit_aggregate;return true;} // exiting ray
        else return false;}
}
//#####################################################################
// Function Closest_Point_On_Boundary
//#####################################################################
template<class TV> TV IMPLICIT_OBJECT<TV>::
Closest_Point_On_Boundary(const TV& location,const T tolerance,const int max_iterations,T* distance) const
{
    T dist=(*this)(location);
    TV new_location=location-dist*Normal(location);
    if(!tolerance){
        if(distance) *distance=dist;
        return new_location;} // only take one iteration
    for(int i=1;i<max_iterations;i++){
        T dist=(*this)(new_location);
        if(abs(dist)<=tolerance) break;
        new_location-=dist*Normal(new_location);}
    if(distance) *distance=(new_location-location).Magnitude();
    return new_location;
}
//#####################################################################
// Function Signed_Distance
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT<TV>::
Signed_Distance(const TV& location) const
{
    return (*this)(location);
}
//#####################################################################
// Function Inside
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT<TV>::
Inside(const TV& location,const T thickness_over_two) const
{
    return box.Inside(location,thickness_over_two) && (*this)(location)<=-thickness_over_two;
}
//#####################################################################
// Function Outside
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT<TV>::
Outside(const TV& location,const T thickness_over_two) const
{
    return box.Outside(location,thickness_over_two) || (*this)(location)>=thickness_over_two;
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT<TV>::
Boundary(const TV& location,const T thickness_over_two) const
{
    return !Inside(location,thickness_over_two) && !Outside(location,thickness_over_two);
}
//#####################################################################
// Function Test_Diff
//#####################################################################
template<class TV> void IMPLICIT_OBJECT<TV>::
Test_Diff(const RANGE<TV>& range,bool test_hess) const
{
    T e=(T)1e-6;
    RANDOM_NUMBERS<T> rand;
    TV X,dX;
    rand.Fill_Uniform(X,range);
    rand.Fill_Uniform(dX,-e,e);

    T p0=Signed_Distance(X),p1=Signed_Distance(X+dX);
    TV n0=Normal(X),n1=Normal(X+dX);
    T diff_pa=(p1-p0)/e,diff_pb=dX.Dot(n0+n1)/(2*e);
    LOG::printf("IMPLICIT_OBJECT test pt %P\n", X);
    LOG::printf("IMPLICIT_OBJECT diff test %g   (%g %g)\n", abs(diff_pa-diff_pb)/maxabs(diff_pa,diff_pb,(T)1e-30),diff_pa,diff_pb);
    if(test_hess){
        SYMMETRIC_MATRIX<T,TV::m> H0=Hessian(X),H1=Hessian(X+dX);
        TV hess_pa=(n1-n0)/e,hess_pb=(H0+H1)*dX/(2*e);
        T mag_pa=hess_pa.Magnitude(),mag_pb=hess_pb.Magnitude();
        LOG::printf("IMPLICIT_OBJECT hess test %g   (%g %g)\n", (hess_pa-hess_pb).Magnitude()/maxabs(mag_pa,mag_pb,(T)1e-30),mag_pa,mag_pb);}
}
//#####################################################################
// Function Principal_Curvatures_Helper
//#####################################################################
template<class T> VECTOR<T,0> inline
Principal_Curvatures_Helper(const VECTOR<T,1>& n,const SYMMETRIC_MATRIX<T,1>& H)
{
    return VECTOR<T,0>();
}
//#####################################################################
// Function Principal_Curvatures_Helper
//#####################################################################
template<class T> VECTOR<T,1> inline
Principal_Curvatures_Helper(const VECTOR<T,2>& n,const SYMMETRIC_MATRIX<T,2>& H)
{
    VECTOR<T,2> t=n.Perpendicular();
    return VECTOR<T,1>(t.Dot(H*t));
}
//#####################################################################
// Function Principal_Curvatures_Helper
//#####################################################################
template<class T> VECTOR<T,2> inline
Principal_Curvatures_Helper(const VECTOR<T,3>& n,const SYMMETRIC_MATRIX<T,3>& H)
{
    SYMMETRIC_MATRIX<T,3> P=(T)1-SYMMETRIC_MATRIX<T,3>::Outer_Product(n),M=SYMMETRIC_MATRIX<T,3>::Conjugate(P,H);
    T trace=M.Trace();
    QUADRATIC<T> quadratic(-1,trace,sqr(M(0,2))+sqr(M(0,1))+sqr(M(1,2))-M(1,1)*M(0,0)-M(2,2)*M(0,0)-M(2,2)*M(1,1));
    quadratic.Compute_Roots();
    if(quadratic.roots == 0) return (T).5*VECTOR<T,2>(trace,trace);
    if(quadratic.roots == 1) return VECTOR<T,2>(quadratic.root[0],quadratic.root[0]);
    return VECTOR<T,2>(quadratic.root[0],quadratic.root[1]);
}
//#####################################################################
// Function Compute_Principal_Curvatures
//#####################################################################
template<class T,int d> VECTOR<T,d-1>
Compute_Principal_Curvatures(const VECTOR<T,d>& n,const SYMMETRIC_MATRIX<T,d>& H)
{
    return Principal_Curvatures_Helper(n,H);
}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class TV> auto IMPLICIT_OBJECT<TV>::
Principal_Curvatures(const TV& X) const -> VECTOR<T,TV::m-1>
{
    return Compute_Principal_Curvatures(Normal(X),Hessian(X));
}
//#####################################################################
// Undefined Functions
//#####################################################################
template<class TV> void IMPLICIT_OBJECT<TV>::Update_Box(){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void IMPLICIT_OBJECT<TV>::Update_Minimum_Cell_Size(const int maximum_depth){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT<TV>::Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT<TV>::operator()(const TV& location) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT<TV>::Extended_Phi(const TV& location) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT<TV>::Phi_Secondary(const TV& location) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> TV IMPLICIT_OBJECT<TV>::Normal(const TV& location,const int aggregate) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> TV IMPLICIT_OBJECT<TV>::Extended_Normal(const TV& location,const int aggregate) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void IMPLICIT_OBJECT<TV>::Compute_Normals(){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void IMPLICIT_OBJECT<TV>::Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void IMPLICIT_OBJECT<TV>::Rescale(const T scaling_factor) {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void IMPLICIT_OBJECT<TV>::Translate(const TV& translation){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void IMPLICIT_OBJECT<TV>::Inflate(const T inflation_distance){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool IMPLICIT_OBJECT<TV>::Lazy_Inside(const TV& location,const T contour_value) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool IMPLICIT_OBJECT<TV>::Lazy_Inside_And_Value(const TV& location,T& phi_value,const T contour_value) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool IMPLICIT_OBJECT<TV>::Lazy_Inside_Extended_Levelset(const TV& unclamped_X,const T contour_value) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool IMPLICIT_OBJECT<TV>::Lazy_Inside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool IMPLICIT_OBJECT<TV>::Lazy_Outside(const TV& location,const T contour_value) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool IMPLICIT_OBJECT<TV>::Lazy_Outside_Extended_Levelset(const TV& unclamped_X,const T contour_value) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool IMPLICIT_OBJECT<TV>::Lazy_Outside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> TV IMPLICIT_OBJECT<TV>::Velocity(const TV& location) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> IMPLICIT_OBJECT<TV>::Hessian(const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT<TV>::Integration_Step(const T phi) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT<TV>::Minimum_Cell_Size() const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();} 
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template IMPLICIT_OBJECT<VECTOR<T,d> >::IMPLICIT_OBJECT(); \
    template IMPLICIT_OBJECT<VECTOR<T,d> >::~IMPLICIT_OBJECT(); \
    template RANGE<VECTOR<T,d> >& IMPLICIT_OBJECT<VECTOR<T,d> >::Box(); \
    template void IMPLICIT_OBJECT<VECTOR<T,d> >::Inflate(const T inflation_distance); \
    template void IMPLICIT_OBJECT<VECTOR<T,d> >::Rescale(const T scaling_factor); \
    template void IMPLICIT_OBJECT<VECTOR<T,d> >::Translate(const VECTOR<T,d> & translation); \
    template T IMPLICIT_OBJECT<VECTOR<T,d> >::Phi_Secondary(const VECTOR<T,d> & location) const; \
    template T IMPLICIT_OBJECT<VECTOR<T,d> >::Signed_Distance(const VECTOR<T,d> & location) const; \
    template T IMPLICIT_OBJECT<VECTOR<T,d> >::Integration_Step(const T phi) const; \
    template VECTOR<T,d>  IMPLICIT_OBJECT<VECTOR<T,d> >::Closest_Point_On_Boundary(const VECTOR<T,d> & location,const T tolerance,const int max_iterations,T* distance) const; \
    template bool IMPLICIT_OBJECT<VECTOR<T,d> >::Inside(const VECTOR<T,d> & location,const T thickness_over_two) const; \
    template bool IMPLICIT_OBJECT<VECTOR<T,d> >::Outside(const VECTOR<T,d> & location,const T thickness_over_two) const; \
    template bool IMPLICIT_OBJECT<VECTOR<T,d> >::Boundary(const VECTOR<T,d> & location,const T thickness_over_two) const; \
    template VECTOR<T,d>  IMPLICIT_OBJECT<VECTOR<T,d> >::Velocity(const VECTOR<T,d> & location) const; \
    template void IMPLICIT_OBJECT<VECTOR<T,d> >::Compute_Normals(); \
    template void IMPLICIT_OBJECT<VECTOR<T,d> >::Update_Minimum_Cell_Size(const int maximum_depth); \
    template void IMPLICIT_OBJECT<VECTOR<T,d> >::Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists); \
    template bool IMPLICIT_OBJECT<VECTOR<T,d> >::Lazy_Inside(const VECTOR<T,d>& X,const T contour_value) const; \
    template bool IMPLICIT_OBJECT<VECTOR<T,d> >::Lazy_Inside_And_Value(const VECTOR<T,d>& X,T& phi_value,const T contour_value) const; \
    template bool IMPLICIT_OBJECT<VECTOR<T,d> >::Lazy_Inside_Extended_Levelset(const VECTOR<T,d>& unclamped_X,const T contour_value) const; \
    template bool IMPLICIT_OBJECT<VECTOR<T,d> >::Lazy_Inside_Extended_Levelset_And_Value(const VECTOR<T,d>& unclamped_X,T& phi_value,const T contour_value) const; \
    template bool IMPLICIT_OBJECT<VECTOR<T,d> >::Lazy_Outside(const VECTOR<T,d>& X,const T contour_value) const; \
    template bool IMPLICIT_OBJECT<VECTOR<T,d> >::Lazy_Outside_Extended_Levelset(const VECTOR<T,d>& unclamped_X,const T contour_value) const; \
    template bool IMPLICIT_OBJECT<VECTOR<T,d> >::Lazy_Outside_Extended_Levelset_And_Value(const VECTOR<T,d>& unclamped_X,T& phi_value,const T contour_value) const; \
    template SYMMETRIC_MATRIX<T,d> IMPLICIT_OBJECT<VECTOR<T,d> >::Hessian(const VECTOR<T,d>& X) const;
INSTANTIATION_HELPER(float,1)
INSTANTIATION_HELPER(float,2)
INSTANTIATION_HELPER(float,3)
template bool IMPLICIT_OBJECT<VECTOR<float,1> >::Intersection(RAY<VECTOR<float,1> >&,const float) const;
template bool IMPLICIT_OBJECT<VECTOR<float,2> >::Intersection(RAY<VECTOR<float,2> >&,const float) const;
template bool IMPLICIT_OBJECT<VECTOR<float,3> >::Intersection(RAY<VECTOR<float,3> >&,const float) const;
template void IMPLICIT_OBJECT<VECTOR<float,1> >::Test_Diff(RANGE<VECTOR<float,1> > const&,bool) const;
template void IMPLICIT_OBJECT<VECTOR<float,2> >::Test_Diff(RANGE<VECTOR<float,2> > const&,bool) const;
template void IMPLICIT_OBJECT<VECTOR<float,3> >::Test_Diff(RANGE<VECTOR<float,3> > const&,bool) const;
INSTANTIATION_HELPER(double,1)
INSTANTIATION_HELPER(double,2)
INSTANTIATION_HELPER(double,3)
template bool IMPLICIT_OBJECT<VECTOR<double,1> >::Intersection(RAY<VECTOR<double,1> >&,const double) const;
template bool IMPLICIT_OBJECT<VECTOR<double,2> >::Intersection(RAY<VECTOR<double,2> >&,const double) const;
template bool IMPLICIT_OBJECT<VECTOR<double,3> >::Intersection(RAY<VECTOR<double,3> >&,const double) const;
template void IMPLICIT_OBJECT<VECTOR<double,1> >::Test_Diff(RANGE<VECTOR<double,1> > const&,bool) const;
template void IMPLICIT_OBJECT<VECTOR<double,2> >::Test_Diff(RANGE<VECTOR<double,2> > const&,bool) const;
template void IMPLICIT_OBJECT<VECTOR<double,3> >::Test_Diff(RANGE<VECTOR<double,3> > const&,bool) const;
template VECTOR<double,0> Compute_Principal_Curvatures<double,1>(VECTOR<double,1> const&,SYMMETRIC_MATRIX<double,1> const&);
template VECTOR<double,1> Compute_Principal_Curvatures<double,2>(VECTOR<double,2> const&,SYMMETRIC_MATRIX<double,2> const&);
template VECTOR<double,2> Compute_Principal_Curvatures<double,3>(VECTOR<double,3> const&,SYMMETRIC_MATRIX<double,3> const&);
template VECTOR<float,0> Compute_Principal_Curvatures<float,1>(VECTOR<float,1> const&,SYMMETRIC_MATRIX<float,1> const&);
template VECTOR<float,1> Compute_Principal_Curvatures<float,2>(VECTOR<float,2> const&,SYMMETRIC_MATRIX<float,2> const&);
template VECTOR<float,2> Compute_Principal_Curvatures<float,3>(VECTOR<float,3> const&,SYMMETRIC_MATRIX<float,3> const&);
}

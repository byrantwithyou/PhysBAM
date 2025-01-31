//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_OBJECT_TRANSFORMED
//#####################################################################
#include <Core/Matrices/FRAME.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Core/Vectors/VECTOR_3D.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Registry/STRUCTURE_REGISTRY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class TRANSFORM> IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<TV>* object_space_implicit_object_input,bool owns_implicit_object_input,const TRANSFORM* transform_input)
    :IMPLICIT_OBJECT_TRANSFORMED_HELPER<TV,TRANSFORM>(transform_input),owns_implicit_object(owns_implicit_object_input),object_space_implicit_object(object_space_implicit_object_input)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class TRANSFORM> IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<TV>* object_space_implicit_object_input,bool owns_implicit_object_input,const TV& center_input,T scale_input)
    :IMPLICIT_OBJECT_TRANSFORMED_HELPER<TV,TRANSFORM>(center_input,scale_input),owns_implicit_object(owns_implicit_object_input),object_space_implicit_object(object_space_implicit_object_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,class TRANSFORM> IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
~IMPLICIT_OBJECT_TRANSFORMED()
{
    if(owns_implicit_object) delete object_space_implicit_object;
}
//#####################################################################
// Function Box
//#####################################################################
template<class TV,class TRANSFORM> RANGE<TV>& IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Box()
{
    Update_Box();
    return box;
}
//#####################################################################
// Function Update_Box
//#####################################################################
template<class TV,class TRANSFORM> void IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Update_Box()
{
    object_space_implicit_object->Update_Box();
    box=World_Space_Box(object_space_implicit_object->box);
}
//#####################################################################
// Function Update_Minimum_Cell_Size
//#####################################################################
template<class TV,class TRANSFORM> void IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Update_Minimum_Cell_Size(const int maximum_depth)
{
    object_space_implicit_object->Update_Minimum_Cell_Size(maximum_depth);
}
//#####################################################################
// Function Minimum_Cell_Size_Within_Box
//#####################################################################
template<class TV,class TRANSFORM> typename TV::SCALAR IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const
{
    return World_Space_Length(object_space_implicit_object->Minimum_Cell_Size_Within_Box(Object_Space_Box(object_space_implicit_object->box)));
}
//#####################################################################
// Function operator()
//#####################################################################
template<class TV,class TRANSFORM> typename TV::SCALAR IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
operator()(const TV& location) const
{
    return World_Space_Length((*object_space_implicit_object)(Object_Space_Point(location)));
}
//#####################################################################
// Function Signed_Distance
//#####################################################################
template<class TV,class TRANSFORM> typename TV::SCALAR IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Signed_Distance(const TV& location) const
{
    return (*this)(location); // to make this class compatible with the geometry classes
}
//#####################################################################
// Function Extended_Phi
//#####################################################################
template<class TV,class TRANSFORM> typename TV::SCALAR IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Extended_Phi(const TV& location) const
{
    return Extended_Value(location);
}
//#####################################################################
// Function Phi_Secondary
//#####################################################################
template<class TV,class TRANSFORM> typename TV::SCALAR IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Phi_Secondary(const TV& location) const
{
    return World_Space_Length(object_space_implicit_object->Phi_Secondary(Object_Space_Point(location)));
}
//#####################################################################
// Function Normal
//#####################################################################
template<class TV,class TRANSFORM> TV IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Normal(const TV& location,const int aggregate) const
{
    return World_Space_Unitless_Vector(object_space_implicit_object->Normal(Object_Space_Point(location),aggregate));
}
//#####################################################################
// Function Extended_Normal
//#####################################################################
template<class TV,class TRANSFORM> TV IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Extended_Normal(const TV& location,const int aggregate) const
{
    return World_Space_Unitless_Vector(object_space_implicit_object->Extended_Normal(Object_Space_Point(location),aggregate));
}
//#####################################################################
// Function Compute_Normals
//#####################################################################
template<class TV,class TRANSFORM> void IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Compute_Normals()
{
    object_space_implicit_object->Compute_Normals();
}
//#####################################################################
// Function Compute_Cell_Minimum_And_Maximum
//#####################################################################
template<class TV,class TRANSFORM> void IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists)
{
    object_space_implicit_object->Compute_Cell_Minimum_And_Maximum(recompute_if_exists);
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class TV,class TRANSFORM> void IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Rescale(const T scaling_factor)
{
    if(!Scale_Transform(scaling_factor)) object_space_implicit_object->Rescale(scaling_factor);
}
//#####################################################################
// Function Inflate
//#####################################################################
template<class TV,class TRANSFORM> void IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Inflate(const T inflation_distance)
{
    object_space_implicit_object->Inflate(Object_Space_Length(inflation_distance));
}
//#####################################################################
// Function Inside
//#####################################################################
template<class TV,class TRANSFORM> bool IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Inside(const TV& location,const T thickness_over_two) const
{
    return object_space_implicit_object->Inside(Object_Space_Point(location),Object_Space_Length(thickness_over_two));
}
//#####################################################################
// Function Lazy_Inside
//#####################################################################
template<class TV,class TRANSFORM> bool IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Lazy_Inside(const TV& location,const T contour_value) const
{
    return object_space_implicit_object->Lazy_Inside(Object_Space_Point(location),Object_Space_Length(contour_value));
}
//#####################################################################
// Function Lazy_Inside_And_Value
//#####################################################################
template<class TV,class TRANSFORM> bool IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Lazy_Inside_And_Value(const TV& location,T& phi,const T contour_value) const
{
    bool result=object_space_implicit_object->Lazy_Inside_And_Value(Object_Space_Point(location),phi,Object_Space_Length(contour_value));
    phi=World_Space_Length(phi);
    return result;
}
//#####################################################################
// Function Lazy_Inside_Extended_Levelset
//#####################################################################
template<class TV,class TRANSFORM> bool IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Lazy_Inside_Extended_Levelset(const TV& location,const T contour_value) const
{
    return object_space_implicit_object->Lazy_Inside_Extended_Levelset(Object_Space_Point(location),Object_Space_Length(contour_value));
}
//#####################################################################
// Function Lazy_Inside_Extended_Levelset_And_Value
//#####################################################################
template<class TV,class TRANSFORM> bool IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Lazy_Inside_Extended_Levelset_And_Value(const TV& location,T& phi_value,const T contour_value) const
{
    bool result=object_space_implicit_object->Lazy_Inside_Extended_Levelset_And_Value(Object_Space_Point(location),phi_value,Object_Space_Length(contour_value));
    phi_value=World_Space_Length(phi_value);
    return result;
}
//#####################################################################
// Function Outside
//#####################################################################
template<class TV,class TRANSFORM> bool IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Outside(const TV& location,const T thickness_over_two) const
{
    return object_space_implicit_object->Outside(Object_Space_Point(location),Object_Space_Length(thickness_over_two));
}
//#####################################################################
// Function Lazy_Outside
//#####################################################################
template<class TV,class TRANSFORM> bool IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Lazy_Outside(const TV& location,const T contour_value) const
{
    return object_space_implicit_object->Lazy_Outside(Object_Space_Point(location),Object_Space_Length(contour_value));
}
//#####################################################################
// Function Lazy_Outside_Extended_Levelset
//#####################################################################
template<class TV,class TRANSFORM> bool IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Lazy_Outside_Extended_Levelset(const TV& location,const T contour_value) const
{
    return object_space_implicit_object->Lazy_Outside_Extended_Levelset(Object_Space_Point(location),Object_Space_Length(contour_value));
}
//#####################################################################
// Function Lazy_Outside_Extended_Levelset_And_Value
//#####################################################################
template<class TV,class TRANSFORM> bool IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Lazy_Outside_Extended_Levelset_And_Value(const TV& location,T& phi_value,const T contour_value) const
{
    bool result=object_space_implicit_object->Lazy_Outside_Extended_Levelset_And_Value(Object_Space_Point(location),phi_value,Object_Space_Length(contour_value));
    phi_value=World_Space_Length(phi_value);
    return result;
}
//#####################################################################
// Function Velocity
//#####################################################################
template<class TV,class TRANSFORM> TV IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Velocity(const TV& location) const
{
    return World_Space_Vector(object_space_implicit_object->Velocity(Object_Space_Point(location)));
}
//#####################################################################
// Function Hessian
//#####################################################################
template<class TV,class TRANSFORM> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Hessian(const TV& X) const
{
    return World_Space_Length_Hessian(object_space_implicit_object->Hessian(Object_Space_Point(X)));
}
//#####################################################################
// Function Integration_Step
//#####################################################################
template<class TV,class TRANSFORM> typename TV::SCALAR IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Integration_Step(const T phi) const
{
    return World_Space_Length(object_space_implicit_object->Integration_Step(Object_Space_Length(phi)));
}
//#####################################################################
// Function Minimum_Cell_Size
//#####################################################################
template<class TV,class TRANSFORM> typename TV::SCALAR IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Minimum_Cell_Size() const
{
    return World_Space_Length(object_space_implicit_object->Minimum_Cell_Size());
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class TV,class TRANSFORM> bool IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Intersection(RAY<TV>& ray,const T thickness) const
{
    RAY<TV> object_space_ray=Object_Space_Ray(ray);
    if(object_space_implicit_object->Intersection(object_space_ray,thickness)){
        ray.semi_infinite=false;
        ray.t_max=World_Space_Length(object_space_ray.t_max);
        ray.aggregate_id=object_space_ray.aggregate_id;
        return true;}
    else return false;
}
//#####################################################################
// Function Closest_Point_On_Boundary
//#####################################################################
template<class TV,class TRANSFORM> TV IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Closest_Point_On_Boundary(const TV& location,const T tolerance,const int max_iterations,T* distance) const
{
    TV result=World_Space_Point(object_space_implicit_object->Closest_Point_On_Boundary(Object_Space_Point(location),Object_Space_Length(tolerance),max_iterations,distance));
    if(distance) *distance=World_Space_Length(*distance);
    return result;
}
//#####################################################################
// Function Name
//#####################################################################
template<class TV,class TRANSFORM> std::string IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Name() const
{
    return Static_Name();
}
//#####################################################################
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<float,1>,FRAME<VECTOR<float,1> > >::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<float,1> >*,bool,FRAME<VECTOR<float,1> > const*);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<float,2>,FRAME<VECTOR<float,2> > >::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<float,2> >*,bool,FRAME<VECTOR<float,2> > const*);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<float,3>,FRAME<VECTOR<float,3> > >::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<float,3> >*,bool,FRAME<VECTOR<float,3> > const*);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<float,1>,float>::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<float,1> >*,bool,VECTOR<float,1> const&,const float);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<float,2>,float>::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<float,2> >*,bool,VECTOR<float,2> const&,const float);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<float,3>,float>::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<float,3> >*,bool,VECTOR<float,3> const&,const float);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<double,1>,FRAME<VECTOR<double,1> > >::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<double,1> >*,bool,FRAME<VECTOR<double,1> > const*);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<double,2>,FRAME<VECTOR<double,2> > >::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<double,2> >*,bool,FRAME<VECTOR<double,2> > const*);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<double,3>,FRAME<VECTOR<double,3> > >::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<double,3> >*,bool,FRAME<VECTOR<double,3> > const*);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<double,1>,double>::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<double,1> >*,bool,VECTOR<double,1> const&,const double);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<double,2>,double>::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<double,2> >*,bool,VECTOR<double,2> const&,const double);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<double,3>,double>::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<double,3> >*,bool,VECTOR<double,3> const&,const double);

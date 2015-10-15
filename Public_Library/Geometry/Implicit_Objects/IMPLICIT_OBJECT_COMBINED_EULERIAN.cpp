//#####################################################################
// Copyright 2007, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
//#####################################################################
// Constructor
//#####################################################################
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_COMBINED_EULERIAN.h>

namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
IMPLICIT_OBJECT_COMBINED_EULERIAN(IMPLICIT_OBJECT<TV>* implicit_object1_input,bool owns_implicit_object1_input,
    IMPLICIT_OBJECT<TV>* implicit_object2_input, bool owns_implicit_object2_input)
    :implicit_object1(implicit_object1_input),implicit_object2(implicit_object2_input),
    owns_implicit_object1(owns_implicit_object1_input),owns_implicit_object2(owns_implicit_object2_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
~IMPLICIT_OBJECT_COMBINED_EULERIAN()
{
    if(owns_implicit_object1) delete implicit_object1;if(owns_implicit_object2) delete implicit_object2;
}
//#####################################################################
// Function Set_Weights
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Set_Weights(T alpha_input)
{
    alpha=alpha_input;
}
//#####################################################################
// Function Update_Box
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Update_Box()
{
    implicit_object1->Update_Box();implicit_object2->Update_Box();box=implicit_object1->box;
    box.Enlarge_Nonempty_Box_To_Include_Points(implicit_object2->box.min_corner,implicit_object2->box.max_corner);
}
//#####################################################################
// Function template<cl
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Update_Minimum_Cell_Size(const int maximum_depth)
{
    implicit_object1->Update_Minimum_Cell_Size(maximum_depth);implicit_object2->Update_Minimum_Cell_Size(maximum_depth);
}
//#####################################################################
// Function Minimum_Cell_Size_Within_Box
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const
{
    return min(implicit_object1->Minimum_Cell_Size_Within_Box(box),implicit_object2->Minimum_Cell_Size_Within_Box(box));
} 
//#####################################################################
// Function Velocity
//#####################################################################
template<class TV> TV IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Velocity(const TV& location) const
{
    return (1-alpha)*implicit_object1->Velocity(location)+alpha*implicit_object2->Velocity(location);
}
//#####################################################################
// Function operator()
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
operator()(const TV& location) const
{
    return Value(location);
}
//#####################################################################
// Function Value
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Value(const TV& location) const
{
    TV V_half=Velocity(location); // half velocity for semi lagrangian
    T phi1=(*implicit_object1)(location-alpha*dt*V_half),phi2=(*implicit_object2)(location+dt*(1-alpha)*V_half); 
    return (1-alpha)*phi1+alpha*phi2;
}
//#####################################################################
// Function Extended_Value
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Extended_Value(const TV& location) const
{
    TV V_half=Velocity(location); // half velocity for semi lagrangian
    T phi1=implicit_object1->Extended_Phi(location-alpha*dt*V_half),phi2=implicit_object2->Extended_Phi(location+dt*(1-alpha)*V_half);
    return (1-alpha)*phi1+alpha*phi2;
}
//#####################################################################
// Function Signed_Distance
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Signed_Distance(const TV& location) const
{
    return (*this)(location);
} // to make this class compatible with the geometry classes
//#####################################################################
// Function Extended_Phi
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Extended_Phi(const TV& location) const
{
    return Extended_Value(location);
}
//#####################################################################
// Function Phi_Secondary
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Phi_Secondary(const TV& location) const
{
    TV V_half=Velocity(location); // half velocity for semi lagrangian
    T phi1=implicit_object1->Phi_Secondary(location-alpha*dt*V_half),phi2=implicit_object2->Phi_Secondary(location+dt*(1-alpha)*V_half); 
    return (1-alpha)*phi1+alpha*phi2;
}
//#####################################################################
// Function Normal
//#####################################################################
template<class TV> TV IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Normal(const TV& location,const int aggregate) const
{
    T dx=Minimum_Cell_Size();TV normal;T one_over_2_dx=(T)1/((T)2*dx);
    for(int i=0;i<TV::m;i++){
        TV offset=dx*TV::Axis_Vector(i);
        normal[i]=(Value(location+offset)-Value(location-offset))*one_over_2_dx;}
    return normal.Normalized();
}
//#####################################################################
// Function Extended_Normal
//#####################################################################
template<class TV> TV IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Extended_Normal(const TV& location,const int aggregate) const
{
    T dx=Minimum_Cell_Size();TV normal;T one_over_2_dx=(T)1/((T)2*dx);
    for(int i=0;i<TV::m;i++){
        TV offset=dx*TV::Axis_Vector(i);
        normal[i]=(Extended_Value(location+offset)-Extended_Value(location-offset))*one_over_2_dx;}
    return normal.Normalized();
}
//#####################################################################
// Function Compute_Normals
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Compute_Normals()
{
} 
//#####################################################################
// Function Compute_Cell_Minimum_And_Maximum
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists)
{
    implicit_object1->Compute_Cell_Minimum_And_Maximum(recompute_if_exists);
    implicit_object2->Compute_Cell_Minimum_And_Maximum(recompute_if_exists);
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Rescale(const T scaling_factor)
{
    implicit_object1->Rescale(scaling_factor);implicit_object2->Rescale(scaling_factor);
}
//#####################################################################
// Function Inflate
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Inflate(const T inflation_distance)
{
    implicit_object1->Inflate(inflation_distance);implicit_object2->Inflate(inflation_distance);
}
//#####################################################################
// Function Lazy_Inside
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Lazy_Inside(const TV& location,const T contour_value) const
{
    return box.Inside(location,-contour_value) && (*this)(location)<=contour_value;
}
//#####################################################################
// Function Lazy_Inside_And_Value
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Lazy_Inside_And_Value(const TV& location,T& phi,const T contour_value) const
{
    if(box.Inside(location,-contour_value)){
        phi=(*this)(location);
        if(phi<=contour_value) return true;}
    return false;
}
//#####################################################################
// Function Lazy_Inside_Extended_Levelset
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Lazy_Inside_Extended_Levelset(const TV& location,const T contour_value) const
{
    return box.Inside(location,-contour_value) && Extended_Value(location)<=contour_value;
}
//#####################################################################
// Function Lazy_Inside_Extended_Levelset_And_Value
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Lazy_Inside_Extended_Levelset_And_Value(const TV& location,T& phi_value,const T contour_value) const
{
    if(box.Inside(location,-contour_value)){
        phi_value=Extended_Phi(location);
        if(phi_value<=contour_value) return true;}
    return false;
}
//#####################################################################
// Function Lazy_Outside
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Lazy_Outside(const TV& location,const T contour_value) const
{
    return !Lazy_Inside(location,contour_value);
}
//#####################################################################
// Function Lazy_Outside_Extended_Levelset
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Lazy_Outside_Extended_Levelset(const TV& location,const T contour_value) const
{
    return !Lazy_Inside_Extended_Levelset(location,contour_value);
}
//#####################################################################
// Function Lazy_Outside_Extended_Levelset_And_Value
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Lazy_Outside_Extended_Levelset_And_Value(const TV& location,T& phi_value,const T contour_value) const
{
    return !Lazy_Inside_Extended_Levelset_And_Value(location,phi_value,contour_value);
}
//#####################################################################
// Function Min_Phi
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Min_Phi() const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Hessian
//#####################################################################
template<class TV> auto IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Hessian(const TV& X) const -> SYMMETRIC_MATRIX<T,TV::m>
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class TV> auto IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Principal_Curvatures(const TV& X) const -> VECTOR<T,d-1>
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Integration_Step
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Integration_Step(const T phi) const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Minimum_Cell_Size
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Minimum_Cell_Size() const
{
    return min(implicit_object1->Minimum_Cell_Size(),implicit_object2->Minimum_Cell_Size());
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Intersection(RAY<TV>& ray,const T thickness) const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Read(TYPED_ISTREAM& input)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>::
Write(TYPED_OSTREAM& output) const
{
    PHYSBAM_FATAL_ERROR();
}
}

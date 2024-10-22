//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_IMPLICIT_OBJECT
//#####################################################################
#include <Core/Matrices/MATRIX_3X3.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Registry/STRUCTURE_REGISTRY.h>
using namespace PhysBAM;
//#####################################################################
template<class TV> LINEAR_INTERPOLATION_UNIFORM<TV,TV> LEVELSET_IMPLICIT_OBJECT<TV>::default_velocity_interpolation;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LEVELSET_IMPLICIT_OBJECT<TV>::
LEVELSET_IMPLICIT_OBJECT(GRID<TV>& grid_input,ARRAY<T,TV_INT>& phi_input)
    :levelset(grid_input,phi_input),V(0),velocity_interpolation(&default_velocity_interpolation),need_destroy_data(false)
{
    Update_Box();Update_Minimum_Cell_Size();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LEVELSET_IMPLICIT_OBJECT<TV>::
~LEVELSET_IMPLICIT_OBJECT()
{
    if(need_destroy_data){delete &levelset.grid;delete &levelset.phi;}
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV> LEVELSET_IMPLICIT_OBJECT<TV>* LEVELSET_IMPLICIT_OBJECT<TV>::
Create()
{
    LEVELSET_IMPLICIT_OBJECT* levelset_implicit_object=new LEVELSET_IMPLICIT_OBJECT(*(new GRID<TV>),*(new ARRAY<T,TV_INT>));
    levelset_implicit_object->need_destroy_data=true;return levelset_implicit_object;
}
//#####################################################################
// Function Update_Minimum_Cell_Size
//#####################################################################
template<class TV> void LEVELSET_IMPLICIT_OBJECT<TV>::
Set_Custom_Secondary_Interpolation(INTERPOLATION_UNIFORM<TV,T>& interpolation)
{
    levelset.Set_Custom_Secondary_Interpolation(interpolation);
}
//#####################################################################
// Function Update_Minimum_Cell_Size
//#####################################################################
template<class TV> void LEVELSET_IMPLICIT_OBJECT<TV>::
Set_Custom_Normal_Interpolation(INTERPOLATION_UNIFORM<TV,TV>& interpolation)
{
    levelset.Set_Custom_Normal_Interpolation(interpolation);
}
//#####################################################################
// Function Update_Minimum_Cell_Size
//#####################################################################
template<class TV> void LEVELSET_IMPLICIT_OBJECT<TV>::
Set_Custom_Velocity_Interpolation(INTERPOLATION_UNIFORM<TV,TV>& interpolation)
{
    velocity_interpolation=&interpolation;
}
//#####################################################################
// Function Update_Box
//#####################################################################
template<class TV> void LEVELSET_IMPLICIT_OBJECT<TV>::
Update_Box()
{
    box=levelset.grid.domain; // workaround bug in gcc 4.2.4 by double assigning - bleurg
    BASE::box=levelset.grid.domain; // BASE:: prefix to work around compiler bug in gcc 4.2.4
}
//#####################################################################
// Function Update_Minimum_Cell_Size
//#####################################################################
template<class TV> void LEVELSET_IMPLICIT_OBJECT<TV>::
Update_Minimum_Cell_Size(const int maximum_depth)
{
    minimum_cell_size=levelset.grid.dX.Min();
}
//#####################################################################
// Function Minimum_Cell_Size_Within_Box
//#####################################################################
template<class TV> typename TV::SCALAR LEVELSET_IMPLICIT_OBJECT<TV>::
Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const
{
    return levelset.grid.dX.Min();
}
//#####################################################################
// Function operator()
//#####################################################################
template<class TV> typename TV::SCALAR LEVELSET_IMPLICIT_OBJECT<TV>::
operator()(const TV& location) const
{
    return levelset.Phi(location);
}
//#####################################################################
// Function Extended_Phi
//#####################################################################
template<class TV> typename TV::SCALAR LEVELSET_IMPLICIT_OBJECT<TV>::
Extended_Phi(const TV& location) const
{
    return levelset.Extended_Phi(location);
}
//#####################################################################
// Function Phi_Secondary
//#####################################################################
template<class TV> typename TV::SCALAR LEVELSET_IMPLICIT_OBJECT<TV>::
Phi_Secondary(const TV& location) const
{
    return levelset.Phi_Secondary(location);
}
//#####################################################################
// Function Normal
//#####################################################################
template<class TV> TV LEVELSET_IMPLICIT_OBJECT<TV>::
Normal(const TV& location,const int aggregate) const
{
    return aggregate==-1?levelset.Normal(location):box.Normal(aggregate);
}
//#####################################################################
// Function Extended_Normal
//#####################################################################
template<class TV> TV LEVELSET_IMPLICIT_OBJECT<TV>::
Extended_Normal(const TV& location,const int aggregate) const
{
    return aggregate==-1?levelset.Extended_Normal(location):box.Normal(aggregate);
}
//#####################################################################
// Function Compute_Normals
//#####################################################################
template<class TV> void LEVELSET_IMPLICIT_OBJECT<TV>::
Compute_Normals()
{
    levelset.Compute_Normals();
}
//#####################################################################
// Function Compute_Cell_Minimum_And_Maximum
//#####################################################################
template<class TV> void LEVELSET_IMPLICIT_OBJECT<TV>::
Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists)
{
    levelset.Compute_Cell_Minimum_And_Maximum(recompute_if_exists);
}
//#####################################################################
// Function Inflate
//#####################################################################
template<class TV> void LEVELSET_IMPLICIT_OBJECT<TV>::
Inflate(const T inflation_distance)
{
    levelset.phi-=inflation_distance;
}
//#####################################################################
// Function Lazy_Inside
//#####################################################################
template<class TV> bool LEVELSET_IMPLICIT_OBJECT<TV>::
Lazy_Inside(const TV& location,const T contour_value) const
{
    return box.Lazy_Inside(location) && levelset.Lazy_Inside(location,contour_value);
}
//#####################################################################
// Function Lazy_Inside_And_Value
//#####################################################################
template<class TV> bool LEVELSET_IMPLICIT_OBJECT<TV>::
Lazy_Inside_And_Value(const TV& location,T& phi_value,const T contour_value) const
{
    return box.Lazy_Inside(location) && levelset.Lazy_Inside_And_Value(location,phi_value,contour_value);
}
//#####################################################################
// Function Lazy_Inside_Extended_Levelset
//#####################################################################
template<class TV> bool LEVELSET_IMPLICIT_OBJECT<TV>::
Lazy_Inside_Extended_Levelset(const TV& unclamped_X,const T contour_value) const
{
    return levelset.Lazy_Inside_Extended_Levelset(unclamped_X,contour_value);
}
//#####################################################################
// Function Lazy_Inside_Extended_Levelset_And_Value
//#####################################################################
template<class TV> bool LEVELSET_IMPLICIT_OBJECT<TV>::
Lazy_Inside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value) const
{
    return levelset.Lazy_Inside_Extended_Levelset_And_Value(unclamped_X,phi_value,contour_value);
}
//#####################################################################
// Function Lazy_Outside
//#####################################################################
template<class TV> bool LEVELSET_IMPLICIT_OBJECT<TV>::
Lazy_Outside(const TV& location,const T contour_value) const
{
    return box.Lazy_Outside(location) || levelset.Lazy_Outside(location,contour_value);
}
//#####################################################################
// Function Lazy_Outside_Extended_Levelset
//#####################################################################
template<class TV> bool LEVELSET_IMPLICIT_OBJECT<TV>::
Lazy_Outside_Extended_Levelset(const TV& unclamped_X,const T contour_value) const
{
    return levelset.Lazy_Outside_Extended_Levelset(unclamped_X,contour_value);
}
//#####################################################################
// Function Lazy_Outside_Extended_Levelset_And_Value
//#####################################################################
template<class TV> bool LEVELSET_IMPLICIT_OBJECT<TV>::
Lazy_Outside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value) const
{
    return levelset.Lazy_Outside_Extended_Levelset_And_Value(unclamped_X,phi_value,contour_value);
}
//#####################################################################
// Function Velocity
//#####################################################################
template<class TV> TV LEVELSET_IMPLICIT_OBJECT<TV>::
Velocity(const TV& location) const
{
    assert(V);return velocity_interpolation->Clamped_To_Array(levelset.grid,*V,location);
}
//#####################################################################
// Function Hessian
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> LEVELSET_IMPLICIT_OBJECT<TV>::
Hessian(const TV& X) const
{
    return levelset.Hessian(X);
}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class TV> typename LEVELSET_IMPLICIT_OBJECT<TV>::T_PRINCIPAL_CURVATURES LEVELSET_IMPLICIT_OBJECT<TV>::
Principal_Curvatures(const TV& X) const
{
    return levelset.Principal_Curvatures(X);
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class TV> void LEVELSET_IMPLICIT_OBJECT<TV>::
Rescale(const T scaling_factor)
{
    PHYSBAM_ASSERT(scaling_factor>0,"Attempting to scale by negative factor");
    levelset.grid.Initialize(levelset.grid.Domain_Indices().Maximum_Corner(),scaling_factor*levelset.grid.domain);
    levelset.phi*=scaling_factor;Update_Box();Update_Minimum_Cell_Size();
}
//#####################################################################
// Function Translate
//#####################################################################
template<class TV> void LEVELSET_IMPLICIT_OBJECT<TV>::
Translate(const TV& translation)
{
    levelset.grid.Initialize(levelset.grid.Domain_Indices().Minimum_Corner(),levelset.grid.domain+translation);
    Update_Box();Update_Minimum_Cell_Size();
}
//#####################################################################
// Function Name
//#####################################################################
template<class TV> std::string LEVELSET_IMPLICIT_OBJECT<TV>::
Static_Name()
{
    return LOG::sprintf("LEVELSET_IMPLICIT_OBJECT<T,VECTOR<T,%d> >",TV::m);
}
//#####################################################################
// Function Extension
//#####################################################################
template<class TV> std::string LEVELSET_IMPLICIT_OBJECT<TV>::
Static_Extension()
{
    return TV::m==2?"phi2d":"phi";
}
//#####################################################################
// Function Integration_Step
//#####################################################################
template<class TV> typename TV::SCALAR LEVELSET_IMPLICIT_OBJECT<TV>::
Integration_Step(const T phi) const
{
    T distance=abs(phi);
    if(distance > 3*minimum_cell_size) return (T).5*distance;
    else if(distance > minimum_cell_size) return (T).25*distance;
    else return (T).1*minimum_cell_size;
}
//#####################################################################
// Function Minimum_Cell_Size
//#####################################################################
template<class TV> typename TV::SCALAR LEVELSET_IMPLICIT_OBJECT<TV>::
Minimum_Cell_Size() const
{
    return minimum_cell_size;
}
//#####################################################################
namespace PhysBAM{
template class LEVELSET_IMPLICIT_OBJECT<VECTOR<float,1> >;
template class LEVELSET_IMPLICIT_OBJECT<VECTOR<float,2> >;
template class LEVELSET_IMPLICIT_OBJECT<VECTOR<float,3> >;
template class LEVELSET_IMPLICIT_OBJECT<VECTOR<double,1> >;
template class LEVELSET_IMPLICIT_OBJECT<VECTOR<double,2> >;
template class LEVELSET_IMPLICIT_OBJECT<VECTOR<double,3> >;
}

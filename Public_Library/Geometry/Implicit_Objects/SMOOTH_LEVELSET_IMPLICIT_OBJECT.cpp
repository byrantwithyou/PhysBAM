//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SMOOTH_LEVELSET_IMPLICIT_OBJECT
//#####################################################################
#include <Core/Matrices/MATRIX_3X3.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Grid_PDE/Interpolation/CUBIC_SPLINE_INTERPOLATION_UNIFORM.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Geometry/Implicit_Objects/SMOOTH_LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Registry/STRUCTURE_REGISTRY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV>::
SMOOTH_LEVELSET_IMPLICIT_OBJECT(GRID<TV>& grid_input,ARRAY<T,TV_INT>& phi_input)
    :BASE(grid_input,phi_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV>::
~SMOOTH_LEVELSET_IMPLICIT_OBJECT()
{
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV> SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV>* SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV>::
Create()
{
    SMOOTH_LEVELSET_IMPLICIT_OBJECT* levelset_implicit_object=new SMOOTH_LEVELSET_IMPLICIT_OBJECT(*(new GRID<TV>),*(new ARRAY<T,TV_INT>));
    levelset_implicit_object->need_destroy_data=true;
    return levelset_implicit_object;
}
//#####################################################################
// Function Normal
//#####################################################################
template<class TV> TV SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV>::
Normal(const TV& location,const int aggregate) const
{
    return levelset.interpolation->Clamped_To_Array_Gradient(levelset.grid,levelset.phi,location);
}
//#####################################################################
// Function Extended_Normal
//#####################################################################
template<class TV> TV SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV>::
Extended_Normal(const TV& location,const int aggregate) const
{
    TV clamped_location(levelset.grid.Clamp(location));
    TV diff=location-clamped_location;
    T magnitude_squared=diff.Magnitude_Squared();
    TV phi_value_grad=levelset.interpolation->Clamped_To_Array_Gradient(levelset.grid,levelset.phi,clamped_location);
    if(!magnitude_squared) return phi_value_grad;
    T phi_value=levelset.interpolation->Clamped_To_Array(levelset.grid,levelset.phi,clamped_location);
    if(phi_value<=0) diff.Normalized();
    for(int i=0;i<TV::m;i++) if(diff(i)) phi_value_grad(i)=0;
    T r=sqrt(magnitude_squared+sqr(phi_value));
    return (diff+phi_value*phi_value_grad)/r;
}
//#####################################################################
// Function Compute_Normals
//#####################################################################
template<class TV> void SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV>::
Compute_Normals()
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Hessian
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV>::
Hessian(const TV& X) const
{
    return levelset.interpolation->Clamped_To_Array_Hessian(levelset.grid,levelset.phi,X);
}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class TV> typename SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV>::T_PRINCIPAL_CURVATURES SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV>::
Principal_Curvatures(const TV& X) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Name
//#####################################################################
template<class TV> std::string SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV>::
Static_Name()
{
    return LOG::sprintf("SMOOTH_LEVELSET_IMPLICIT_OBJECT<T,VECTOR<T,%d> >",TV::dimension);
}
//#####################################################################
// Function Extension
//#####################################################################
template<class TV> std::string SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV>::
Static_Extension()
{
    return TV::dimension==2?"phi2d":"phi";
}
//#####################################################################
namespace PhysBAM{
template class SMOOTH_LEVELSET_IMPLICIT_OBJECT<VECTOR<float,1> >;
template class SMOOTH_LEVELSET_IMPLICIT_OBJECT<VECTOR<float,2> >;
template class SMOOTH_LEVELSET_IMPLICIT_OBJECT<VECTOR<float,3> >;
template class SMOOTH_LEVELSET_IMPLICIT_OBJECT<VECTOR<double,1> >;
template class SMOOTH_LEVELSET_IMPLICIT_OBJECT<VECTOR<double,2> >;
template class SMOOTH_LEVELSET_IMPLICIT_OBJECT<VECTOR<double,3> >;
}

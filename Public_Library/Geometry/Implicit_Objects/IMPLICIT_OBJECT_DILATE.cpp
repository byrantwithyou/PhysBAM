//#####################################################################
// Copyright 2016
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Math_Tools/RANGE.h>
#include <Core/Matrices/MATRIX_3X3.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_DILATE.h>
#include <Geometry/Intersections/RAY_BOX_INTERSECTION.h>
#include <Geometry/Level_Sets/IMPLICIT_OBJECT_ON_A_RAY.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
namespace PhysBAM{
//#####################################################################
// Destructor
//#####################################################################
template<class TV> IMPLICIT_OBJECT_DILATE<TV>::
~IMPLICIT_OBJECT_DILATE()
{
    if(owns_io) delete io;
}
//#####################################################################
// Function Update_Box
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_DILATE<TV>::
Update_Box()
{
    io->Update_Box();
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV> IMPLICIT_OBJECT_DILATE<TV>* IMPLICIT_OBJECT_DILATE<TV>::
Create()
{
    return new IMPLICIT_OBJECT_DILATE(0,0);
}
//#####################################################################
// Function Update_Minimum_Cell_Size
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_DILATE<TV>::
Update_Minimum_Cell_Size(const int maximum_depth)
{
    io->Update_Minimum_Cell_Size(maximum_depth);
}
//#####################################################################
// Function Minimum_Cell_Size_Within_Box
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_DILATE<TV>::
Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const
{
    return io->Minimum_Cell_Size_Within_Box(box);
}
//#####################################################################
// Function operator
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_DILATE<TV>::
operator()(const TV& X) const
{
    return (*io)(X)-dilation;
}
//#####################################################################
// Function Extended_Phi
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_DILATE<TV>::
Extended_Phi(const TV& X) const
{
    return io->Extended_Phi(X)-dilation;
}
//#####################################################################
// Function Phi_Secondary
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_DILATE<TV>::
Phi_Secondary(const TV& X) const
{
    return io->Phi_Secondary(X)-dilation;
}
//#####################################################################
// Function Normal
//#####################################################################
template<class TV> TV IMPLICIT_OBJECT_DILATE<TV>::
Normal(const TV& X,const int aggregate) const
{
    return io->Normal(X,aggregate);
}
//#####################################################################
// Function Extended_Normal
//#####################################################################
template<class TV> TV IMPLICIT_OBJECT_DILATE<TV>::
Extended_Normal(const TV& X,const int aggregate) const
{
    return io->Extended_Normal(X,aggregate);
}
//#####################################################################
// Function Compute_Normals
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_DILATE<TV>::
Compute_Normals()
{
    io->Compute_Normals();
}
//#####################################################################
// Function Compute_Cell_Minimum_And_Maximum
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_DILATE<TV>::
Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists)
{
    io->Compute_Cell_Minimum_And_Maximum(recompute_if_exists);
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_DILATE<TV>::
Rescale(const T scaling_factor)
{
    io->Rescale(scaling_factor);
    dilation*=scaling_factor;
}
//#####################################################################
// Function Translate
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_DILATE<TV>::
Translate(const TV& translation)
{
    io->Translate(translation);
}
//#####################################################################
// Function Inflate
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_DILATE<TV>::
Inflate(const T inflation_distance)
{
    io->Inflate(inflation_distance);
}
//#####################################################################
// Function Lazy_Inside
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_DILATE<TV>::
Lazy_Inside(const TV& X,const T contour_value) const
{
    return this->Inside(X,contour_value);
}
//#####################################################################
// Function Lazy_Inside_And_Value
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_DILATE<TV>::
Lazy_Inside_And_Value(const TV& X,T& phi_value,const T contour_value) const
{
    bool b=io->Lazy_Inside_And_Value(X,phi_value,contour_value-dilation);
    phi_value-=dilation;
    return b;
}
//#####################################################################
// Function Lazy_Inside_Extended_Levelset
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_DILATE<TV>::
Lazy_Inside_Extended_Levelset(const TV& X,const T contour_value) const
{
    return io->Lazy_Inside_Extended_Levelset(X,contour_value-dilation);
}
//#####################################################################
// Function Lazy_Inside_Extended_Levelset_And_Value
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_DILATE<TV>::
Lazy_Inside_Extended_Levelset_And_Value(const TV& X,T& phi_value,const T contour_value) const
{
    bool b=io->Lazy_Inside_Extended_Levelset_And_Value(X,phi_value,contour_value-dilation);
    phi_value-=dilation;
    return b;
}
//#####################################################################
// Function Lazy_Outside
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_DILATE<TV>::
Lazy_Outside(const TV& X,const T contour_value) const
{
    return this->Inside(X,contour_value);
}
//#####################################################################
// Function Lazy_Outside_Extended_Levelset
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_DILATE<TV>::
Lazy_Outside_Extended_Levelset(const TV& X,const T contour_value) const
{
    return io->Lazy_Outside_Extended_Levelset(X,contour_value-dilation);
}
//#####################################################################
// Function Lazy_Outside_Extended_Levelset_And_Value
//#####################################################################
template<class TV> bool IMPLICIT_OBJECT_DILATE<TV>::
Lazy_Outside_Extended_Levelset_And_Value(const TV& X,T& phi_value,const T contour_value) const
{
    bool b=io->Lazy_Outside_Extended_Levelset_And_Value(X,phi_value,contour_value-dilation);
    phi_value-=dilation;
    return b;
}
//#####################################################################
// Function Velocity
//#####################################################################
template<class TV> TV IMPLICIT_OBJECT_DILATE<TV>::
Velocity(const TV& X) const
{
    return io->Velocity(X);
}
//#####################################################################
// Function Velocity
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> IMPLICIT_OBJECT_DILATE<TV>::
Hessian(const TV& X) const
{
    return io->Hessian(X);
}
//#####################################################################
// Function Velocity
//#####################################################################
template<class TV> auto IMPLICIT_OBJECT_DILATE<TV>::
Principal_Curvatures(const TV& X) const -> T_CURVATURES
{
    return io->Principal_Curvatures(X);
}
//#####################################################################
// Function Integration_Step
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_DILATE<TV>::
Integration_Step(const T phi) const
{
    return io->Integration_Step(phi);
}
//#####################################################################
// Function Minimum_Cell_Size
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_DILATE<TV>::
Minimum_Cell_Size() const
{
    return io->Minimum_Cell_Size();
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_DILATE<TV>::
Read(TYPED_ISTREAM& input)
{
    Read_Binary(input,dilation);
    if(io) io->Read_Structure(input);
    else io=dynamic_cast<IMPLICIT_OBJECT<TV>*>(Create_Structure(input));
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_DILATE<TV>::
Write(TYPED_OSTREAM& output) const
{
    Write_Binary(output,dilation);
    io->Write_Structure(output);
}
//#####################################################################
// Function Name
//#####################################################################
template<class TV> std::string IMPLICIT_OBJECT_DILATE<TV>::
Static_Name()
{
    return LOG::sprintf("IMPLICIT_OBJECT_DILATE<T,VECTOR<T,%d> >",TV::m);
}
//#####################################################################
// Function Extension
//#####################################################################
template<class TV> std::string IMPLICIT_OBJECT_DILATE<TV>::
Static_Extension()
{
    return TV::m==2?"dilate_phi2d":"dilate_phi";
}
template class IMPLICIT_OBJECT_DILATE<VECTOR<float,1> >;
template class IMPLICIT_OBJECT_DILATE<VECTOR<float,2> >;
template class IMPLICIT_OBJECT_DILATE<VECTOR<float,3> >;
template class IMPLICIT_OBJECT_DILATE<VECTOR<double,1> >;
template class IMPLICIT_OBJECT_DILATE<VECTOR<double,2> >;
template class IMPLICIT_OBJECT_DILATE<VECTOR<double,3> >;

}

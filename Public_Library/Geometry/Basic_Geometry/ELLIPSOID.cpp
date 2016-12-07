//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Neil Molino, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Math_Tools/max.h>
#include <Core/Math_Tools/min.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Geometry/Basic_Geometry/ELLIPSOID.h>
using namespace PhysBAM;
//#####################################################################
// Function Normal
//#####################################################################
template<class TV> TV ELLIPSOID<TV>::
Normal(const TV& location) const   
{
    return orientation.Rotate(sqr(radii).Inverse_Times(orientation.Inverse_Rotate(location-center))).Normalized();
}
//#####################################################################
// Function Inside
//#####################################################################
template<class TV> bool ELLIPSOID<TV>::
Inside(const TV& location,const T thickness_over_two) const
{
    TV scaled_offset=radii.Inverse_Times(orientation.Inverse_Rotate(location-center));
    return scaled_offset.Magnitude_Squared() <= sqr(1-thickness_over_two);
}
//#####################################################################
// Function Outside
//#####################################################################
template<class TV> bool ELLIPSOID<TV>::
Outside(const TV& location,const T thickness_over_two) const
{
    TV scaled_offset=radii.Inverse_Times(orientation.Inverse_Rotate(location-center));
    return scaled_offset.Magnitude_Squared() >= sqr(1+thickness_over_two);
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class TV> bool ELLIPSOID<TV>::
Boundary(const TV& location,const T thickness_over_two) const
{
    return !Inside(location,thickness_over_two) && !Outside(location,thickness_over_two);
}
//#####################################################################
// Function Approximate_Surface
//#####################################################################
// not necessarily the closest point on the surface, but approximates it in some sense
template<class TV> TV ELLIPSOID<TV>::
Approximate_Surface(const TV& location) const 
{
    TV scaled_offset=radii.Inverse_Times(orientation.Inverse_Rotate(location-center));
    return orientation.Rotate(radii*scaled_offset.Normalized())+center;
}
//#####################################################################
// Function Approximate_Surface
//#####################################################################
template<class TV> typename TV::SCALAR ELLIPSOID<TV>::
Approximate_Signed_Distance(const TV& location) const       
{
    TV offset=orientation.Inverse_Rotate(location-center),scaled_offset=radii.Inverse_Times(offset);
    T scaled_magnitude=scaled_offset.Normalize();
    TV closest_offset=radii*scaled_offset;
    T distance=(closest_offset-offset).Magnitude();
    return scaled_magnitude<1?-distance:distance;
}
//#####################################################################
// Function Calculate_Approximate_Signed_Distance
//#####################################################################
// better to reinitialize after this because uses the Approximate_Signed_Distance() function
template<class TV> void ELLIPSOID<TV>::
Calculate_Approximate_Signed_Distance(const GRID<TV>& grid,ARRAY<T,TV_INT>& phi) const
{  
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next())
        phi(iterator.Cell_Index())=Approximate_Signed_Distance(iterator.Location());
}
//#####################################################################
// Function Covariance_Ellipsoid
//#####################################################################
template<class TV> template<class T_ARRAY_TV> ELLIPSOID<TV> ELLIPSOID<TV>::
Covariance_Ellipsoid(const T_ARRAY_TV& points)
{
    TV average=points.Average();
    SYMMETRIC_MATRIX<T,3> covariance;DIAGONAL_MATRIX<T,3> eigenvalues;MATRIX<T,3> eigenvectors;
    for(int p=0;p<points.m;p++){
        TV variance=points(p)-average;
        covariance.x00+=variance.x*variance.x;covariance.x10+=variance.y*variance.x;covariance.x20+=variance.z*variance.x;
        covariance.x11+=variance.y*variance.y;covariance.x21+=variance.z*variance.y;covariance.x22+=variance.z*variance.z;}
    covariance/=points.m-(T)1;covariance.Fast_Solve_Eigenproblem(eigenvalues,eigenvectors);
    return ELLIPSOID<TV>(average,eigenvalues.Sqrt(),eigenvectors);
}
//#####################################################################
namespace PhysBAM{
template class ELLIPSOID<VECTOR<double,2> >;
template class ELLIPSOID<VECTOR<double,3> >;
template class ELLIPSOID<VECTOR<float,2> >;
template class ELLIPSOID<VECTOR<float,3> >;
template ELLIPSOID<VECTOR<double,3> > ELLIPSOID<VECTOR<double,3> >::Covariance_Ellipsoid<ARRAY<VECTOR<double,3>,int> >(ARRAY<VECTOR<double,3>,int> const&);
template ELLIPSOID<VECTOR<float,3> > ELLIPSOID<VECTOR<float,3> >::Covariance_Ellipsoid<ARRAY<VECTOR<float,3>,int> >(ARRAY<VECTOR<float,3>,int> const&);
}

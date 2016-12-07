//#####################################################################
// Copyright 2002-2009, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ELLIPSOID
//#####################################################################
#ifndef __ELLIPSOID__
#define __ELLIPSOID__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Matrices/MATRIX_3X3.h>
#include <Core/Matrices/ROTATION.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class TV>
class ELLIPSOID
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    enum WORKAROUND {d=TV::m};
public:
    TV center;
    DIAGONAL_MATRIX<T,d> radii;
    ROTATION<TV> orientation;

    ELLIPSOID()
        :radii(TV::All_Ones_Vector())
    {}

    ELLIPSOID(const TV& center,const DIAGONAL_MATRIX<T,d>& radii,const ROTATION<TV>& orientation=ROTATION<TV>())
        :center(center),radii(radii),orientation(orientation.Normalized())
    {}

    ELLIPSOID(const TV& center,const MATRIX<T,TV::m>& scaled_axes) // assumes axes are orthogonal
        :center(center)
    {
        MATRIX<T,TV::m> axis_normalized(scaled_axes);
        radii.x=axis_normalized.Normalize_Columns();
        orientation=ROTATION<TV>(axis_normalized).Normalized();
    }

    ELLIPSOID(const TV& center,const DIAGONAL_MATRIX<T,TV::m>& radii,const MATRIX<T,TV::m>& axes) // assumes axes are orthonormal
        :center(center),radii(radii)
    {
        MATRIX<T,TV::m> axis_normalized(axes);
        axis_normalized.Normalize_Columns();
        orientation=ROTATION<TV>(axis_normalized).Normalized();
    }

    ORIENTED_BOX<TV> Oriented_Bounding_Box() const
    {MATRIX<T,d> axes(orientation.Rotation_Matrix()*radii);return ORIENTED_BOX<TV>(center-axes.Column_Sum(),(T)2*axes);}

    RANGE<TV> Bounding_Box() const
    {return Oriented_Bounding_Box().Axis_Aligned_Bounding_Box();}

    SYMMETRIC_MATRIX<T,d> Metric_Tensor() const
    {return (orientation.Rotation_Matrix()*radii.Inverse()).Outer_Product_Matrix();}

    T Lipschitz_Constant() const
    {return 1/radii.Min();}

    T Volume() const
    {return (T)unit_sphere_size[d]*radii.Determinant();}

//#####################################################################
    TV Normal(const TV& location) const;
    bool Inside(const TV& location,const T thickness_over_two=0) const;
    bool Outside(const TV& location,const T thickness_over_two=0) const;
    bool Boundary(const TV& location,const T thickness_over_two) const;
    TV Approximate_Surface(const TV& location) const; // not necessarily the closest point on the surface, but approximates it
    T Approximate_Signed_Distance(const TV& location) const; // not the true signed distance, but has correct inside/outside
    void Calculate_Approximate_Signed_Distance(const GRID<TV>& grid,ARRAY<T,TV_INT>& phi) const; // better to reinitialize after this because uses approx. surface
    template<class T_ARRAY_TV> static ELLIPSOID<TV> Covariance_Ellipsoid(const T_ARRAY_TV& points);
//#####################################################################
};
}
#endif

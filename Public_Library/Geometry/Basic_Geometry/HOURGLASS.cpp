//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Math_Tools/constants.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/ROTATION.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/HOURGLASS.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> HOURGLASS<VECTOR<T,2> >::
HOURGLASS(const TV& axis,const TV& center,T bulb_radius,T neck_radius,T height,T neck_width)
    :center(center)
{
    TV ax=axis.Normalized();
    R=ROTATION<TV>::From_Rotated_Vector(ax,TV(0,1)).Rotation_Matrix();
    T c1=neck_width/2+neck_radius;
    r1=neck_radius;
    T c0=height/2-bulb_radius;
    r0=bulb_radius;
    T r=r0+r1;
    T b=c0*c0+c1*c1;
    T a=b-r*r;
    PHYSBAM_ASSERT(a>0);
    T s=sqrt(a);
    T d=(s*c0+c1*r)/b;
    T x0=r*d;
    T x1=s*d;
    T z=c1*(c0-x1)/x0;
    v0=TV(x0-c1,x1);
    C0=TV(0,c0);
    C1=TV(c1,0);
    N=TV(0,z)-C1;
    N.Normalize();
    o=C1.Dot(v0);
    q=C0.Dot(v0);

    TV off=ax*(height/2-bulb_radius);
    TV w_c0=center+off;
    TV w_c1=center-off;
    RANGE<TV> b0(w_c0-bulb_radius,w_c0+bulb_radius);
    RANGE<TV> b1(w_c1-bulb_radius,w_c1+bulb_radius);
    bounding_box=b0.Unite(b1);
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> HOURGLASS<VECTOR<T,3> >::
HOURGLASS(const TV& axis,const TV& center,T bulb_radius,T neck_radius,T height,T neck_width)
    :axis(axis.Normalized()),center(center),hg2(VECTOR<T,2>(0,1),VECTOR<T,2>(),bulb_radius,neck_radius,height,neck_width)
{
    TV oa=axis.Unit_Orthogonal_Vector();
    M.Set_Row(0,oa);
    M.Set_Row(1,axis.Cross(oa));
    TV off=axis*(height/2-bulb_radius);
    TV w_c0=center+off;
    TV w_c1=center-off;
    RANGE<TV> b0(w_c0-bulb_radius,w_c0+bulb_radius);
    RANGE<TV> b1(w_c1-bulb_radius,w_c1+bulb_radius);
    bounding_box=b0.Unite(b1);
}
template class HOURGLASS<VECTOR<float,2> >;
template class HOURGLASS<VECTOR<double,2> >;
template class HOURGLASS<VECTOR<float,3> >;
template class HOURGLASS<VECTOR<double,3> >;
}

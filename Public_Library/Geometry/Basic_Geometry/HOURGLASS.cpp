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
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> HOURGLASS<VECTOR<T,2> >::
HOURGLASS(const TV& axis,const TV& center,T bulb_radius,T neck_radius,T height,T neck_width)
    :axis(axis.Normalized()),center(center)
{
    R=ROTATION<TV>::From_Rotated_Vector(axis,TV(0,1)).Rotation_Matrix();
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
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,2> > HOURGLASS<VECTOR<T,2> >::
Bounding_Box() const
{
    TV off=axis*C0.y;
    TV w_c0=center+off;
    TV w_c1=center-off;
    RANGE<TV> b0(w_c0-r0,w_c0+r0);
    RANGE<TV> b1(w_c1-r0,w_c1+r0);
    return b0.Unite(b1);
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
namespace TESSELLATION{
//#####################################################################
// Function Tessellate_Boundary
//#####################################################################
template<class T> SEGMENTED_CURVE_2D<T>*
Tessellate_Boundary(const HOURGLASS<VECTOR<T,2> >& h,int axis_div)
{
    typedef VECTOR<T,2> TV;
    SEGMENTED_CURVE_2D<T>& sc=*SEGMENTED_CURVE_2D<T>::Create();
    sc.particles.Add_Elements(axis_div*2);

    TV Y0=h.C0-h.r0*h.N;
    TV Y1=Y0-h.v0;

    T t1=acos(-h.N.y),len1=t1*h.r0;
//    X=C0+r0*TV(sin(t/r0),cos(t/r0));
//    t = 0 .. len1

    T s=h.v0.Magnitude();
//    X=Y0-v0/s*t
//    t = 0 .. s

    T t2=3*pi/2-t1,len3=(pi-t2)*h.r1;

    T len=2*(len1+s+len3);

    for(int i=0;i<axis_div;i++){
        TV X;
        T t=i*len/axis_div;
        T scale=1;
        if(t>len/2){t=len-t;scale=-1;}
        if(t<len1) X=h.C0+h.r0*TV(sin(t/h.r0),cos(t/h.r0));
        else if(t-len1<s) X=Y0-h.v0/s*(t-len1);
        else X=h.C1+h.r1*TV(cos((t-len1-s)/h.r1+t2),sin((t-len1-s)/h.r1+t2));
        X.y*=scale;
        sc.particles.X(i)=X;
        sc.particles.X(i+axis_div)=-X;}
    for(int i=0;i<sc.particles.X.m;i++)
        sc.particles.X(i)=h.center+h.R.Transpose_Times(sc.particles.X(i));
    sc.mesh.Initialize_Straight_Mesh(sc.particles.X.m,true);
    return &sc;
}
//#####################################################################
// Function Tessellate_Boundary
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>*
Tessellate_Boundary(const HOURGLASS<VECTOR<T,3> >& h,int axis_div,int circ_div)
{
    typedef VECTOR<T,3> TV;
    TRIANGULATED_SURFACE<T>& ts=*TRIANGULATED_SURFACE<T>::Create();
    SEGMENTED_CURVE_2D<T>* sc=Tessellate_Boundary(h.hg2,axis_div);
    ts.mesh.Initialize_Cylinder_Mesh(axis_div-1,circ_div,true);
    ts.particles.Add_Elements((axis_div-1)*circ_div+2);
    T dtheta=(T)pi*2/circ_div;
    for(int i=0,p=0;i<axis_div-1;i++)
        for(int j=0;j<circ_div;j++){
        T theta=j*dtheta;
        VECTOR<T,2> Y(sc->particles.X(i+1));
        ts.particles.X(p++)=h.center+h.axis*Y.y+Y.x*h.M.Transpose_Times(VECTOR<T,2>(sin(theta),cos(theta)));}
    VECTOR<T,2> Y0(sc->particles.X(0));
    VECTOR<T,2> Y1(sc->particles.X(axis_div));
    ts.particles.X(ts.particles.X.m-2)=h.center+h.axis*Y0.y;
    ts.particles.X(ts.particles.X.m-1)=h.center+h.axis*Y1.y;
    delete sc;
    return &ts;
}
template SEGMENTED_CURVE_2D<double>* Tessellate_Boundary<double>(HOURGLASS<VECTOR<double,2> > const&,int);
template SEGMENTED_CURVE_2D<float>* Tessellate_Boundary<float>(HOURGLASS<VECTOR<float,2> > const&,int);
template TRIANGULATED_SURFACE<double>* Tessellate_Boundary<double>(HOURGLASS<VECTOR<double, 3> > const&, int, int);
template TRIANGULATED_SURFACE<float>* Tessellate_Boundary<float>(HOURGLASS<VECTOR<float, 3> > const&, int, int);
}
template class HOURGLASS<VECTOR<float,2> >;
template class HOURGLASS<VECTOR<double,2> >;
template class HOURGLASS<VECTOR<float,3> >;
template class HOURGLASS<VECTOR<double,3> >;
}

//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Avi Robinson-Mosher, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPHERE
//##################################################################### 
#include <Core/Math_Tools/choice.h>
#include <Core/Math_Tools/cube.h>
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
using namespace PhysBAM;
//#####################################################################
// Function Sector_Volumes
//#####################################################################
// Also see http://mathworld.wolfram.com/CircularSegment.html for fast approximate formulas
template<class T> void Sector_Volumes_Helper(const SPHERE<VECTOR<T,2> >& circle,const VECTOR<T,2>& origin,T volumes[4],const T thickness_over_two)
{
    typedef VECTOR<T,2> TV;
    T radius=circle.radius;
    TV center=circle.center;

    // sectors in usual order (left to right, bottom to top)
    if(circle.Inside(origin,thickness_over_two)){
        T radius_squared=sqr(radius),y_discriminant=sqrt(radius_squared-sqr(origin.y-center.y)),x_discriminant=sqrt(radius_squared-sqr(origin.x-center.x));
        T x0=center.x+y_discriminant,x1=center.x-y_discriminant,y0=center.y+x_discriminant,y1=center.y-x_discriminant; // intersection points
        T xl1=x0-origin.x,xl2=origin.x-x1,yl1=y0-origin.y,yl2=origin.y-y1; // the lengths of the axes inside the circle
        T ry_times_x_discriminant=(center.y-origin.y)*x_discriminant,rx_times_y_discriminant=(center.x-origin.x)*y_discriminant;
        volumes[0]=(T).5*xl2*yl2+circle.Circular_Segment_Area(max((T)0,radius-sqrt((T).5*(radius_squared+ry_times_x_discriminant+rx_times_y_discriminant))));
        volumes[1]=(T).5*xl1*yl2+circle.Circular_Segment_Area(max((T)0,radius-sqrt((T).5*(radius_squared+ry_times_x_discriminant-rx_times_y_discriminant))));
        volumes[2]=(T).5*xl2*yl1+circle.Circular_Segment_Area(max((T)0,radius-sqrt((T).5*(radius_squared-ry_times_x_discriminant+rx_times_y_discriminant))));
        volumes[3]=(T).5*xl1*yl1+circle.Circular_Segment_Area(max((T)0,radius-sqrt((T).5*(radius_squared-ry_times_x_discriminant-rx_times_y_discriminant))));}
    else if(circle.Bounding_Box().Inside(origin,thickness_over_two)){
        T horizontal_area=circle.Circular_Segment_Area(max((T)0,radius-abs(center.x-origin.x))),vertical_area=circle.Circular_Segment_Area(max((T)0,radius-abs(center.y-origin.y)));
        T remaining_area=circle.Size()-horizontal_area-vertical_area;
        if(origin.x>center.x){
            if(origin.y>center.y){volumes[0]=remaining_area;volumes[1]=horizontal_area;volumes[2]=vertical_area;volumes[3]=0;}
            else{volumes[0]=vertical_area;volumes[1]=0;volumes[2]=remaining_area;volumes[3]=horizontal_area;}}
        else{
            if(origin.y>center.y){volumes[0]=horizontal_area;volumes[1]=remaining_area;volumes[2]=0;volumes[3]=vertical_area;}
            else{volumes[0]=0;volumes[1]=vertical_area;volumes[2]=horizontal_area;volumes[3]=remaining_area;}}}
    else if(circle.Bounding_Box().Remove_Dimension(1).Inside(VECTOR<T,1>(origin.x),thickness_over_two)){ // cuts through vertically
        T left,right;
        T h=origin.x-center.x;
        if(h>0){right=circle.Circular_Segment_Area(max((T)0,radius-h));left=circle.Size()-right;}
        else{left=circle.Circular_Segment_Area(max((T)0,radius+h));right=circle.Size()-left;}
        if(origin.y>center.y){volumes[0]=left;volumes[1]=right;volumes[2]=volumes[3]=0;}
        else{volumes[2]=left;volumes[3]=right;volumes[0]=volumes[1]=0;}}
    else if(circle.Bounding_Box().Remove_Dimension(0).Inside(VECTOR<T,1>(origin.y),thickness_over_two)){
        T top,bottom;
        T h=origin.y-center.y;
        if(h>0){top=circle.Circular_Segment_Area(max((T)0,radius-h));bottom=circle.Size()-top;}
        else{bottom=circle.Circular_Segment_Area(max((T)0,radius+h));top=circle.Size()-bottom;}
        if(origin.x>center.x){volumes[0]=bottom;volumes[2]=top;volumes[1]=volumes[3]=0;}
        else{volumes[1]=bottom;volumes[3]=top;volumes[0]=volumes[2]=0;}}
    else{
        for(int i=0;i<4;i++) volumes[i]=0;
        volumes[(origin.y>center.y?0:2)+(origin.x>center.x?0:1)]=circle.Size();}
}
template<class T> void Sector_Volumes_Helper(const SPHERE<VECTOR<T,3> >& sphere,const VECTOR<T,3>& origin,T volumes[8],const T thickness_over_two)
{
    typedef VECTOR<T,3> TV;
    // TODO: this is temporary; gives box sectors, not sphere sectors
    const RANGE<TV> box=sphere.Bounding_Box();
    TV positive_lengths,max_corner=box.Maximum_Corner(),edge_lengths=box.Edge_Lengths();
    for(int i=0;i<3;i++) positive_lengths(i)=clamp(max_corner(i)-origin(i),(T)0,edge_lengths(i));
    for(int i=0;i<8;i++){volumes[i]=1;for(int j=0;j<3;j++) volumes[i]*=(i&(1<<j))?positive_lengths(j):edge_lengths(j)-positive_lengths(j);}
}
template<class TV> void SPHERE<TV>::
Sector_Volumes(const TV& origin,T volumes[1<<d],const T thickness_over_two) const
{
    Sector_Volumes_Helper(*this,origin,volumes,thickness_over_two);
}
//#####################################################################
// Function Octant_Volume_Helper
//#####################################################################
template<class T> T Octant_Volume_Helper(T x,T x1,T y,T y1)
{
    T s2=max(1-x1-y1,(T)0),s=sqrt(s2),twoxys=2*x*y*s;return ((T)1/6)*(atan2(-twoxys,s2-x1*y1)+(3-y1)*y*(atan2(x,s)-(T)pi/4)+(3-x1)*x*(atan2(y,s)-(T)pi/4)+twoxys);
}
//#####################################################################
// Function Octant_Volume_Internal
//#####################################################################
template<class T> T Octant_Volume_Internal(const VECTOR<T,3>& p)
{
    T x2=sqr(p.x),y2=sqr(p.y),z2=sqr(p.z);return Octant_Volume_Helper(p.x,x2,p.y,y2)+Octant_Volume_Helper(p.x,x2,p.z,z2)+Octant_Volume_Helper(p.y,y2,p.z,z2)-p.x*p.y*p.z+pi/6;
}
//#####################################################################
// Function Octant_Volume_Wedge
//#####################################################################
template<class T> T Octant_Volume_Wedge(T x,T y)
{
    T x2=sqr(x),y2=sqr(y);return 2*(Octant_Volume_Helper(x,x2,y,y2)-((T)pi/24)*((3-x2)*x+(3-y2)*y)+pi/6);
}
//#####################################################################
// Function Octant_Volume_Cap
//#####################################################################
template<class T> T Octant_Volume_Cap(T z)
{
    return ((T)pi/3)*sqr(z-1)*(z+2);
}
//#####################################################################
// Function Octant_Volume
//#####################################################################
template<class TV> typename TV::SCALAR SPHERE<TV>::
Octant_Volume(const VECTOR<T,3>& min_corner) const
{
    STATIC_ASSERT(d==3);
    TV p((min_corner-center)/radius);T r3=cube(radius);
    exchange_sort(p.x,p.y,p.z); // x <= y <= z
    // Trivial & cheap cases (none, all, cap)
    if(p.z>=1) return 0;
    if(p.z<=-1) return (T)pi*4/3*r3;
    if(p.y<=-1) return Octant_Volume_Cap(p.z)*r3;
    // Inside
    if(p.Magnitude_Squared()<=1) return Octant_Volume_Internal(p)*r3;
    // all coordinates positive, not inside - no intersection
    if(p.x>=0) return 0;
    // two coordinates positive, not inside - wedge or no intersection
    if(p.y>=0){if(sqr(p.y)+sqr(p.z)>=1) return 0;return Octant_Volume_Wedge(p.y,p.z)*r3;}
    // one coordinates positive - cap with zero, one, or two wedges removed
    T vol=Octant_Volume_Cap(p.z)*r3;
    if(p.z>=0){
        T rz=1-sqr(p.z),ryz=rz-sqr(p.y);if(ryz<=0) return vol;
        T rxz=rz-sqr(p.x);vol-=Octant_Volume_Wedge(-p.y,p.z)*r3;if(rxz<=0) return vol;
        return vol-Octant_Volume_Wedge(-p.x,p.z)*r3;}
    // all coordinates negative; no octant edges intersect sphere - sphere with zero, one, or two caps removed
    vol-=Octant_Volume_Cap(-p.y)*r3;T ry=1-sqr(p.y),ryz=ry-sqr(p.z);if(ryz<=0 && p.x<=-1) return vol;
    if(p.x>-1) vol-=Octant_Volume_Cap(-p.x)*r3;
    if(ryz<=0) return vol;
    // all coordinates negative; two or three caps removed, one, two, or three wedges added back in to prevent duplicate removal
    vol+=Octant_Volume_Wedge(-p.y,-p.z)*r3;T rxz=1-sqr(p.x)-sqr(p.z);if(rxz<=0) return vol;
    vol+=Octant_Volume_Wedge(-p.x,-p.z)*r3;T rxy=ry-sqr(p.y);if(rxy<=0) return vol;
    return vol+Octant_Volume_Wedge(-p.x,-p.y)*r3;
}
//#####################################################################
// Function Name
//#####################################################################
template<class TV> std::string SPHERE<TV>::
Name()
{
    if(TV::m==2) return "CIRCLE<T>";
    else if(TV::m==3) return "SPHERE<T>";
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
#define P(...) __VA_ARGS__
#define INSTANTIATION_HELPER_T_TV(T,TV,d) \
    template void SPHERE<TV>::Sector_Volumes(const TV& origin,T volumes[1<<d],const T thickness_over_two) const; \
    template std::string SPHERE<TV>::Name();
#define INSTANTIATION_HELPER(T,d) INSTANTIATION_HELPER_T_TV(T,P(VECTOR<T,d> ),d)
INSTANTIATION_HELPER(float,2)
INSTANTIATION_HELPER(float,3)
template float SPHERE<VECTOR<float,3> >::Octant_Volume(const VECTOR<float,3>& min_corner) const;
INSTANTIATION_HELPER(double,2)
INSTANTIATION_HELPER(double,3)
template double SPHERE<VECTOR<double,3> >::Octant_Volume(const VECTOR<double,3>& min_corner) const;

//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Neil Molino, Avi Robinson-Mosher, Andrew Selle, Eftychios Sifakis, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_3D
//##################################################################### 
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <Tools/Math_Tools/INTERVAL.h>
#include <Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <Tools/Polynomials/CUBIC.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Function Change_Size
//#####################################################################
// Enlarges the triangle by pushing out the triangle edges by distance 'delta' orthogonally to the edges.
// This keeps the incenter fixed.  If the triangle is degenerate, it will not be changed.
template<class T> void TRIANGLE_3D<T>::
Change_Size(const T delta)
{
    TV edge_lengths((X.z-X.y).Magnitude(),(X.x-X.z).Magnitude(),(X.y-X.x).Magnitude());
    T perimeter=edge_lengths.Sum(),area=Area();
    if(!perimeter || !area) return; // don't know which direction to enlarge a degenerate triangle, so do nothing
    T scale=1+delta*(T).5*perimeter/area;
    TV incenter=Point_From_Barycentric_Coordinates(edge_lengths/perimeter);
    X.x=incenter+(X.x-incenter)*scale;X.y=incenter+(X.y-incenter)*scale;X.z=incenter+(X.z-incenter)*scale;
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool TRIANGLE_3D<T>::
Inside(const TV& point,const T thickness_over_two) const
{
    return Point_Inside_Triangle(point,thickness_over_two);
}
//#####################################################################
// Function Point_Inside_Triangle
//#####################################################################
template<class T> bool TRIANGLE_3D<T>::
Point_Inside_Triangle(const TV& point,const T thickness_over_2) const
{
    return Plane().Boundary(point,thickness_over_2) && Planar_Point_Inside_Triangle(point,thickness_over_2);
}
//#####################################################################
// Function Planar_Point_Inside_Triangle
//#####################################################################
template<class T> bool TRIANGLE_3D<T>::
Planar_Point_Inside_Triangle(const TV& point,const T thickness_over_2) const
{
    TV normal=Normal();
    PLANE<T> edge_plane(TV::Cross_Product(X.y-X.x,normal).Normalized(),X.x);
    if(edge_plane.Outside(point,thickness_over_2)) return false;
    edge_plane.normal=TV::Cross_Product(X.x-X.z,normal).Normalized();
    if(edge_plane.Outside(point,thickness_over_2)) return false;
    edge_plane.normal=TV::Cross_Product(X.z-X.y,normal).Normalized();
    edge_plane.x0=X.y;
    if(edge_plane.Outside(point,thickness_over_2)) return false;
    return true;
}
//#####################################################################
// Function Lazy_Planar_Point_Inside_Triangle
//#####################################################################
template<class T> bool TRIANGLE_3D<T>::
Lazy_Planar_Point_Inside_Triangle(const TV& point) const
{
    TV normal=Normal();
    TV edge_normal_1=TV::Cross_Product(X.y-X.x,normal),point_minus_x1=point-X.x;
    if(TV::Dot_Product(edge_normal_1,point_minus_x1) > 0) return false;
    TV edge_normal_2=TV::Cross_Product(X.x-X.z,normal);
    if(TV::Dot_Product(edge_normal_2,point_minus_x1) > 0) return false;
    edge_normal_1+=edge_normal_2;
    if(TV::Dot_Product(edge_normal_1,point-X.y) < 0) return false; // this equals X.y-X.z (== -edge_normal_3)
    return true;
}
//#####################################################################
// Function Minimum_Edge_Length
//#####################################################################
template<class T> T TRIANGLE_3D<T>::
Minimum_Edge_Length() const
{
    return min((X.y-X.x).Magnitude(),(X.z-X.y).Magnitude(),(X.x-X.z).Magnitude());
}
//#####################################################################
// Function Maximum_Edge_Length
//#####################################################################
template<class T> T TRIANGLE_3D<T>::
Maximum_Edge_Length() const
{      
    return max((X.y-X.x).Magnitude(),(X.z-X.y).Magnitude(),(X.x-X.z).Magnitude());
}
//#####################################################################
// Function Region
//#####################################################################
// returns the region one is near in priority of 1=vertex, 2=edge, 3=face based on distance, region_id differentiates which point or edge
template<class T> int TRIANGLE_3D<T>::
Region(const TV& location,int& region_id,const T distance) const
{
    SEGMENT_3D<T> segment(X.x,X.y);
    T d1=segment.Distance_From_Point_To_Segment(location);
    segment.X.x=X.z;
    T d2=segment.Distance_From_Point_To_Segment(location);
    segment.X.y=X.x;
    T d3=segment.Distance_From_Point_To_Segment(location);
    if(d1 > distance && d2 > distance && d3 > distance) return 2; // face is closest
    else if(d1 <= distance){
        if(d2 <= distance){region_id=2;return 0;} // vertex 2
        else if(d3 <= distance){region_id=1;return 0;} // vertex 1
        else{region_id=1;return 1;}} // edge 1
    else if(d2 <= distance){
        if(d3 <= distance){region_id=3;return 0;} // vertex 3
        else{region_id=2;return 1;}} // edge 2
    else{region_id=3;return 1;} // edge 3
}
//#####################################################################
// Function Closest_Point
//#####################################################################
template<class T> VECTOR<T,3> TRIANGLE_3D<T>::
Closest_Point(const TV& location,TV& weights) const
{
    weights=Barycentric_Coordinates(location);
    // project closest point to the triangle if it's not already inside it
    if(weights.x<0){
        T a12=SEGMENT_3D<T>::Interpolation_Fraction(location,X.y,X.z); // Check edge X.y--X.z
        if(a12<0){
            if(weights.z<0){ // Closest point is on edge X.x--X.y
                T a01=clamp<T>(SEGMENT_3D<T>::Interpolation_Fraction(location,X.x,X.y),0,1);weights=TV(1-a01,a01,0);return weights.x*X.x+weights.y*X.y;}
            else{weights=TV(0,1,0);return X.y;}} // Closest point is X.y
        else if(a12>1){
            if(weights.y<0){ // Closest point is on edge X.x--X.z
                T a02=clamp<T>(SEGMENT_3D<T>::Interpolation_Fraction(location,X.x,X.z),0,1);weights=TV(1-a02,0,a02);return weights.x*X.x+weights.z*X.z;}
            else{weights=TV(0,0,1);return X.z;}} // Closest point is X.z
        else{weights=TV(0,1-a12,a12);return weights.y*X.y+weights.z*X.z;}} // Closest point is on edge X.y--X.z
    else if(weights.y<0){
        T a02=SEGMENT_3D<T>::Interpolation_Fraction(location,X.x,X.z); // Check edge X.x--X.z
        if(a02<0){
            if(weights.z<0){ // Closest point is on edge X.x--X.y
                T a01=clamp<T>(SEGMENT_3D<T>::Interpolation_Fraction(location,X.x,X.y),0,1);weights=TV(1-a01,a01,0);return weights.x*X.x+weights.y*X.y;}
            else{weights=TV(1,0,0);return X.x;}} // Closest point is X.x
        else if(a02>1){weights=TV(0,0,1);return X.z;} // Closest point is X.z
        else{weights=TV(1-a02,0,a02);return weights.x*X.x+weights.z*X.z;}} // Closest point is on edge X.x--X.z
    else if(weights.z<0){ // Closest point is on edge X.x--X.y
        T a01=clamp<T>(SEGMENT_3D<T>::Interpolation_Fraction(location,X.x,X.y),0,1);weights=TV(1-a01,a01,0);return weights.x*X.x+weights.y*X.y;}
    return weights.x*X.x+weights.y*X.y+weights.z*X.z; // Point is interior to the triangle
}
//#####################################################################
// Function Distance_To_Triangle
//#####################################################################
template<class T> T TRIANGLE_3D<T>::
Distance_To_Triangle(const TV& location) const
{   
    TV weights,projected_point=Closest_Point(location,weights);
    return (location-projected_point).Magnitude();
}
//#####################################################################
// Function Minimum_Angle
//#####################################################################
template<class T> T TRIANGLE_3D<T>::
Minimum_Angle() const
{
    TV s1=(X.x-X.y).Normalized(),s2=(X.y-X.z).Normalized(),s3=(X.z-X.x).Normalized();
    return acos(max(TV::Dot_Product(s1,-s2),TV::Dot_Product(-s1,s3),TV::Dot_Product(s2,-s3)));
}
//#####################################################################
// Function Maximum_Angle
//#####################################################################
template<class T> T TRIANGLE_3D<T>::
Maximum_Angle() const
{
    TV s1=(X.x-X.y).Normalized(),s2=(X.y-X.z).Normalized(),s3=(X.z-X.x).Normalized();
    return acos(min(TV::Dot_Product(s1,-s2),TV::Dot_Product(-s1,s3),TV::Dot_Product(s2,-s3)));
}
//#####################################################################
// Function Signed_Solid_Angle
//#####################################################################
// positive for normals that point away from the center - not reliable if center is too close to the triangle face
template<class T> T TRIANGLE_3D<T>::
Signed_Solid_Angle(const TV& center) const
{
    TV r=(X.x-center).Normalized(),u=X.y-X.x,v=X.z-X.x;u-=TV::Dot_Product(u,r)*r;v-=TV::Dot_Product(v,r)*r;
    T solid_angle=-(T)pi+TV::Angle_Between(u,v);
    r=(X.y-center).Normalized();u=X.x-X.y,v=X.z-X.y;u-=TV::Dot_Product(u,r)*r;v-=TV::Dot_Product(v,r)*r;
    solid_angle+=TV::Angle_Between(u,v);
    r=(X.z-center).Normalized();u=X.x-X.z,v=X.y-X.z;u-=TV::Dot_Product(u,r)*r;v-=TV::Dot_Product(v,r)*r;
    solid_angle+=TV::Angle_Between(u,v);
    solid_angle=max(T(0),min((T)(2*pi),solid_angle));
    if(TV::Dot_Product(r,Normal()) < 0) solid_angle*=(-1);
    return solid_angle;
}
//#####################################################################
// Function Point_Face_Interaction
//#####################################################################
// outputs unsigned distance
template<class T> bool TRIANGLE_3D<T>::
Point_Face_Interaction(const TV& x,const T interaction_distance,const bool allow_negative_weights,T& distance) const                       
{
    distance=Signed_Distance(x);
    return abs(distance)<=interaction_distance && Planar_Point_Inside_Triangle(x,allow_negative_weights?interaction_distance:0);
}
//#####################################################################
// Function Point_Face_Interaction_Data
//#####################################################################
template<class T> void TRIANGLE_3D<T>::
Point_Face_Interaction_Data(const TV& x,T& distance,TV& interaction_normal,VECTOR<T,TV::m+1>& weights,const bool perform_attractions) const                       
{
    interaction_normal=Normal();
    weights=Barycentric_Coordinates(x).Insert(-1,0);
    if(!perform_attractions && distance<0){distance*=-1;interaction_normal*=-1;} // distance > 0, interaction_normal points from the triangle to the point
}
//#####################################################################
// Function Point_Face_Interaction
//#####################################################################
template<class T> bool TRIANGLE_3D<T>::
Point_Face_Interaction(const TV& x,const TV& v,const TV& v1,const TV& v2,const TV& v3,const T interaction_distance,T& distance,
    TV& interaction_normal,VECTOR<T,TV::m+1>& weights,const bool allow_negative_weights,const bool exit_early) const
{      
    if(!Point_Face_Interaction(x,interaction_distance,allow_negative_weights,distance)) return false;
    if(!exit_early) Point_Face_Interaction_Data(x,distance,interaction_normal,weights,false); // relative speed is in the normal direction
    return true;
}
//#####################################################################
// Function Robust_Point_Triangle_Collision
//#####################################################################
template<class T> POINT_SIMPLEX_COLLISION_TYPE TRIANGLE_3D<T>::
Robust_Point_Face_Collision(const TRIANGLE_3D<T>& initial_triangle,const TRIANGLE_3D<T>& final_triangle,const TV& x,const TV& final_x,const T dt,
    const T collision_thickness,T& collision_time,TV& normal,VECTOR<T,TV::m+1>& weights)
{
    if(final_triangle.Point_Inside_Triangle(final_x,collision_thickness)){
        collision_time=dt;
        weights=final_triangle.Barycentric_Coordinates(final_x).Insert(-1,0);
        normal=final_triangle.Normal();
        return POINT_SIMPLEX_COLLISION_ENDS_INSIDE;}
    if(initial_triangle.Point_Inside_Triangle(x,collision_thickness)){
        collision_time=0;
        weights=initial_triangle.Barycentric_Coordinates(x).Insert(-1,0);
        normal=initial_triangle.Normal();
        return POINT_SIMPLEX_COLLISION_ENDS_OUTSIDE;}
    TV v1=(final_triangle.X.x-initial_triangle.X.x)/dt,v2=(final_triangle.X.y-initial_triangle.X.y)/dt,v3=(final_triangle.X.z-initial_triangle.X.z)/dt,v=(final_x-x)/dt;
    if(initial_triangle.Point_Face_Collision(x,v,v1,v2,v3,dt,collision_thickness,collision_time,normal,weights,false))
        return POINT_SIMPLEX_COLLISION_ENDS_OUTSIDE;
    return POINT_SIMPLEX_NO_COLLISION;
}
//#####################################################################
// Function Point_Triangle_Collision
//#####################################################################
template<class T> bool TRIANGLE_3D<T>::
Point_Face_Collision(const TV& x,const TV& v,const TV& v1,const TV& v2,const TV& v3,const T dt,const T collision_thickness,
    T& collision_time,TV& normal,VECTOR<T,TV::m+1>& weights,const bool exit_early) const
{
    // find cubic and compute the roots as possible collision times
    TV ABo=X.y-X.x,ABv=dt*(v2-v1),ACo=X.z-X.x,ACv=dt*(v3-v1);
    TV No=TV::Cross_Product(ABo,ACo),Nv=TV::Cross_Product(ABo,ACv)+TV::Cross_Product(ABv,ACo),Na=TV::Cross_Product(ABv,ACv);
    TV APo=x-X.x,APv=dt*(v-v1);

    CUBIC<double> cubic((double)TV::Dot_Product(Na,APv),(double)TV::Dot_Product(Nv,APv)+TV::Dot_Product(Na,APo),
                                       (double)TV::Dot_Product(No,APv)+TV::Dot_Product(Nv,APo),(double)TV::Dot_Product(No,APo));
    double xmin=0,xmax=1.000001;
    int num_intervals=0;VECTOR<INTERVAL<double>,3> intervals;
    cubic.Compute_Intervals(xmin,xmax,num_intervals,intervals(0),intervals(1),intervals(2));
    if(!num_intervals) return false;
  
    // find and check roots
    T distance;
    ITERATIVE_SOLVER<double> iterative_solver;
    iterative_solver.tolerance=1e-14;
    for(int k=0;k<num_intervals;k++){
        collision_time=dt*(T)iterative_solver.Bisection_Secant_Root(cubic,intervals(k).min_corner,intervals(k).max_corner);
        TRIANGLE_3D<T> triangle(X.x+collision_time*v1,X.y+collision_time*v2,X.z+collision_time*v3);
        if(triangle.Point_Face_Interaction(x+collision_time*v,v,v1,v2,v3,collision_thickness,distance,normal,weights,true,exit_early)) return true;}

    return false;
}
//#####################################################################
// Function Clip_To_Box
//#####################################################################
template<class T> void TRIANGLE_3D<T>::
Clip_To_Box(const RANGE<TV>& box,ARRAY<TRIANGLE_3D<T> >& clipped_simplices) const
{
    // cut with all sides of box
    clipped_simplices.Remove_All();
    clipped_simplices.Append(*this);
    for(int axis=0;axis<TV::dimension;axis++){
        for(int i=clipped_simplices.m-1;i>=0;i--){
            Cut_With_Hyperplane_And_Discard_Outside_Simplices(clipped_simplices(i),PLANE<T>(-TV::Axis_Vector(axis),box.min_corner),clipped_simplices);
            clipped_simplices.Remove_Index_Lazy(i);}
        for(int i=clipped_simplices.m-1;i>=0;i--){
            // TODO: make this more efficient by not removing a triangle that is fully inside
            Cut_With_Hyperplane_And_Discard_Outside_Simplices(clipped_simplices(i),PLANE<T>(TV::Axis_Vector(axis),box.max_corner),clipped_simplices);
            clipped_simplices.Remove_Index_Lazy(i);}}
}
//#####################################################################
// Function Cut_With_Hyperplane_And_Discard_Outside_Simplices
//#####################################################################
template<class T> void TRIANGLE_3D<T>::
Cut_With_Hyperplane_And_Discard_Outside_Simplices(const TRIANGLE_3D<T>& triangle,const PLANE<T>& cutting_plane,ARRAY<TRIANGLE_3D<T> >& negative_triangles)
{
    TV phi_nodes;
    VECTOR<TV,3> X_nodes;
    X_nodes(0)=triangle.X.x;
    X_nodes(1)=triangle.X.y;
    X_nodes(2)=triangle.X.z;
    for(int i=0;i<3;i++){phi_nodes(i)=cutting_plane.Signed_Distance(X_nodes(i));}

    // left simplices are in the negative halfspace, right simplices in the positive halfspace
    int positive_count=0,single_node_sign;
    for(int i=0;i<3;i++) if(phi_nodes(i)>0) positive_count++;
    switch(positive_count){
        case 0: // in negative halfspace
            negative_triangles.Append(triangle);break;
        case 1:
        case 2:
            single_node_sign=positive_count==1?1:-1;
            // draw positive triangle. has correct positive/negative area based on whether triangle is backwards or not
            for(int i=0;i<3;i++)if(LEVELSET_UTILITIES<T>::Sign(phi_nodes(i))==single_node_sign){
                VECTOR<TV,2> interface_locations;int index=(i+1)%3;
                VECTOR<int,2> other_locations;
                for(int j=0;j<2;j++,index=(index+1)%3){
                    other_locations(j)=index;
                    interface_locations(j)=LINEAR_INTERPOLATION<T,TV>::Linear(X_nodes(i),X_nodes(index),LEVELSET_UTILITIES<T>::Theta(phi_nodes(i),phi_nodes(index)));}
                if(positive_count==1){ // add two triangles to negative triangles
                    negative_triangles.Append(TRIANGLE_3D<T>(interface_locations(0),X_nodes(other_locations(0)),X_nodes(other_locations(1))));
                    negative_triangles.Append(TRIANGLE_3D<T>(interface_locations(0),X_nodes(other_locations(1)),interface_locations(1)));}
                else // add triangle to negative_triangles
                    negative_triangles.Append(TRIANGLE_3D<T>(X_nodes(i),interface_locations(0),interface_locations(1)));
                return;}
        case 3: // in positive halfspace
            break;}
}
//#####################################################################
// Function Cut_With_Hyperplane
//#####################################################################
template<class T> void TRIANGLE_3D<T>::
Cut_With_Hyperplane(const TRIANGLE_3D<T>& triangle,const PLANE<T>& cutting_plane,ARRAY<TRIANGLE_3D<T> >& negative_triangles,ARRAY<TRIANGLE_3D<T> >& positive_triangles,T tol)
{
    TV phi_nodes;
    VECTOR<TV,3> X_nodes;
    X_nodes(0)=triangle.X.x;
    X_nodes(1)=triangle.X.y;
    X_nodes(2)=triangle.X.z;
    for(int i=0;i<3;i++) phi_nodes(i)=cutting_plane.Signed_Distance(X_nodes(i));

    int positive_count=0,single_node_sign;
    for(int i=0;i<3;i++) if(phi_nodes(i)>0) positive_count++;
    switch(positive_count){
        case 0: // in negative halfspace
            negative_triangles.Append(triangle);break;
        case 1:
        case 2:
            single_node_sign=positive_count==1?1:-1;
            // draw positive triangle. has correct positive/negative area based on whether triangle is backwards or not
            for(int i=0;i<3;i++)if(LEVELSET_UTILITIES<T>::Sign(phi_nodes(i))==single_node_sign){
                VECTOR<TV,2> interface_locations;int index=(i+1)%3;
                VECTOR<T,2> theta;
                VECTOR<int,2> other_locations;
                for(int j=0;j<2;j++,index=(index+1)%3){
                    other_locations(j)=index;
                    theta(j)=LEVELSET_UTILITIES<T>::Theta(phi_nodes(i),phi_nodes(index));
                    interface_locations(j)=LINEAR_INTERPOLATION<T,TV>::Linear(X_nodes(i),X_nodes(index),theta(j));}
                if(positive_count==1){ // add 2 negative and 1 positive
                    if(1-theta(0)>tol) negative_triangles.Append(TRIANGLE_3D<T>(interface_locations(0),X_nodes(other_locations(0)),X_nodes(other_locations(1))));
                    if(1-theta(1)>tol) negative_triangles.Append(TRIANGLE_3D<T>(interface_locations(0),X_nodes(other_locations(1)),interface_locations(1)));
                    if(theta(0)>tol && theta(1)>tol) positive_triangles.Append(TRIANGLE_3D<T>(X_nodes(i),interface_locations(0),interface_locations(1)));}
                else{ // add 2 positive and 1 negative
                    if(1-theta(0)>tol) positive_triangles.Append(TRIANGLE_3D<T>(interface_locations(0),X_nodes(other_locations(0)),X_nodes(other_locations(1))));
                    if(1-theta(1)>tol) positive_triangles.Append(TRIANGLE_3D<T>(interface_locations(0),X_nodes(other_locations(1)),interface_locations(1)));
                    if(theta(0)>tol && theta(1)>tol) negative_triangles.Append(TRIANGLE_3D<T>(X_nodes(i),interface_locations(0),interface_locations(1)));}
                return;}
        case 3: // in positive halfspace
            positive_triangles.Append(triangle);break;
            break;}
}
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool TRIANGLE_3D<T>::
Intersects(const TRIANGLE_3D<T>& triangle,T theta_tol,INTERSECTS_HELPER* ih) const
{
    INTERSECTS_HELPER lh[2],*h=ih?ih:lh; // h[0].w is a z, h[1].w is a y.
    VECTOR<TV,3> Y[2]={X,triangle.X};
    for(int i=0;i<2;i++) h[i].n=(Y[i][0]-Y[i][2]).Cross(Y[i][1]-Y[i][2]);
    TV r=h[0].n.Cross(h[1].n);
    for(int k=0;k<2;k++){
        h[k].neg=0;
        h[k].pos=0;
        for(int i=0;i<3;i++){
            h[k].x(i)=r.Dot(Y[k][i]);
            h[k].w(i)=h[1-k].n.Dot(Y[k][i]-Y[1-k][0]);
            h[k].neg|=(h[k].w(i)<0)<<i;
            h[k].pos|=(h[k].w(i)>0)<<i;}}
    for(int k=0;k<2;k++) if(!h[k].neg || !h[k].pos) return false;
    for(int k=0;k<2;k++){
        for(h[k].is=0;h[k].is<3;h[k].is++) if(h[k].pos==(1<<h[k].is) || h[k].neg==(1<<h[k].is)) break;
        h[k].i[0]=h[k].is+1;
        if(h[k].i[0]>2) h[k].i[0]=0;
        h[k].i[1]=3-h[k].is-h[k].i[0];
        for(int j=0;j<2;j++){
            h[k].th[j]=h[k].w(h[k].is)/(h[k].w(h[k].is)-h[k].w(h[k].i[j]));
            h[k].t[j]=h[k].x(h[k].is)+h[k].th[j]*(h[k].x(h[k].i[j])-h[k].x(h[k].is));}}

    if(theta_tol>0){
        for(int k=0;k<2;k++)
            for(int j=0;j<2;j++)
                if(h[k].th[j]<theta_tol){
                    h[k].th[j]=0;
                    h[k].th[1-j]=0;
                    h[k].neg&=~(1<<h[k].is);
                    h[k].pos&=~(1<<h[k].is);}
        for(int k=0;k<2;k++)
            for(int j=0;j<2;j++)
                if(1-h[k].th[j]<theta_tol){
                    h[k].th[j]=1;
                    h[k].neg&=~(1<<h[k].i[j]);
                    h[k].pos&=~(1<<h[k].i[j]);}
        for(int k=0;k<2;k++) if(!h[k].neg || !h[k].pos) return false;}

    for(int k=0;k<2;k++) if(min(h[k].t[0],h[k].t[1])>=max(h[1-k].t[0],h[1-k].t[1])) return false;
    return true;
}
//#####################################################################
// Function Barycentric_Coordinates
//#####################################################################
template<class T> VECTOR<T,3> TRIANGLE_3D<T>::
Barycentric_Coordinates(const TV& location,const TV& x0,const TV& x1,const TV& x2) // clockwise vertices
{
    TV u=x1-x0,v=x2-x0,w=location-x0;
    T u_dot_u=TV::Dot_Product(u,u),v_dot_v=TV::Dot_Product(v,v),u_dot_v=TV::Dot_Product(u,v),
        u_dot_w=TV::Dot_Product(u,w),v_dot_w=TV::Dot_Product(v,w);
    T denominator=u_dot_u*v_dot_v-sqr(u_dot_v),one_over_denominator;
    if(abs(denominator)>(T)1e-16) one_over_denominator=1/denominator;else one_over_denominator=(T)1e16;
    T a=(v_dot_v*u_dot_w-u_dot_v*v_dot_w)*one_over_denominator,b=(u_dot_u*v_dot_w-u_dot_v*u_dot_w)*one_over_denominator;
    return TV(1-a-b,a,b);
}
//#####################################################################
// Function Clamped_Barycentric_Coordinates
//#####################################################################
template<class T> VECTOR<T,3> TRIANGLE_3D<T>::
Clamped_Barycentric_Coordinates(const TV& location,const TV& x0,const TV& x1,const TV& x2,const T tolerance) // clockwise vertices
{
    TV u=x1-x0,v=x2-x0,w=location-x0;
    T u_dot_u=TV::Dot_Product(u,u),v_dot_v=TV::Dot_Product(v,v),u_dot_v=TV::Dot_Product(u,v),
       u_dot_w=TV::Dot_Product(u,w),v_dot_w=TV::Dot_Product(v,w);
    if(abs(u_dot_u)<tolerance){
        if(abs(v_dot_v)<tolerance) return TV((T)one_third,(T)one_third,(T)one_third); // single point
        T c=clamp(v_dot_w/v_dot_v,(T)0,(T)1);T a_and_b=(T).5*(1-c);return TV(a_and_b,a_and_b,c);} // x0 and x1 are a single point
    else if(abs(v_dot_v)<tolerance){
        T b=clamp(u_dot_w/u_dot_u,(T)0,(T)1);T a_and_c=(T).5*(1-b);return TV(a_and_c,b,a_and_c);} // x0 and x2 are a single point
    else{
        T denominator=u_dot_u*v_dot_v-sqr(u_dot_v); 
        if(abs(denominator)<tolerance){
            if(u_dot_v>0){ // u and v point in the same direction
                if(u_dot_u>u_dot_v){T b=clamp(u_dot_w/u_dot_u,(T)0,(T)1);return TV(1-b,b,0);}
                else{T c=clamp(v_dot_w/v_dot_v,(T)0,(T)1);return TV(1-c,0,c);}}
            else if(u_dot_w>0){T b=clamp(u_dot_w/u_dot_u,(T)0,(T)1);return TV(1-b,b,0);} // u and v point in opposite directions, and w is on the u segment
            else{T c=clamp(v_dot_w/v_dot_v,(T)0,(T)1);return TV(1-c,0,c);}} // u and v point in opposite directions, and w is on the v segment
        T one_over_denominator=1/denominator;
        T a=clamp((v_dot_v*u_dot_w-u_dot_v*v_dot_w)*one_over_denominator,(T)0,(T)1),b=clamp((u_dot_u*v_dot_w-u_dot_v*u_dot_w)*one_over_denominator,(T)0,(T)1);
        return TV(1-a-b,a,b);}
}
//#####################################################################
namespace PhysBAM{
template class TRIANGLE_3D<float>;
template class TRIANGLE_3D<double>;
}

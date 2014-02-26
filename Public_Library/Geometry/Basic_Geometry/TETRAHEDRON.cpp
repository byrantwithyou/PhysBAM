//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Neil Molino, Avi Robinson-Mosher, Andrew Selle, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TETRAHEDRON
//##################################################################### 
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Arrays/ARRAY_VIEW.h>
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <Tools/Vectors/VECTOR_3D.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> TETRAHEDRON<T>::
TETRAHEDRON()
    :X(TV(),TV(0,1,0),TV(1,0,0),TV(0,0,-1))
{
    Create_Triangles();
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> TETRAHEDRON<T>::
TETRAHEDRON(const TV& x1_input,const TV& x2_input,const TV& x3_input,const TV& x4_input)
    :X(x1_input,x2_input,x3_input,x4_input)
{
    Create_Triangles();
}
//#####################################################################
// Function Create_Triangles
//#####################################################################
template<class T> void TETRAHEDRON<T>::
Create_Triangles()
{
    int t=(TV::Dot_Product(TV::Cross_Product(X[1]-X[0],X[2]-X[0]),X[3]-X[0])<=0)+1,s=3-t;
    triangle(0).X(0)=X(0);
    triangle(0).X(s)=X(1);
    triangle(0).X(t)=X(2);
    triangle(1).X(0)=X(0);
    triangle(1).X(s)=X(3);
    triangle(1).X(t)=X(1);
    triangle(2).X(0)=X(0);
    triangle(2).X(s)=X(2);
    triangle(2).X(t)=X(3);
    triangle(3).X(0)=X(1);
    triangle(3).X(s)=X(3);
    triangle(3).X(t)=X(2);
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> TETRAHEDRON<T>::
Normal(const VECTOR<T,3>& location,const int aggregate) const
{
    return triangle(aggregate).Normal();
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool TETRAHEDRON<T>::
Inside(const VECTOR<T,3>& location,const T thickness_over_two) const
{
    for(int i=0;i<4;i++)
        if(!triangle(i).Inside_Plane(location,thickness_over_two))
            return false;
    return true;
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool TETRAHEDRON<T>::
Outside(const VECTOR<T,3>& location,const T thickness_over_two) const
{
    for(int i=0;i<4;i++)
        if(triangle(i).Outside_Plane(location,thickness_over_two))
            return true;
    return false;
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class T> bool TETRAHEDRON<T>::
Boundary(const VECTOR<T,3>& location,const T thickness_over_two) const
{
    return !Inside(location,thickness_over_two) && !Outside(location,thickness_over_two);
}
//#####################################################################
// Function Thickened
//#####################################################################
template<class T> TETRAHEDRON<T> TETRAHEDRON<T>::
Thickened(const T thickness_over_two) const
{
    assert(Signed_Volume()>0);
    T m1=thickness_over_two/triangle(3).Signed_Distance(X[0]),m2=thickness_over_two/triangle(2).Signed_Distance(X[1]),m3=thickness_over_two/triangle(1).Signed_Distance(X[2]),m4=thickness_over_two/triangle(0).Signed_Distance(X[3]);
    return TETRAHEDRON<T>(Point_From_Barycentric_Coordinates(VECTOR<T,3>(m1,m2,m3)),Point_From_Barycentric_Coordinates(VECTOR<T,3>((T)1-m2-m3-sqrt((T)3)*m4,m2,m3)),
                          Point_From_Barycentric_Coordinates(VECTOR<T,3>(m1,(T)1-m1-m3-sqrt((T)3)*m4,m3)),Point_From_Barycentric_Coordinates(VECTOR<T,3>(m1,m2,(T)1-m1-m2-sqrt((T)3)*m4)));
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > TETRAHEDRON<T>::
Bounding_Box() const
{      
    return RANGE<TV>::Bounding_Box(X);
}
//#####################################################################
// Function Surface
//#####################################################################
template<class T> VECTOR<T,3> TETRAHEDRON<T>::
Surface(const VECTOR<T,3>& location) const
{      
    if(Inside(location)){
        int t=0;T distance=FLT_MAX;
        for(int i=0;i<4;i++){
            T d=triangle(i).Normal().Dot(triangle(i).X.x-location);
            if(d<distance){t=i;distance=d;}}
        return location+distance*triangle(t).Normal();}
    else{
        VECTOR<T,3> surface_point(location);
        for(int i=0;i<4;i++)
            if(triangle(i).Lazy_Outside_Plane(surface_point))
                surface_point-=(location-triangle(i).X.x).Projected_On_Unit_Direction(triangle(i).Normal());
        return surface_point;}
}
//#####################################################################
// Function Closest_Point
//#####################################################################
template<class T> VECTOR<T,3> TETRAHEDRON<T>::
Closest_Point(const VECTOR<T,3>& location,VECTOR<T,3>& weights) const
{      
    if(Inside(location)){weights=First_Three_Barycentric_Coordinates(location);return location;}
    
    VECTOR<T,3> triangle_weights,triangle_closest_point=TRIANGLE_3D<T>(X[0],X[1],X[2]).Closest_Point(location,triangle_weights),closest_point=triangle_closest_point;
    weights=VECTOR<T,3>(triangle_weights.x,triangle_weights.y,triangle_weights.z);
    T triangle_distance_squared=(closest_point-location).Magnitude_Squared(),distance_squared=triangle_distance_squared;
    
    triangle_closest_point=TRIANGLE_3D<T>(X[0],X[1],X[3]).Closest_Point(location,triangle_weights);triangle_distance_squared=(triangle_closest_point-location).Magnitude_Squared();
    if(triangle_distance_squared<distance_squared){
        weights=VECTOR<T,3>(triangle_weights.x,triangle_weights.y,0);closest_point=triangle_closest_point;distance_squared=triangle_distance_squared;}
    
    triangle_closest_point=TRIANGLE_3D<T>(X[0],X[2],X[3]).Closest_Point(location,triangle_weights);triangle_distance_squared=(triangle_closest_point-location).Magnitude_Squared();
    if(triangle_distance_squared<distance_squared){
        weights=VECTOR<T,3>(triangle_weights.x,0,triangle_weights.y);closest_point=triangle_closest_point;distance_squared=triangle_distance_squared;}
    
    triangle_closest_point=TRIANGLE_3D<T>(X[1],X[2],X[3]).Closest_Point(location,triangle_weights);triangle_distance_squared=(triangle_closest_point-location).Magnitude_Squared();
    if(triangle_distance_squared<distance_squared){
        weights=VECTOR<T,3>(0,triangle_weights.x,triangle_weights.y);return triangle_closest_point;}
        
    return closest_point;
}
//#####################################################################
// Function Volume
//#####################################################################
template<class T> T TETRAHEDRON<T>::
Volume() const
{
    return ((T)1/6)*abs(VECTOR<T,3>::Dot_Product(VECTOR<T,3>::Cross_Product(X[1]-X[0],X[2]-X[0]),X[3]-X[0]));
}
//#####################################################################
// Function Signed_Volume
//#####################################################################
template<class T> T TETRAHEDRON<T>::
Signed_Volume() const
{  
    return ((T)1/6)*VECTOR<T,3>::Dot_Product(VECTOR<T,3>::Cross_Product(X[1]-X[0],X[2]-X[0]),X[3]-X[0]); 
}
//#####################################################################
// Function Minimum_Angle
//#####################################################################
template<class T> T TETRAHEDRON<T>::
Minimum_Angle() const
{  
    return min(triangle(0).Minimum_Angle(),triangle(1).Minimum_Angle(),triangle(2).Minimum_Angle(),triangle(3).Minimum_Angle());
}
//#####################################################################
// Function Maximum_Angle
//#####################################################################
template<class T> T TETRAHEDRON<T>::
Maximum_Angle() const
{  
    return max(triangle(0).Maximum_Angle(),triangle(1).Maximum_Angle(),triangle(2).Maximum_Angle(),triangle(3).Maximum_Angle());
}
//#####################################################################
// Function Minimum_Altitude
//#####################################################################
template<class T> T TETRAHEDRON<T>::
Minimum_Altitude() const
{  
    return min(abs(triangle(0).Signed_Distance(X[3])),abs(triangle(1).Signed_Distance(X[2])),abs(triangle(2).Signed_Distance(X[1])),abs(triangle(3).Signed_Distance(X[0])));
}
//#####################################################################
// Function Aspect_Ratio
//#####################################################################
template<class T> T TETRAHEDRON<T>::
Aspect_Ratio() const
{  
    return Maximum_Edge_Length()/Minimum_Altitude();
}
//#####################################################################
// Function Minimum_Dihedral_Angle
//#####################################################################
template<class T> T TETRAHEDRON<T>::
Minimum_Dihedral_Angle() const
{  
    return acos(-min(VECTOR<T,3>::Dot_Product(triangle(0).Normal(),triangle(1).Normal()),VECTOR<T,3>::Dot_Product(triangle(0).Normal(),triangle(2).Normal()),
                               VECTOR<T,3>::Dot_Product(triangle(0).Normal(),triangle(3).Normal()),VECTOR<T,3>::Dot_Product(triangle(1).Normal(),triangle(2).Normal()),
                               VECTOR<T,3>::Dot_Product(triangle(1).Normal(),triangle(3).Normal()),VECTOR<T,3>::Dot_Product(triangle(2).Normal(),triangle(3).Normal()))); 
}
//#####################################################################
// Function Maximum_Dihedral_Angle
//#####################################################################
template<class T> T TETRAHEDRON<T>::
Maximum_Dihedral_Angle() const
{  
    return acos(-max(VECTOR<T,3>::Dot_Product(triangle(0).Normal(),triangle(1).Normal()),VECTOR<T,3>::Dot_Product(triangle(0).Normal(),triangle(2).Normal()),
                               VECTOR<T,3>::Dot_Product(triangle(0).Normal(),triangle(3).Normal()),VECTOR<T,3>::Dot_Product(triangle(1).Normal(),triangle(2).Normal()),
                               VECTOR<T,3>::Dot_Product(triangle(1).Normal(),triangle(3).Normal()),VECTOR<T,3>::Dot_Product(triangle(2).Normal(),triangle(3).Normal())));
}
//#####################################################################
// Function Signed_Reciprocal_Aspect_Ratio
//#####################################################################
template<class T> T TETRAHEDRON<T>::
Signed_Reciprocal_Aspect_Ratio(const VECTOR<T,3>& x0,const VECTOR<T,3>& x1,const VECTOR<T,3>& x2,const VECTOR<T,3>& x3)
{
    return min(VECTOR<T,3>::Dot_Product(x3-x0,PLANE<T>::Normal(x0,x1,x2)),VECTOR<T,3>::Dot_Product(x2-x0,PLANE<T>::Normal(x0,x3,x1)),
                     VECTOR<T,3>::Dot_Product(x1-x0,PLANE<T>::Normal(x0,x2,x3)),VECTOR<T,3>::Dot_Product(x0-x1,PLANE<T>::Normal(x1,x3,x2)))/
              max((x0-x1).Magnitude(),(x0-x2).Magnitude(),(x0-x3).Magnitude(),(x1-x2).Magnitude(),(x1-x3).Magnitude(),(x2-x3).Magnitude());
}
//#####################################################################
// Function Negative_Material
//#####################################################################
template<class T> T TETRAHEDRON<T>::
Negative_Material(const ARRAY<VECTOR<T,3> >& X,const ARRAY<T>& phis,const VECTOR<int,4>& indices)
{
    VECTOR<T,4> local_phi;
    int positive_count=0;for(int i=0;i<4;i++) positive_count+=((local_phi[i]=phis(indices[i]))>0);
    switch(positive_count){
      case 0: return Signed_Volume(X(indices[0]),X(indices[1]),X(indices[2]),X(indices[3]));
      case 1:
        for(int i=0;i<4;i++)if(local_phi[i]>0){
            VECTOR<VECTOR<T,3>,3> interface_locations;int index=(i+1)%4;
            for(int j=0;j<3;j++,index=(index+1)%4)
                interface_locations[j]=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(X(indices[i]),X(indices[index]),LEVELSET_UTILITIES<T>::Theta(local_phi[i],local_phi[index]));
            if(i%2 == 0) exchange(interface_locations[0],interface_locations[2]);
            return Signed_Volume(X(indices[0]),X(indices[1]),X(indices[2]),X(indices[3]))+TETRAHEDRON<T>::Signed_Volume(X(indices[i]),interface_locations[0],interface_locations[1],interface_locations[2]);}
      case 2:{
        positive_count=0;int negative_count=0;int negative_indices[2],positive_indices[2];
        for(int i=0;i<4;i++){if(local_phi[i]<=0) negative_indices[negative_count++]=i;else positive_indices[positive_count++]=i;}
        if((negative_indices[1]-negative_indices[0])%2 == 1) exchange(positive_indices[0],positive_indices[1]);  // odd wrong, even right (odd=swap)
        VECTOR<T,3> interface_locations[2][2];
        for(int j=0;j<2;j++)for(int k=0;k<2;k++){
            int n=negative_indices[j],p=positive_indices[k];
            interface_locations[j][k]=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(X(indices[n]),X(indices[p]),LEVELSET_UTILITIES<T>::Theta(local_phi[n],local_phi[p]));}
        // the tets are: (1-,1+2-,2+2-,1-2+) and (1-,1-1+,1+2-,1-2+), and (1-,1+2-,2-,2-2+)
        VECTOR<T,3> n0=X(indices[negative_indices[0]]);
        return TETRAHEDRON<T>::Signed_Volume(n0,interface_locations[1][0],X(indices[negative_indices[1]]),interface_locations[1][1])+
            TETRAHEDRON<T>::Signed_Volume(n0,interface_locations[1][0],interface_locations[1][1],interface_locations[0][1])+
            TETRAHEDRON<T>::Signed_Volume(n0,interface_locations[0][0],interface_locations[1][0],interface_locations[0][1]);}
      case 3:
        for(int i=0;i<4;i++)if(local_phi[i]<=0){
            VECTOR<VECTOR<T,3>,3> interface_locations;int index=(i+1)%4;
            for(int j=0;j<3;j++,index=(index+1)%4)
                interface_locations[j]=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(X(indices[i]),X(indices[index]),LEVELSET_UTILITIES<T>::Theta(local_phi[i],local_phi[index]));
            if(i%2 == 0) exchange(interface_locations[0],interface_locations[2]);
            return -TETRAHEDRON<T>::Signed_Volume(X(indices[i]),interface_locations[0],interface_locations[2],interface_locations[3]);}
      default:return 0;}

    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Cut_With_Hyperplane
//#####################################################################
template<class T> void TETRAHEDRON<T>::
Cut_With_Hyperplane(ARRAY<VECTOR<T,3> >& X,const PLANE<T>& cutting_plane,const VECTOR<int,4>& indices,ARRAY<VECTOR<int,4> >& left_simplices,ARRAY<VECTOR<int,4> >& right_simplices)
{
    VECTOR<T,4> phis;
    VECTOR<VECTOR<T,3>,4> X_nodes;
    for(int i=0;i<4;i++){X_nodes[i]=X(indices[i]);phis[i]=cutting_plane.Signed_Distance(X_nodes[i]);}
    Cut_Simplex(X,indices,X_nodes,phis,left_simplices,right_simplices);
}
//#####################################################################
// Function Clip_To_Box
//#####################################################################
template<class T> void TETRAHEDRON<T>::
Clip_To_Box(const RANGE<TV>& box,ARRAY<TETRAHEDRON<T> >& clipped_simplices) const
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
template<class T> void TETRAHEDRON<T>::
Cut_With_Hyperplane_And_Discard_Outside_Simplices(const TETRAHEDRON<T>& tetrahedron,const PLANE<T>& cutting_plane,ARRAY<TETRAHEDRON<T> >& negative_tetrahedra) const
{
    // left simplices are in the negative halfspace, right simplices in the positive halfspace
    VECTOR<T,4> phi_nodes;
    VECTOR<VECTOR<T,3>,4> X_nodes;X_nodes(0)=tetrahedron.X[0];X_nodes(1)=tetrahedron.X[1];X_nodes(2)=tetrahedron.X[2];X_nodes(3)=tetrahedron.X[3];
    for(int i=0;i<4;i++){phi_nodes(i)=cutting_plane.Signed_Distance(X_nodes(i));}

    int positive_count=0;
    for(int i=0;i<4;i++) if(phi_nodes[i]>0) positive_count++;
    switch(positive_count){
        case 0: // we are in the negative halfspace
            negative_tetrahedra.Append(tetrahedron);
            break;
        case 1: // tet in positive halfspace, three tets in negative
            for(int i=0;i<4;i++)if(phi_nodes[i]>0){
                VECTOR<VECTOR<T,3>,3> interface_locations;int index=(i+1)%4;
                VECTOR<int,3> other_indices;
                for(int j=0;j<3;j++,index=(index+1)%4){
                    other_indices[j]=index;
                    interface_locations[j]=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(X_nodes[i],X_nodes[index],LEVELSET_UTILITIES<T>::Theta(phi_nodes[i],phi_nodes[index]));}
                if(i%2 == 0){exchange(interface_locations[0],interface_locations[2]);exchange(other_indices[1],other_indices[3]);}
                // (i1,o1,o2,o3), (i2,i1,o2,o3), (i3,i1,i2,o3)
                negative_tetrahedra.Append(TETRAHEDRON<T>(interface_locations[0],X_nodes(other_indices[0]),X_nodes(other_indices[1]),X_nodes(other_indices[2])));
                negative_tetrahedra.Append(TETRAHEDRON<T>(interface_locations[1],interface_locations[0],X_nodes(other_indices[1]),X_nodes(other_indices[2])));
                negative_tetrahedra.Append(TETRAHEDRON<T>(interface_locations[2],interface_locations[0],interface_locations[1],X_nodes(other_indices[2])));
                return;}
        case 2:{ // three tets in each halfspace
            positive_count=0;int negative_count=0;int negative_indices[2],positive_indices[2];
            for(int i=0;i<4;i++){if(phi_nodes[i]<=0) negative_indices[negative_count++]=i;else positive_indices[positive_count++]=i;}
            if((negative_indices[1]-negative_indices[0])%2 == 0) exchange(positive_indices[0],positive_indices[1]);  // odd wrong, even right (odd=swap)
            VECTOR<T,3> interface_locations[2][2];
            for(int j=0;j<2;j++)for(int k=0;k<2;k++){
                int n=negative_indices[j],p=positive_indices[k];
                interface_locations[j][k]=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(X_nodes[n],X_nodes[p],LEVELSET_UTILITIES<T>::Theta(phi_nodes[n],phi_nodes[p]));}
            // the tets are: (1-,1+2-,2+2-,1-2+) and (1-,1-1+,1+2-,1-2+), and (1-,1+2-,2-,2-2+)
            VECTOR<T,3> n0=X_nodes(negative_indices[0]);
            negative_tetrahedra.Append(TETRAHEDRON<T>(n0,interface_locations[1][0],X_nodes(negative_indices[1]),interface_locations[1][1]));
            negative_tetrahedra.Append(TETRAHEDRON<T>(n0,interface_locations[1][0],interface_locations[1][1],interface_locations[0][1]));
            negative_tetrahedra.Append(TETRAHEDRON<T>(n0,interface_locations[0][0],interface_locations[1][0],interface_locations[0][1]));
            break;}
        case 3: // tet in negative halfspace, three tets in positive
            for(int i=0;i<4;i++)if(phi_nodes[i]<=0){
                VECTOR<VECTOR<T,3>,3> interface_locations;int index=(i+1)%4;
                VECTOR<int,3> other_indices;
                for(int j=0;j<3;j++,index=(index+1)%4){
                    other_indices[j]=index;
                    interface_locations[j]=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(X_nodes[i],X_nodes[index],LEVELSET_UTILITIES<T>::Theta(phi_nodes[i],phi_nodes[index]));}
                if(i%2 == 0){exchange(interface_locations[0],interface_locations[2]);exchange(other_indices[0],other_indices[2]);}
                negative_tetrahedra.Append(TETRAHEDRON<T>(X_nodes(i),interface_locations[0],interface_locations[1],interface_locations[2]));
                return;}
        case 4:
            break;}
}
//#####################################################################
// Function Cut_Simplex
//#####################################################################
template<class T> static inline void Add_Points_As_Tetrahedron(const ARRAY<VECTOR<T,3> >& X,ARRAY<VECTOR<int,4> >& tets,const int x0,const int x1,const int x2,const int x3)
{
    /*if(TETRAHEDRON<T>::Signed_Volume(X(x0),X(x1),X(x2),X(x3)))*/ tets.Append(VECTOR<int,4>(x0,x1,x2,x3));
}
template<class T> void TETRAHEDRON<T>::
Cut_Simplex(ARRAY<VECTOR<T,3> >& X,const VECTOR<int,4>& indices,const VECTOR<VECTOR<T,3>,4>& X_nodes,const VECTOR<T,4>& phi_nodes,ARRAY<VECTOR<int,4> >& left_simplices,
    ARRAY<VECTOR<int,4> >& right_simplices)
{
    // left simplices are in the negative halfspace, right simplices in the positive halfspace
    int positive_count=0;
    for(int i=0;i<4;i++) if(phi_nodes[i]>0) positive_count++;
    switch(positive_count){
      case 0: // we are in the negative halfspace
        Add_Points_As_Tetrahedron(X,left_simplices,indices[0],indices[1],indices[2],indices[3]);
        break;
      case 1: // tet in positive halfspace, three tets in negative
        for(int i=0;i<4;i++)if(phi_nodes[i]>0){
            VECTOR<int,3> interface_locations;int index=(i+1)%4;
            VECTOR<int,3> other_indices;
            for(int j=0;j<3;j++,index=(index+1)%4){
                other_indices[j]=indices[index];
                interface_locations[j]=X.Append(LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(X_nodes[i],X_nodes[index],LEVELSET_UTILITIES<T>::Theta(phi_nodes[i],phi_nodes[index])));}
            if(i%2 == 0){exchange(interface_locations[0],interface_locations[2]);exchange(other_indices[0],other_indices[2]);}
            Add_Points_As_Tetrahedron(X,right_simplices,indices[i],interface_locations[0],interface_locations[1],interface_locations[2]);
            // other three in left simplices
            // (i1,o1,o2,o3), (i2,i1,o2,o3), (i3,i1,i2,o3)
            Add_Points_As_Tetrahedron(X,left_simplices,interface_locations[0],other_indices[0],other_indices[1],other_indices[2]);
            Add_Points_As_Tetrahedron(X,left_simplices,interface_locations[1],interface_locations[0],other_indices[1],other_indices[2]);
            Add_Points_As_Tetrahedron(X,left_simplices,interface_locations[2],interface_locations[0],interface_locations[1],other_indices[2]);
            return;}
      case 2:{ // three tets in each halfspace
        positive_count=0;int negative_count=0;int negative_indices[2],positive_indices[2];
        for(int i=0;i<4;i++){if(phi_nodes[i]<=0) negative_indices[negative_count++]=i;else positive_indices[positive_count++]=i;}
        if((negative_indices[1]-negative_indices[0])%2 == 1) exchange(positive_indices[0],positive_indices[1]);  // odd wrong, even right (odd=swap)
        int interface_locations[2][2];
        for(int j=0;j<2;j++)for(int k=0;k<2;k++){
            int n=negative_indices[j],p=positive_indices[k];
            interface_locations[j][k]=X.Append(LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(X_nodes[n],X_nodes[p],LEVELSET_UTILITIES<T>::Theta(phi_nodes[n],phi_nodes[p])));}
        // the tets are: (1-,1+2-,2+2-,1-2+) and (1-,0+,1+2-,1-2+), and (1-,1+2-,2-,2-2+)
        // the tets are: (1+,1-2+,2-2+,1+2-) and (1+,1+1-,1-2+,1+2-), and (1+,1-2+,2+,2+2-)
        int n0=indices[negative_indices[0]];int p0=indices[positive_indices[0]];
        Add_Points_As_Tetrahedron(X,left_simplices,n0,interface_locations[1][0],indices[negative_indices[1]],interface_locations[1][1]);
        Add_Points_As_Tetrahedron(X,left_simplices,n0,interface_locations[1][0],interface_locations[1][1],interface_locations[0][1]);
        Add_Points_As_Tetrahedron(X,left_simplices,n0,interface_locations[0][0],interface_locations[1][0],interface_locations[0][1]);
        Add_Points_As_Tetrahedron(X,right_simplices,p0,interface_locations[0][1],indices[positive_indices[1]],interface_locations[1][1]);
        Add_Points_As_Tetrahedron(X,right_simplices,p0,interface_locations[0][1],interface_locations[1][1],interface_locations[1][0]);
        Add_Points_As_Tetrahedron(X,right_simplices,p0,interface_locations[0][0],interface_locations[0][1],interface_locations[1][0]);}
        break;
      case 3: // tet in negative halfspace, three tets in positive
        for(int i=0;i<4;i++)if(phi_nodes[i]<=0){
            VECTOR<int,3> interface_locations;int index=(i+1)%4;
            VECTOR<int,3> other_indices;
            for(int j=0;j<3;j++,index=(index+1)%4){
                other_indices[j]=indices[index];
                interface_locations[j]=X.Append(LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(X_nodes[i],X_nodes[index],LEVELSET_UTILITIES<T>::Theta(phi_nodes[i],phi_nodes[index])));}
            if(i%2 == 0){exchange(interface_locations[0],interface_locations[2]);exchange(other_indices[0],other_indices[2]);}
            Add_Points_As_Tetrahedron(X,left_simplices,indices[i],interface_locations[0],interface_locations[1],interface_locations[2]);
            Add_Points_As_Tetrahedron(X,right_simplices,interface_locations[0],other_indices[0],other_indices[1],other_indices[2]);
            Add_Points_As_Tetrahedron(X,right_simplices,interface_locations[1],interface_locations[0],other_indices[1],other_indices[2]);
            Add_Points_As_Tetrahedron(X,right_simplices,interface_locations[2],interface_locations[0],interface_locations[1],other_indices[2]);
            return;}
      case 4:
        Add_Points_As_Tetrahedron(X,right_simplices,indices[0],indices[1],indices[2],indices[3]);
        break;}
}
//#####################################################################
namespace PhysBAM{
template class TETRAHEDRON<float>;
template class TETRAHEDRON<double>;
}

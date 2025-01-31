//#####################################################################
// Copyright 2006, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_2D
//##################################################################### 
#include <Core/Arrays/ARRAY.h>
#include <Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <Geometry/Basic_Geometry/LINE_2D.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Function Negative_Material
//#####################################################################
template<class T> T TRIANGLE_2D<T>::
Negative_Material(const ARRAY<TV>& X,const ARRAY<T>& phis,const VECTOR<int,3>& indices)
{
    int positive_count=0;
    VECTOR<T,3> local_phi;
    T area=Signed_Area(X(indices[0]),X(indices[1]),X(indices[2])); // doesn't matter what's inside; also avoids duplicate node case (which messes up theta)
    // special case: if two phis are zero, that's an interface 
    if(!area) return 0;
    for(int i=0;i<3;i++) positive_count+=((local_phi[i]=phis(indices[i]))>0);
    switch(positive_count){
        case 0: return area;
        case 1:
            // draw positive triangle. has correct positive/negative area based on whether triangle is backwards or not
            for(int i=0;i<3;i++)if(local_phi[i]>0){
                VECTOR<TV,2> interface_locations;int index=(i+1)%3;
                for(int j=0;j<2;j++,index=(index+1)%3)
                    interface_locations[j]=LINEAR_INTERPOLATION<T,TV>::Linear(X(indices[i]),X(indices[index]),LEVELSET_UTILITIES<T>::Theta(local_phi[i],local_phi[index]));
                return area-TRIANGLE_2D<T>::Signed_Area(X(indices[i]),interface_locations[0],interface_locations[1]);}
        case 2:
            // draw negative triangle
            for(int i=0;i<3;i++)if(local_phi[i]<=0){
                VECTOR<TV,2> interface_locations;int index=(i+1)%3;
                for(int j=0;j<2;j++,index=(index+1)%3)
                    interface_locations[j]=LINEAR_INTERPOLATION<T,TV>::Linear(X(indices[i]),X(indices[index]),LEVELSET_UTILITIES<T>::Theta(local_phi[i],local_phi[index]));
                return TRIANGLE_2D<T>::Signed_Area(X(indices[i]),interface_locations[0],interface_locations[1]);}
        case 3: return (T)0;}

    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Cut_With_Hyperplane
//#####################################################################
template<class T> void TRIANGLE_2D<T>::
Cut_With_Hyperplane(ARRAY<TV>& X,const LINE_2D<T>& cutting_plane,const VECTOR<int,3>& indices,ARRAY<VECTOR<int,3> >& left_tris,ARRAY<VECTOR<int,3> >& right_tris)
{
    VECTOR<T,3> phis;
    VECTOR<TV,3> X_nodes;
    for(int i=0;i<3;i++){X_nodes[i]=X(indices[i]);phis[i]=cutting_plane.Signed_Distance(X_nodes[i]);}
    Cut_Simplex(X,indices,X_nodes,phis,left_tris,right_tris);
}
//#####################################################################
// Function Cut_Simplex
//#####################################################################
template<class T> static inline void Add_Points_As_Triangle(const ARRAY<VECTOR<T,2> >& X,ARRAY<VECTOR<int,3> >& tris,const int x0,const int x1,const int x2)
{
    tris.Append(VECTOR<int,3>(x0,x1,x2));
}
template<class T> void TRIANGLE_2D<T>::
Cut_Simplex(ARRAY<TV>& X,const VECTOR<int,3>& indices,const VECTOR<TV,3>& X_nodes,const VECTOR<T,3>& phi_nodes,ARRAY<VECTOR<int,3> >& left_tris,ARRAY<VECTOR<int,3> >& right_tris)
{
    int positive_count=0;
    for(int i=0;i<3;i++) if(phi_nodes[i]>0) positive_count++;
    switch(positive_count){
      case 0: // place in left list
        Add_Points_As_Triangle(X,left_tris,indices[0],indices[1],indices[2]);break;
      case 1:
        // draw positive triangle. has correct positive/negative area based on whether triangle is backwards or not
        for(int i=0;i<3;i++)if(phi_nodes[i]>0){
            VECTOR<int,2> interface_locations;int index=(i+1)%3;
            VECTOR<int,2> other_locations;
            for(int j=0;j<2;j++,index=(index+1)%3){
                other_locations[j]=indices[index];
                interface_locations[j]=X.Append(LINEAR_INTERPOLATION<T,TV>::Linear(X_nodes[i],X_nodes[index],LEVELSET_UTILITIES<T>::Theta(phi_nodes[i],phi_nodes[index])));}
            // add triangle to right tris
            Add_Points_As_Triangle(X,right_tris,indices[i],interface_locations[0],interface_locations[1]);
            // add two triangles to left tris
            Add_Points_As_Triangle(X,left_tris,interface_locations[0],other_locations[0],other_locations[1]);
            Add_Points_As_Triangle(X,left_tris,interface_locations[0],other_locations[1],interface_locations[1]);
            return;}
      case 2:
        // draw negative triangle
        for(int i=0;i<3;i++)if(phi_nodes[i]<=0){
            VECTOR<int,2> interface_locations;int index=(i+1)%3;
            VECTOR<int,2> other_locations;
            for(int j=0;j<2;j++,index=(index+1)%3){
                other_locations[j]=indices[index];
                interface_locations[j]=X.Append(LINEAR_INTERPOLATION<T,TV>::Linear(X_nodes[i],X_nodes[index],LEVELSET_UTILITIES<T>::Theta(phi_nodes[i],phi_nodes[index])));}
            // add triangle to left tris
            Add_Points_As_Triangle(X,left_tris,indices[i],interface_locations[0],interface_locations[1]);
            // add two triangles to right tris
            Add_Points_As_Triangle(X,right_tris,interface_locations[0],other_locations[0],other_locations[1]);
            Add_Points_As_Triangle(X,right_tris,interface_locations[0],other_locations[1],interface_locations[1]);
            return;}
      case 3: // place in right list
        Add_Points_As_Triangle(X,right_tris,indices[0],indices[1],indices[2]);break;}
}
//#####################################################################
// Function Clip_To_Box
//#####################################################################
template<class T> void TRIANGLE_2D<T>::
Clip_To_Box(const RANGE<TV>& box,ARRAY<TRIANGLE_2D<T> >& clipped_simplices) const
{
    // cut with all sides of box
    clipped_simplices.Remove_All();
    clipped_simplices.Append(*this);
    for(int axis=0;axis<TV::m;axis++){
        for(int i=clipped_simplices.m-1;i>=0;i--){
            Cut_With_Hyperplane_And_Discard_Outside_Simplices(clipped_simplices(i),LINE_2D<T>(-TV::Axis_Vector(axis),box.min_corner),clipped_simplices);
            clipped_simplices.Remove_Index_Lazy(i);}
        for(int i=clipped_simplices.m-1;i>=0;i--){
            // TODO: make this more efficient by not removing a triangle that is fully inside
            Cut_With_Hyperplane_And_Discard_Outside_Simplices(clipped_simplices(i),LINE_2D<T>(TV::Axis_Vector(axis),box.max_corner),clipped_simplices);
            clipped_simplices.Remove_Index_Lazy(i);}}
}
//#####################################################################
// Function Cut_With_Hyperplane_And_Discard_Outside_Simplices
//#####################################################################
template<class T> void TRIANGLE_2D<T>::
Cut_With_Hyperplane_And_Discard_Outside_Simplices(const TRIANGLE_2D<T>& triangle,const LINE_2D<T>& cutting_plane,ARRAY<TRIANGLE_2D<T> >& negative_triangles) const
{
    VECTOR<T,3> phi_nodes;
    VECTOR<VECTOR<T,2>,3> X_nodes;X_nodes(0)=triangle.X[0];X_nodes(1)=triangle.X[1];X_nodes(2)=triangle.X[2];
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
                VECTOR<VECTOR<T,2>,2> interface_locations;int index=(i+1)%3;
                VECTOR<int,2> other_locations;
                for(int j=0;j<2;j++,index=(index+1)%3){
                    other_locations(j)=index;
                    interface_locations(j)=LINEAR_INTERPOLATION<T,VECTOR<T,2> >::Linear(X_nodes(i),X_nodes(index),LEVELSET_UTILITIES<T>::Theta(phi_nodes(i),phi_nodes(index)));}
                if(positive_count==1){ // add two triangles to negative triangles
                    negative_triangles.Append(TRIANGLE_2D<T>(interface_locations(0),X_nodes(other_locations(0)),X_nodes(other_locations(1))));
                    negative_triangles.Append(TRIANGLE_2D<T>(interface_locations(0),X_nodes(other_locations(1)),interface_locations(1)));}
                else // add triangle to negative_triangles
                    negative_triangles.Append(TRIANGLE_2D<T>(X_nodes(i),interface_locations(0),interface_locations(1)));
                return;}
        case 3: // in positive halfspace
            break;}
}
//#####################################################################
// Function Cut_With_Hyperplane_And_Discard_Outside_Simplices
//#####################################################################
template<class T> bool TRIANGLE_2D<T>::
Intersects(const TRIANGLE_2D& tri) const
{
    const TRIANGLE_2D* tris[2]={this,&tri};
    int cs=0;

    for(int t=0;t<2;t++)
        for(int v=0;v<3;v++){
            for(int e=0;e<3;e++)
                cs|=(Signed_Area(tris[t]->X(v),tris[1-t]->X(e),tris[1-t]->X((e+1)%3))>=0)<<(9*t+3*v+e);
            int x=(cs>>(9*t+3*v))&7;
            if(x==7 || x==0) return true;} // Vertex inside

    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++){
            int in=(i+1)%3,jn=(j+1)%3;
            if(((cs>>(3*j+i))^(cs>>(3*jn+i)))&1)
                if(((cs>>(9+3*i+j))^(cs>>(9+3*in+j)))&1)
                    return true;} // Edge-edge intersection

    return false;
}
//#####################################################################
#define INSTANTIATION_HELPER(T) \
    template T TRIANGLE_2D<T>::Negative_Material(const ARRAY<VECTOR<T,2> >& X,const ARRAY<T>& phis,const VECTOR<int,3>& indices);\
    template void TRIANGLE_2D<T>::Cut_With_Hyperplane(ARRAY<VECTOR<T,2> >& X,const LINE_2D<T>& cutting_plane,const VECTOR<int,3>& indices,ARRAY<VECTOR<int,3> >& left_tris, \
        ARRAY<VECTOR<int,3> >& right_tris); \
    template void TRIANGLE_2D<T>::Cut_Simplex(ARRAY<VECTOR<T,2> >& X,const VECTOR<int,3>& indices,const VECTOR<VECTOR<T,2>,3>& X_nodes,const VECTOR<T,3>& phi_nodes, \
        ARRAY<VECTOR<int,3> >& left_tris,ARRAY<VECTOR<int,3> >& right_tris); \
    template void TRIANGLE_2D<T>::Clip_To_Box(const RANGE<VECTOR<T,2> >& box,ARRAY<TRIANGLE_2D<T> >& clipped_simplices) const; \
    template void TRIANGLE_2D<T>::Cut_With_Hyperplane_And_Discard_Outside_Simplices(const TRIANGLE_2D<T>& triangle,const LINE_2D<T>& cutting_plane,ARRAY<TRIANGLE_2D<T> >& negative_triangles) const;

INSTANTIATION_HELPER(float)
template bool TRIANGLE_2D<float>::Intersects(TRIANGLE_2D<float> const&) const;
INSTANTIATION_HELPER(double)
template bool TRIANGLE_2D<double>::Intersects(TRIANGLE_2D<double> const&) const;

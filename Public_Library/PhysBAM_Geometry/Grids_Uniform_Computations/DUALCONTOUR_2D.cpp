//#####################################################################
// Copyright 2004-2008, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DUALCONTOUR_2D
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Math_Tools/sign.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/DUALCONTOUR_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
using namespace PhysBAM;
//#####################################################################
// Function Generate_Topology
//#####################################################################
template<class T> void DUALCONTOUR_2D<T>::
Generate_Topology()
{
    int m=levelset.grid.counts.x,n=levelset.grid.counts.y;
    topology.Preallocate(5000);
    vertices.Resize(0,m-1,0,n-1);
    TV_INT i;
    for(i.x=0;i.x<m-1;i.x++)for(i.y=1;i.y<n-1;i.y++){ // generate one segment per x edge
        T phi1=levelset.phi(i),phi2=levelset.phi(i.x+1,i.y);
        int index=vertices.Standard_Index(i),index_j=index-1;
        if(phi1<=contour_value&&phi2>contour_value) topology.Append(VECTOR<int,2>(index_j,index));
        if(phi1>contour_value&&phi2<=contour_value) topology.Append(VECTOR<int,2>(index,index_j));}
    for(i.x=1;i.x<m-1;i.x++)for(i.y=0;i.y<n-1;i.y++){ // generate one segment per y edge
        T phi1=levelset.phi(i),phi2=levelset.phi(i.x,i.y+1);
        int index=vertices.Standard_Index(i),index_i=index-vertices.counts.y;
        if(phi1<=contour_value&&phi2>contour_value) topology.Append(VECTOR<int,2>(index_i,index));
        if(phi1>contour_value && phi2<=contour_value) topology.Append(VECTOR<int,2>(index,index_i));}
    topology.Compact();
    for(int t=0;t<topology.m;t++){
        int i,j;topology(t).Get(i,j);
        vertices.array(i)=vertices.array(j)=0;}
}
//#####################################################################
// Function Generate_Vertices
//#####################################################################
template<class T> void DUALCONTOUR_2D<T>::
Generate_Vertices()
{
    int m=grid.counts.x,n=grid.counts.y;
    int count=vertices.array.Count_Matches(0);
    geometry.Resize(count);
    normals.Resize(count);
    levelset.Compute_Normals();
    int vertex=0;
    for(int i=0;i<m-1;i++)for(int j=0;j<n-1;j++) if(vertices(i,j)>=0){ // generate vertices where needed
        vertices(i,j)=vertex++;
        TV position=grid.Center(TV_INT(i,j));
        TV position_guess;
        TV normal=levelset.Normal(position);
        T phi=levelset.Phi(position);
        T delta_phi;
        T dx=grid.dX.Min();
        bool failed_in_positive_normal=false;
        int iterations=0;
        if(is_distance_field) while(abs(phi-contour_value)>1e-5*grid.dX.Min() && (iterations++)<10){
            position-=(phi-contour_value)*normal;
            phi=levelset.Phi(position);
            normal=levelset.Normal(position);}
        else while(abs(phi-contour_value)>1e-5*grid.dX.Min() && (iterations++)<10){
            delta_phi=levelset.Phi(position+normal*dx)-phi;
            if(abs(delta_phi)>1e-8*grid.dX.Min()) position_guess=position-dx*((phi-contour_value)/delta_phi)*normal; // Newton-Raphson
            if(abs(levelset.Phi(position_guess)-contour_value)>abs(phi-contour_value)){
                if(!failed_in_positive_normal){failed_in_positive_normal=true;normal=-normal;} // Look away from a shock
                else break;}
            else{
                failed_in_positive_normal=false;
                position=position_guess;
                phi=levelset.Phi(position);
                normal=levelset.Normal(position);}}
        geometry(vertices(i,j))=position;normals(vertices(i,j))=normal;}
}
//#####################################################################
// Function Ensure_Vertices_In_Correct_Cells
//#####################################################################
template<class T> void DUALCONTOUR_2D<T>::
Ensure_Vertices_In_Correct_Cells()
{
    int vertex=0;
    TV_INT i;
    for(i.x=0;i.x<grid.counts.x-1;i.x++) for(i.y=0;i.y<grid.counts.y-1;i.y++) if(vertices(i)>=0){
        TV_INT v=grid.Cell(geometry(vertex),0);
        if(i!=v){
            TV cell_center=grid.Center(i);TV offset=(T).5*grid.dX;
            geometry(vertex)=RANGE<TV>(cell_center-offset,cell_center+offset).Surface(geometry(vertex));}
        vertex++;}
}
//#####################################################################
// Function Get_Segmented_Curve
//#####################################################################
template<class T> SEGMENTED_CURVE_2D<T>* DUALCONTOUR_2D<T>::
Get_Segmented_Curve()
{
    SEGMENTED_CURVE_2D<T>* curve=SEGMENTED_CURVE_2D<T>::Create();
    curve->particles.Add_Elements(geometry.m);
    curve->particles.X=geometry;
    curve->mesh.number_nodes=geometry.m;
    curve->mesh.elements.Exact_Resize(topology.m);
    for(int t=0;t<topology.m;t++){
        int i,j;topology(t).Get(i,j);i=vertices.array(i);j=vertices.array(j);
        curve->mesh.elements(t).Set(i,j);}
    curve->Update_Segment_List();
    return curve;
}
//#####################################################################
// Function Get_Triangulated_Surface
//#####################################################################
template<class T> TRIANGULATED_AREA<T>* DUALCONTOUR_2D<T>::
Get_Triangulated_Area(const int sign)
{
    TRIANGULATED_AREA<T>* area=TRIANGULATED_AREA<T>::Create();
    GEOMETRY_PARTICLES<TV>& particles=area->particles;
    particles.Preallocate(grid.counts.x*grid.counts.y*4);
    TRIANGLE_MESH& mesh=area->mesh;
    int triangle_count=0;for(int i=1;i<grid.counts.x-1;i++) for(int j=1;j<grid.counts.y-1;j++) if(levelset.phi(i,j)*sign>=0) triangle_count+=4;
    mesh.elements.Exact_Resize(triangle_count);
    int current_triangle=0;
    GRID<TV> mac_grid=grid.Get_MAC_Grid();
    for(int i=1;i<grid.counts.x-1;i++) for(int j=1;j<grid.counts.y-1;j++) if(levelset.phi(i,j)*sign>=0){
        int v0=particles.Add_Element(),v1=particles.Add_Element(),v2=particles.Add_Element(),v3=particles.Add_Element(),v4=particles.Add_Element();
        particles.X(v0)=grid.X(TV_INT(i,j));
        particles.X(v1)=vertices(i-1,j-1)>=0?geometry(vertices(i-1,j-1)):mac_grid.X(TV_INT(i-1,j-1));
        particles.X(v2)=vertices(i,j-1)>=0?geometry(vertices(i,j-1)):mac_grid.X(TV_INT(i,j-1));
        particles.X(v3)=vertices(i,j)>=0?geometry(vertices(i,j)):mac_grid.X(TV_INT(i,j));
        particles.X(v4)=vertices(i-1,j)>=0?geometry(vertices(i-1,j)):mac_grid.X(TV_INT(i-1,j));
        mesh.elements(current_triangle++).Set(v2,v1,v0);
        mesh.elements(current_triangle++).Set(v3,v2,v0);
        mesh.elements(current_triangle++).Set(v4,v3,v0);
        mesh.elements(current_triangle++).Set(v1,v4,v0);}
    mesh.number_nodes=particles.Size();
    return area;
}
//#####################################################################
template class DUALCONTOUR_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DUALCONTOUR_2D<double>;
#endif

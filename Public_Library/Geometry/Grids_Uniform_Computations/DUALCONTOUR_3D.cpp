//#####################################################################
// Copyright 2004-2012, Alexey, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DUALCONTOUR_3D
//#####################################################################
#include <Tools/Math_Tools/RANGE.h>
#include <Geometry/Grids_Uniform_Computations/DUALCONTOUR_3D.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Function Generate_Topology
//#####################################################################
template<class T> void DUALCONTOUR_3D<T>::
Generate_Topology()
{
    int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z;
    topology.Preallocate(5000);
    vertices.Resize(0,m-1,0,n-1,0,mn-1);
    int z=vertices.stride.y,yz=vertices.stride.x;
    TV_INT i;
    for(i.x=0;i.x<m-1;i.x++)for(i.y=1;i.y<n-1;i.y++)for(i.z=1;i.z<mn-1;i.z++){ // generate one triangle pair per x edge
        T phi1=levelset.phi(i),phi2=levelset.phi(i.x+1,i.y,i.z);
        int index=vertices.Standard_Index(i),index_j=index-z,index_ij=index-1,index_j_ij=index_j-1;
        if(phi1<=0&&phi2>0)topology.Append(VECTOR<int,4>(index_j_ij,index_ij,index,index_j));
        else if(phi1>0&&phi2<=0)topology.Append(VECTOR<int,4>(index_j,index,index_ij,index_j_ij));}
    for(i.x=1;i.x<m-1;i.x++)for(i.y=0;i.y<n-1;i.y++)for(i.z=1;i.z<mn-1;i.z++){ // generate one triangle pair per y edge
        T phi1=levelset.phi(i),phi2=levelset.phi(i.x,i.y+1,i.z);
        int index=vertices.Standard_Index(i),index_i=index-yz,index_ij=index-1,index_i_ij=index_i-1;
        if(phi1<=0&&phi2>0)topology.Append(VECTOR<int,4>(index_i_ij,index_i,index,index_ij));
        else if(phi1>0&&phi2<=0)topology.Append(VECTOR<int,4>(index_ij,index,index_i,index_i_ij));}
    for(i.x=1;i.x<m-1;i.x++)for(i.y=1;i.y<n-1;i.y++)for(i.z=0;i.z<mn-1;i.z++){ // generate one triangle pair per z edge
        T phi1=levelset.phi(i),phi2=levelset.phi(i.x,i.y,i.z+1);
        int index=vertices.Standard_Index(i),index_i=index-yz,index_j=index-z,index_i_j=index_i-z;
        if(phi1<=0&&phi2>0)topology.Append(VECTOR<int,4>(index_i_j,index_j,index,index_i));
        else if(phi1>0&&phi2<=0)topology.Append(VECTOR<int,4>(index_i,index,index_j,index_i_j));}
    topology.Compact();
    for(int t=0;t<topology.m;t++){
        int i,j,k,l;topology(t).Get(i,j,k,l);
        vertices.array(i)=vertices.array(j)=vertices.array(k)=vertices.array(l)=1;}
}
//#####################################################################
// Function Generate_Vertices
//#####################################################################
template<class T> void DUALCONTOUR_3D<T>::
Generate_Vertices()
{
    int count=vertices.Sum();
    geometry.Resize(count);
    normals.Resize(count);
    levelset.Compute_Normals();
    int vertex=-1;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),grid.counts-1));it.Valid();it.Next())
        if(vertices(it.index)){ // generate vertices where needed
            vertices(it.index)=++vertex;
            TV position=grid.Center(it.index);
            TV normal=levelset.Normal(position);
            T phi=levelset.Phi(position);
            int iterations=0;
            while(abs(phi)>1e-5*grid.dX.Min() && (iterations++)<10){
                position-=phi*normal;
                phi=levelset.Phi(position);
                normal=levelset.Normal(position);}
            geometry(vertex)=position;normals(vertex)=normal;}
        else vertices(it.index)=-1;
}
//#####################################################################
// Function Ensure_Vertices_In_Correct_Cells
//#####################################################################
template<class T> void DUALCONTOUR_3D<T>::
Ensure_Vertices_In_Correct_Cells()
{
    int vertex=-1;
    TV_INT i;
    for(i.x=0;i.x<grid.counts.x-1;i.x++)for(i.y=0;i.y<grid.counts.y-1;i.y++)for(i.z=0;i.z<grid.counts.z-1;i.z++)if(vertices(i)>=0){
        vertex++;TV_INT v=grid.Clamp_To_Cell(geometry(vertex));
        if(i!=v){
            TV cell_center=grid.Center(i);TV offset((T).5*grid.dX);
            geometry(vertex)=RANGE<TV>(cell_center-offset,cell_center+offset).Surface(geometry(vertex));
        }}
}
//#####################################################################
// Function Get_Triangulated_Surface
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* DUALCONTOUR_3D<T>::
Get_Triangulated_Surface()
{
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    GEOMETRY_PARTICLES<TV>& particles=surface->particles;
    particles.Add_Elements(geometry.m);
    particles.X=geometry;
    ARRAY<TV>* vertex_normals=new ARRAY<TV>(normals);
    TRIANGLE_MESH& mesh=surface->mesh;
    mesh.number_nodes=geometry.m;
    mesh.elements.Exact_Resize(2*topology.m);
    int current_triangle=0;
    for(int t=0;t<topology.m;t++){
        int i,j,k,l;topology(t).Get(i,j,k,l);i=vertices.array(i);j=vertices.array(j);k=vertices.array(k);l=vertices.array(l);
        mesh.elements(current_triangle++).Set(i,j,l);mesh.elements(current_triangle++).Set(j,k,l);}
    surface->Update_Triangle_List();surface->Use_Vertex_Normals();surface->vertex_normals=vertex_normals;
    return surface;
}
//#####################################################################
namespace PhysBAM{
template class DUALCONTOUR_3D<float>;
template class DUALCONTOUR_3D<double>;
}

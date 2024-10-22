//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// TESSELLATION
//##################################################################### 
#include <Core/Math_Tools/constants.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

namespace PhysBAM{
namespace TESSELLATION{
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const SPHERE<VECTOR<T,3> >& sphere,int levels)
{
    typedef VECTOR<T,3> TV;
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    GEOMETRY_PARTICLES<TV>& particles=surface->particles;particles.Add_Elements(6);
    particles.X(0)=TV(-1,0,0);particles.X(1)=TV(1,0,0);particles.X(2)=TV(0,-1,0);
    particles.X(3)=TV(0,1,0);particles.X(4)=TV(0,0,-1);particles.X(5)=TV(0,0,1);
    ARRAY<VECTOR<int,3> >& triangles=surface->mesh.elements;triangles.Exact_Resize(8);
    triangles(0).Set(0,5,3);triangles(1).Set(0,2,5);triangles(2).Set(5,1,3);triangles(3).Set(5,2,1);
    triangles(4).Set(4,0,3);triangles(5).Set(4,2,0);triangles(6).Set(1,2,4);triangles(7).Set(1,4,3);
    surface->mesh.number_nodes=6;
    surface->mesh.Initialize_Neighbor_Nodes();
    for(int i=0;i<levels;i++) surface->Root_Three_Subdivide();
    for(int p=0;p<particles.Size();p++) particles.X(p)=sphere.center+sphere.radius*particles.X(p).Normalized();
    return surface;
}
template<class T> TRIANGULATED_AREA<T>* Generate_Triangles(const SPHERE<VECTOR<T,2> >& circle,int levels)
{
    assert(levels>=1);
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,3> E;
    TRIANGULATED_AREA<T>* area=TRIANGULATED_AREA<T>::Create();
    GEOMETRY_PARTICLES<TV>& particles=area->particles;particles.Add_Elements(1+3*levels*(levels+1));
    particles.X(0)=TV(0,0);
    for(int i=0,k=1;i<levels;i++)
        for(int j=0;j<i*6+6;j++){
            T a=(T)pi/3*(j+1-(T).5*i)/(i+1);
            particles.X(k++)=(T)(i+1)*TV(cos(a),sin(a));}
    for(int i=1;i<7;i++) area->mesh.elements.Append(E(0,i,i%6+1));
    for(int i=1,p=1;i<levels;i++){
        int n=i*6+6;
        int c=p+n-6;
        area->mesh.elements.Append(E(c-1,c+n-1,c));
        area->mesh.elements.Append(E(c,p,c-1));
        for(int j=0;j<n-1;j++){
            area->mesh.elements.Append(E(p,c+j,c+j+1));
            if(j%(i+1)!=i/2){
                area->mesh.elements.Append(E(c+j+1,p+1,p));
                p++;}}
        p=c;}
    particles.X=particles.X*(circle.radius/levels)+circle.center;
    area->Update_Number_Nodes();
    return area;
}
template<class T> SEGMENTED_CURVE_2D<T>* Tessellate_Boundary(const SPHERE<VECTOR<T,2> >& sphere,int levels)
{
    assert(levels>=1);
    SEGMENTED_CURVE_2D<T>* curve=SEGMENTED_CURVE_2D<T>::Create();
    int n=1<<levels;
    curve->particles.Add_Elements(n);
    for(int i=0;i<n;i++) curve->particles.X(i)=VECTOR<T,2>((T)cos((T)pi*2*i/n),(T)sin((T)pi*2*i/n))*sphere.radius+sphere.center;
    for(int i=0;i<n-1;i++) curve->mesh.elements.Append(VECTOR<int,2>(i,i+1));
    curve->mesh.elements.Append(VECTOR<int,2>(n-1,0));
    curve->Update_Number_Nodes();
    return curve;
}
//#####################################################################
template TRIANGULATED_SURFACE<float>* Generate_Triangles(const SPHERE<VECTOR<float,3> >&,int);
template TRIANGULATED_AREA<float>* Generate_Triangles(const SPHERE<VECTOR<float,2> >&,int);
template SEGMENTED_CURVE_2D<float>* Tessellate_Boundary(const SPHERE<VECTOR<float,2> >& sphere,int levels);
template TRIANGULATED_SURFACE<double>* Generate_Triangles(const SPHERE<VECTOR<double,3> >&,int);
template TRIANGULATED_AREA<double>* Generate_Triangles(const SPHERE<VECTOR<double,2> >&,int);
template SEGMENTED_CURVE_2D<double>* Tessellate_Boundary(const SPHERE<VECTOR<double,2> >& sphere,int levels);
}
}

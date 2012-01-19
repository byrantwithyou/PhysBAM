//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// TESSELLATION
//##################################################################### 
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Geometry/Basic_Geometry/SMOOTH_GEAR.h>
#include <PhysBAM_Geometry/Tessellation/GEAR_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{template<class TV> void Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color);}
namespace PhysBAM{
namespace TESSELLATION{
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> void Boundary_Points(ARRAY<VECTOR<T,2> >& pts,const SMOOTH_GEAR<VECTOR<T,2> >& gear,int n)
{
    typedef VECTOR<T,2> TV;
    TV dc=gear.Co-gear.Ci;
    T hai=atan2(dc.y,dc.x),hao=hai-gear.den/2;
    T cs=cos(gear.den),sn=sin(gear.den);
    MATRIX<T,2> R(cs,sn,-sn,cs);
    int i=0,m=(int)ceil(n*hai/(hai+hao));
    T ds=2*(hai+hao)/n;
    for(;i<m;i++){T a=i*ds-hai;pts.Append(gear.Ci+gear.s*TV(cos(a),sin(a)));}
    for(;i<n;i++){T a=pi+gear.den/2-(i*ds-2*hai-hao);pts.Append(gear.Co+gear.s*TV(cos(a),sin(a)));}
    for(int j=1;j<=(gear.n-1)*n;j++) pts.Append(R*pts(j));
}

template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const SMOOTH_GEAR<VECTOR<T,3> >& gear,int n)
{
    typedef VECTOR<T,3> TV;
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    typedef VECTOR<int,3> E;
    ARRAY<VECTOR<T,2> > pts;
    Boundary_Points(pts, gear.g, n);
    GEOMETRY_PARTICLES<TV>& particles=surface->particles;particles.array_collection->Add_Elements(pts.m*2+2);
    for(int i=0;i<pts.m;i++){
        particles.X(i)=pts(i).Append(-gear.w);
        particles.X(i+pts.m)=pts(i).Append(gear.w);}
    particles.X(2*pts.m+1)=TV(0,0,-gear.w);
    particles.X(2*pts.m+2)=TV(0,0,gear.w);
    for(int i=1;i<pts.m;i++){
        surface->mesh.elements.Append(E(i,i+1,2*pts.m+1));
        surface->mesh.elements.Append(E(pts.m+i+1,pts.m+i,2*pts.m+2));
        surface->mesh.elements.Append(E(i+1,i,pts.m+i+1));
        surface->mesh.elements.Append(E(pts.m+i+1,i,pts.m+i));}
    surface->mesh.elements.Append(E(pts.m,1,2*pts.m+1));
    surface->mesh.elements.Append(E(2*pts.m,pts.m+1,2*pts.m+2));
    surface->mesh.elements.Append(E(1,pts.m,2*pts.m));
    surface->mesh.elements.Append(E(2*pts.m,pts.m,pts.m+1));
    surface->Update_Number_Nodes();
    return surface;
}
template<class T> TRIANGULATED_AREA<T>* Generate_Triangles(const SMOOTH_GEAR<VECTOR<T,2> >& gear,int n)
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,3> E;
    TRIANGULATED_AREA<T>* area=TRIANGULATED_AREA<T>::Create();
    ARRAY<VECTOR<T,2> > pts;
    Boundary_Points(pts, gear, n);
    GEOMETRY_PARTICLES<TV>& particles=area->particles;particles.array_collection->Add_Elements(pts.m);
    area->particles.X=pts;
    particles.array_collection->Add_Element();
    area->particles.X.Last()=TV();
    for(int i=1;i<pts.m;i++) area->mesh.elements.Append(E(i,i+1,area->particles.X.m));
    area->mesh.elements.Append(E(pts.m,1,area->particles.X.m));
    area->Update_Number_Nodes();
    return area;
}
template<class T> SEGMENTED_CURVE_2D<T>* Tessellate_Boundary(const SMOOTH_GEAR<VECTOR<T,2> >& gear,int n)
{
    SEGMENTED_CURVE_2D<T>* curve=SEGMENTED_CURVE_2D<T>::Create();
    ARRAY<VECTOR<T,2> > pts;
    Boundary_Points(pts, gear, n);
    curve->particles.array_collection->Add_Elements(pts.m);
    curve->particles.X=pts;
    for(int i=1;i<pts.m;i++) curve->mesh.elements.Append(VECTOR<int,2>(i,i+1));
    curve->mesh.elements.Append(VECTOR<int,2>(pts.m,1));
    curve->Update_Number_Nodes();
    return curve;
}
//#####################################################################
template TRIANGULATED_SURFACE<float>* Generate_Triangles(const SMOOTH_GEAR<VECTOR<float,3> >&,int);
template TRIANGULATED_AREA<float>* Generate_Triangles(const SMOOTH_GEAR<VECTOR<float,2> >&,int);
template SEGMENTED_CURVE_2D<float>* Tessellate_Boundary(const SMOOTH_GEAR<VECTOR<float,2> >& gear,int levels);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template TRIANGULATED_SURFACE<double>* Generate_Triangles(const SMOOTH_GEAR<VECTOR<double,3> >&,int);
template TRIANGULATED_AREA<double>* Generate_Triangles(const SMOOTH_GEAR<VECTOR<double,2> >&,int);
template SEGMENTED_CURVE_2D<double>* Tessellate_Boundary(const SMOOTH_GEAR<VECTOR<double,2> >& gear,int levels);
#endif
}
}

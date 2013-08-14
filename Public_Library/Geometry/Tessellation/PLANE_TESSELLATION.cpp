//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLANE_TESSELLATION
//#####################################################################
#include <Geometry/Basic_Geometry/LINE_2D.h>
#include <Geometry/Basic_Geometry/PLANE.h>
#include <Geometry/Tessellation/PLANE_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

namespace PhysBAM{
namespace TESSELLATION{
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>*
Generate_Triangles(const PLANE<T>& plane,T tesselated_size,int elements_per_side)
{
    typedef VECTOR<T,3> TV;
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    int m=elements_per_side+1,n=elements_per_side;
    surface->particles.Add_Elements(m*m);
    surface->mesh.Initialize_Herring_Bone_Mesh(m,m);
    TV u=tesselated_size*plane.normal.Unit_Orthogonal_Vector(),v=TV::Cross_Product(plane.normal,u);
    for(int j=0,p=0;j<=n;j++)
        for(int i=0;i<=n;i++)
            surface->particles.X(p++)=plane.x1+u*((T)2*j/n-1)-v*((T)2*i/n-1);
    surface->Update_Number_Nodes();
    return surface;
}
//#####################################################################
// Function Tessellate_Boundary
//#####################################################################
template<class T> SEGMENTED_CURVE_2D<T>*
Tessellate_Boundary(const LINE_2D<T>& line,T tesselated_size,int elements_per_side)
{
    typedef VECTOR<T,2> TV;
    SEGMENTED_CURVE_2D<T>* curve=SEGMENTED_CURVE_2D<T>::Create();
    int m=elements_per_side+1,n=elements_per_side;
    curve->mesh.Initialize_Straight_Mesh(m);
    curve->particles.Add_Elements(m);
    TV u=tesselated_size*line.normal.Orthogonal_Vector();
    for(int i=0;i<=n;i++)
        curve->particles.X(i)=line.x1+u*((T)2*i/n-1);
    curve->Update_Number_Nodes();
    return curve;
}
//#####################################################################
template TRIANGULATED_SURFACE<double>* Generate_Triangles<double>(PLANE<double> const&,double,int);
template TRIANGULATED_SURFACE<float>* Generate_Triangles<float>(PLANE<float> const&,float,int);
template SEGMENTED_CURVE_2D<double>* Tessellate_Boundary<double>(LINE_2D<double> const&,double,int);
template SEGMENTED_CURVE_2D<float>* Tessellate_Boundary<float>(LINE_2D<float> const&,float,int);
}
}

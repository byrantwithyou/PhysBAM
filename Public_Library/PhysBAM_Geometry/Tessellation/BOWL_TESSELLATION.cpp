//#####################################################################
// Copyright 2011, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOWL_TESSELLATION
//##################################################################### 
#include <PhysBAM_Tools/Matrices/MATRIX_3X2.h>
#include <PhysBAM_Tools/Vectors/COMPLEX.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOWL.h>
#include <PhysBAM_Geometry/Tessellation/BOWL_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

namespace PhysBAM{
namespace TESSELLATION{
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const BOWL<T>& bowl,const int n)
{
    typedef VECTOR<T,3> TV;
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    surface->mesh.Initialize_Torus_Mesh(4,n);
    GEOMETRY_PARTICLES<TV>& particles=surface->particles;particles.array_collection->Add_Elements(4*n);
    MATRIX<T,3,2> radial_basis;
    radial_basis.Column(1)=TV(0,1,0);
    radial_basis.Column(2)=TV(0,0,1);
    for(int i=1,p=0;i<=n;++i){
        TV radial=radial_basis*COMPLEX<T>::Unit_Polar(T(2*pi/n)*i).Vector();
        for (int j=0; j<=n; j++){
            TV spherical=radial_basis*COMPLEX<T>::Unit_Polar(T(0.5*pi/n)*j).Vector();
            particles.X(++p)=bowl.outer_radius*radial*spherical.y;
            particles.X(p).x=-bowl.outer_radius*spherical.z;}
        for (int j=n; j>=0; j--){
            TV spherical=radial_basis*COMPLEX<T>::Unit_Polar(T(0.5*pi/n)*j).Vector();
            particles.X(++p)=bowl.inner_radius*radial*spherical.y;
            particles.X(p).x=-bowl.inner_radius*spherical.z;}
    }
    return surface;
}
//#####################################################################
template TRIANGULATED_SURFACE<float>* Generate_Triangles(const BOWL<float>&,const int);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template TRIANGULATED_SURFACE<double>* Generate_Triangles(const BOWL<double>&,const int);
#endif
}
}

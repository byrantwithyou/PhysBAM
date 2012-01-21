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
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const BOWL<T>& bowl,const int n_radial,const int n_vertical)
{
    typedef VECTOR<T,3> TV;
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    surface->mesh.Initialize_Torus_Mesh(2*(n_vertical+1),n_radial);
    GEOMETRY_PARTICLES<TV>& particles=surface->particles;particles.array_collection->Add_Elements(2*(n_vertical+1)*n_radial);
    MATRIX<T,3,2> radial_basis;
    radial_basis.Column(0)=TV(0,0,1);
    radial_basis.Column(1)=TV(1,0,0);
    for(int i=0,p=0;i<n_radial;++i){
        TV radial=radial_basis*COMPLEX<T>::Unit_Polar(T(2*pi/n_radial)*i).Vector();
        for (int j=n_vertical; j>=0; j--){
            TV spherical=radial_basis*COMPLEX<T>::Unit_Polar(T(0.5*pi/n_vertical)*j).Vector();
            particles.X(p++)=radial*(bowl.hole_radius+spherical.z*bowl.height);
            particles.X(p).y=bowl.height*spherical.x;}
        for (int j=0; j<=n_vertical; j++){
            TV spherical=radial_basis*COMPLEX<T>::Unit_Polar(T(0.5*pi/n_vertical)*j).Vector();
            particles.X(p++)=radial*(bowl.hole_radius+spherical.z*bowl.depth);
            particles.X(p).y=bowl.depth*spherical.x;}
    }
    return surface;
}
//#####################################################################
template TRIANGULATED_SURFACE<float>* Generate_Triangles(const BOWL<float>&,const int,const int);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template TRIANGULATED_SURFACE<double>* Generate_Triangles(const BOWL<double>&,const int,const int);
#endif
}
}

//#####################################################################
// Copyright 2011, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOWL_TESSELLATION
//##################################################################### 
#include <Core/Matrices/MATRIX_3X2.h>
#include <Geometry/Basic_Geometry/BOWL.h>
#include <Geometry/Tessellation/BOWL_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <complex>

namespace PhysBAM{
namespace TESSELLATION{
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const BOWL<T>& bowl,const int n_radial,const int n_vertical)
{
    typedef VECTOR<T,3> TV;
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    int hres=2*n_vertical,m=0;
    if(!bowl.hole_radius){
        hres-=2;
        m=1;}
    surface->mesh.Initialize_Cylinder_Mesh(hres,n_radial);
    GEOMETRY_PARTICLES<TV>& particles=surface->particles;
    particles.Add_Elements(hres*n_radial+2);
    MATRIX<T,3,2> radial_basis;
    radial_basis.Set_Column(0,TV(0,0,1));
    radial_basis.Set_Column(1,TV(1,0,0));
    TV depth_axis(0,1,0);
    int p=0;
    for(int j=n_vertical-1-m;j>=0;j--){
        std::complex<T> c=std::polar((T)1,T(0.5*pi/(n_vertical-1))*j);
        TV spherical=radial_basis*VECTOR<T,2>(c.real(),c.imag());
        for(int i=0;i<n_radial;i++){
            std::complex<T> c=std::polar((T)1,T(2*pi/n_radial)*i);
            TV radial=radial_basis*VECTOR<T,2>(c.real(),c.imag());
            particles.X(p++)=radial*(bowl.hole_radius+spherical.z*bowl.height)+bowl.height*spherical.x*depth_axis;}}
    for(int j=0;j<=n_vertical-1-m;j++){
        std::complex<T> c=std::polar((T)1,T(0.5*pi/(n_vertical-1))*j);
        TV spherical=radial_basis*VECTOR<T,2>(c.real(),c.imag());
        for(int i=0;i<n_radial;i++){
            std::complex<T> c=std::polar((T)1,T(2*pi/n_radial)*i);
            TV radial=radial_basis*VECTOR<T,2>(c.real(),c.imag());
            particles.X(p++)=radial*(bowl.hole_radius+spherical.z*bowl.depth)+bowl.depth*spherical.x*depth_axis;}}
    particles.X(p++)=bowl.height*depth_axis;
    particles.X(p++)=bowl.depth*depth_axis;
    PHYSBAM_ASSERT(p==particles.number);
    PHYSBAM_ASSERT(p==surface->mesh.number_nodes);
    return surface;
}
//#####################################################################
template TRIANGULATED_SURFACE<float>* Generate_Triangles(const BOWL<float>&,const int,const int);
template TRIANGULATED_SURFACE<double>* Generate_Triangles(const BOWL<double>&,const int,const int);
}
}

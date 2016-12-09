//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RING_TESSELLATION
//##################################################################### 
#include <Core/Matrices/MATRIX_3X2.h>
#include <Geometry/Basic_Geometry/RING.h>
#include <Geometry/Tessellation/RING_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <complex>

namespace PhysBAM{
namespace TESSELLATION{
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const RING<T>& ring,const int n)
{
    typedef VECTOR<T,3> TV;
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    surface->mesh.Initialize_Torus_Mesh(4,n);
    GEOMETRY_PARTICLES<TV>& particles=surface->particles;particles.Add_Elements(4*n);
    MATRIX<T,3,2> radial_basis;
    radial_basis.Set_Column(0,ring.plane1.normal.Orthogonal_Vector());
    radial_basis.Set_Column(1,TV::Cross_Product(ring.plane1.normal,radial_basis.Column(0)));
    for(int i=0,p=0;i<n;++i){
        std::complex<T> c=std::polar((T)1,T(2*pi/n)*i);
        TV radial=radial_basis*VECTOR<T,2>(c.real(),c.imag());
        particles.X(p++)=ring.plane1.x0+ring.inner_radius*radial;
        particles.X(p++)=ring.plane1.x0+ring.outer_radius*radial;
        particles.X(p++)=ring.plane2.x0+ring.outer_radius*radial;
        particles.X(p++)=ring.plane2.x0+ring.inner_radius*radial;}
    return surface;
}
//#####################################################################
template TRIANGULATED_SURFACE<float>* Generate_Triangles(const RING<float>&,const int);
template TRIANGULATED_SURFACE<double>* Generate_Triangles(const RING<double>&,const int);
}
}

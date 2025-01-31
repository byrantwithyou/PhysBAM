//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#include <Core/Matrices/MATRIX_3X3.h>
#include <Geometry/Basic_Geometry/TORUS.h>
#include <Geometry/Tessellation/TORUS_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

namespace PhysBAM{
namespace TESSELLATION{
//#####################################################################
// Function Generate_Triangles
//#####################################################################
// m is inner_resolution, n is outer_resolution
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const TORUS<T>& torus,const int m,const int n)
{
    typedef VECTOR<T,3> TV;
    T mscale=(T)pi*2/m,nscale=(T)pi*2/n;
    TV y=torus.axis.Unit_Orthogonal_Vector();
    MATRIX<T,3> transform(y,TV::Cross_Product(torus.axis,y),torus.axis);
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    GEOMETRY_PARTICLES<TV>& particles=surface->particles;particles.Add_Elements(m*n);
    TRIANGLE_MESH& mesh=surface->mesh;mesh.number_nodes=m*n;mesh.elements.Preallocate(2*m*n);
    for(int i=0;i<m;++i){
        T axial=torus.inner_radius*cos(i*mscale),inplane=torus.inner_radius*sin(i*mscale)+torus.outer_radius;
        for(int j=0;j<n;j++){
            T cj=cos(j*nscale),sj=sin(j*nscale);
            particles.X(i+j*m)=transform*TV(cj*inplane,sj*inplane,axial);
            int ni=(i+1)%m,nj=(j+1)%n,p00=i+j*m,p01=i+nj*m,p10=ni+j*m,p11=p10+p01-p00;
            mesh.elements.Append(VECTOR<int,3>(p00,p10,p01));
            mesh.elements.Append(VECTOR<int,3>(p11,p01,p10));}}
    return surface;
}
//#####################################################################
template TRIANGULATED_SURFACE<float>* Generate_Triangles(const TORUS<float>&,const int,const int);
template TRIANGULATED_SURFACE<double>* Generate_Triangles(const TORUS<double>&,const int,const int);
}
}

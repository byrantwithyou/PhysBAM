//#####################################################################
// Copyright 2005, Geoffrey Irving, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STRAIN_MEASURE_HEXAHEDRONS  
//#####################################################################
#ifndef __STRAIN_MEASURE_HEXAHEDRONS__
#define __STRAIN_MEASURE_HEXAHEDRONS__

#include <Core/Matrices/MATRIX_3X3.h>
#include <Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
namespace PhysBAM{

template<class T>
class STRAIN_MEASURE_HEXAHEDRONS
{
    typedef VECTOR<T,3> TV;
public:
    HEXAHEDRALIZED_VOLUME<T>& mesh_object;
    HEXAHEDRON_MESH& mesh;
    DEFORMABLE_PARTICLES<TV>& particles;
    ARRAY<VECTOR<ARRAY<TV>,8>  > H_DmH_inverse; // 8x3 matrix per gauss point per element
    ARRAY<VECTOR<T,8> > DmH_determinant;
    T DmH_minimum_altitude; // this is only vaguely an altitude
private:
    ARRAY<VECTOR<ARRAY<TV>,8> >* H_DmH_inverse_save;
    ARRAY<ARRAY<TV> > H; // 8x3 matrix per gauss point
public:
    
    STRAIN_MEASURE_HEXAHEDRONS(HEXAHEDRALIZED_VOLUME<T>& mesh_object_input)
        :mesh_object(mesh_object_input),mesh(mesh_object.mesh),particles(dynamic_cast<DEFORMABLE_PARTICLES<TV>&>(mesh_object.particles)),
            H_DmH_inverse_save(0)
    {
        Initialize_H();
        Initialize_H_DmH_Inverse_From_Current_Positions();
    }

    STRAIN_MEASURE_HEXAHEDRONS(const STRAIN_MEASURE_HEXAHEDRONS&) = delete;
    void operator=(const STRAIN_MEASURE_HEXAHEDRONS&) = delete;
    
    ~STRAIN_MEASURE_HEXAHEDRONS()
    {delete H_DmH_inverse_save;}

    void Initialize_H_DmH_Inverse_From_Current_Positions()
    {Initialize_H_DmH_Inverse(particles.X);}

    void Initialize_H_DmH_Inverse(ARRAY_VIEW<const TV> X)
    {H_DmH_inverse.Resize(8,mesh.elements.m);DmH_determinant.Resize(8,mesh.elements.m);DmH_minimum_altitude=FLT_MAX;
    for(int h=0;h<mesh.elements.m;h++)for(int g=0;g<8;g++){
        H_DmH_inverse(g,h).Resize(8);
        MATRIX<T,3> DmH;for(int k=0;k<8;k++)DmH+=MATRIX<T,3>::Outer_Product(X(mesh.elements(h)(k)),H(g)(k));
        DmH_determinant(g,h)=DmH.Determinant();DmH_minimum_altitude=min(DmH_minimum_altitude,DmH.Tetrahedron_Minimum_Altitude());
        MATRIX<T,3> DmH_inverse=DmH.Inverse();for(int k=0;k<8;k++)H_DmH_inverse(g,h)(k)=H(g)(k)*DmH_inverse;}}
    
    void Initialize_H_DmH_Inverse_Save()
    {H_DmH_inverse_save=new ARRAY<VECTOR<TV,8> >(mesh.elements.m);
    ARRAY<VECTOR<TV,8> >::copy(H_DmH_inverse,*H_DmH_inverse_save);}

    void Copy_H_DmH_Inverse_Save_Into_H_DmH_Inverse(const ARRAY<int>& map)
    {H_DmH_inverse.Resize(map.m);ARRAY<VECTOR<ARRAY<TV>,8> >::permute(*H_DmH_inverse_save,H_DmH_inverse,map);}

    MATRIX<T,3> Gradient(ARRAY_VIEW<const TV> X,const ARRAY<VECTOR<ARRAY<TV>,8> >& De_inverse,const int gauss_index,const int hexahedron_index) const
    {MATRIX<T,3> gradient;for(int k=0;k<8;k++)gradient+=MATRIX<T,3>::Outer_Product(X(mesh.elements(hexahedron_index)(k)),De_inverse(hexahedron_index)(gauss_index)(k));return gradient;}

    MATRIX<T,3> F(const int gauss_index,const int hexahedron_index) const
    {return Gradient(particles.X,H_DmH_inverse,gauss_index,hexahedron_index);}

    MATRIX<T,3> dF(ARRAY_VIEW<const TV> dX,const int gauss_index,const int hexahedron_index) const
    {return Gradient(dX,H_DmH_inverse,gauss_index,hexahedron_index);}

    MATRIX<T,3> Velocity_Gradient(const int gauss_index,const int hexahedron_index) const
    {return Gradient(particles.V,H_DmH_inverse,gauss_index,hexahedron_index);}

private:
    void Initialize_H()
    {H.Resize(8);
    for(int i=0;i<2;i++)for(int j=0;j<2;j++)for(int ij=0;ij<2;ij++){
        int gauss_index=4*i+2*j+ij+1;H(gauss_index).Resize(8);
        TV xi=1/sqrt((T)3)*TV((T)(2*i-1),(T)(2*j-1),(T)(2*ij-1));
        for(int ii=0;ii<2;ii++)for(int jj=0;jj<2;jj++)for(int ijij=0;ijij<2;ijij++){
            int node=4*ii+2*jj+ijij+1;TV s((T)(2*ii-1),(T)(2*jj-1),(T)(2*ijij-1)),sxi=s*xi;
            H(gauss_index)(node)=(T).125*s*TV((1+sxi.y)*(1+sxi.z),(1+sxi.x)*(1+sxi.z),(1+sxi.x)*(1+sxi.y));}}}

//#####################################################################
};
}
#endif

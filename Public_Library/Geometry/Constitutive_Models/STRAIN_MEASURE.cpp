//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Igor Neverov, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STRAIN_MEASURE
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Matrices/MATRIX_2X2.h>
#include <Core/Matrices/MATRIX_3X2.h>
#include <Core/Matrices/MATRIX_3X3.h>
#include <Core/Matrices/UPPER_TRIANGULAR_MATRIX_2X2.h>
#include <Core/Matrices/UPPER_TRIANGULAR_MATRIX_3X3.h>
#include <Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> STRAIN_MEASURE<TV,d>::
STRAIN_MEASURE(T_MESH_OBJECT& mesh_object)
    :mesh_object(mesh_object),mesh(mesh_object.mesh),particles(dynamic_cast<GEOMETRY_PARTICLES<TV>&>(mesh_object.particles)),
    Dm_inverse_save(0)
{
    Initialize_Dm_Inverse(mesh_object.particles.X);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int d> STRAIN_MEASURE<TV,d>::
~STRAIN_MEASURE()
{
    delete Dm_inverse_save;
}
//#####################################################################
// Function Initialize_Dm_Inverse
//#####################################################################
template<class TV,int d> void STRAIN_MEASURE<TV,d>::
Initialize_Dm_Inverse(ARRAY_VIEW<const TV> X)
{
    Dm_inverse.Resize(mesh.elements.m,false,false);
    for(int t=0;t<mesh.elements.m;t++){
        UPPER_TRIANGULAR_MATRIX<T,d> R=Ds(X,t).R_From_QR_Factorization();
        if(R.Determinant()<=0) PHYSBAM_FATAL_ERROR("Inverted or degenerate rest state");
        Dm_inverse(t)=R.Inverse();}
}
//#####################################################################
// Function Initialize_Dm_Inverse_Save
//#####################################################################
template<class TV,int d> void STRAIN_MEASURE<TV,d>::
Initialize_Dm_Inverse_Save()
{
    Dm_inverse_save=new ARRAY<UPPER_TRIANGULAR_MATRIX<T,d> >;
    *Dm_inverse_save=Dm_inverse;
}
//#####################################################################
// Function Copy_Dm_Inverse_Save_Into_Dm_Inverse
//#####################################################################
template<class TV,int d> void STRAIN_MEASURE<TV,d>::
Copy_Dm_Inverse_Save_Into_Dm_Inverse(const ARRAY<int>& map)
{
    Dm_inverse.Resize(map.m,false,false);
    Dm_inverse=Dm_inverse_save->Subset(map);
}
//#####################################################################
// Function Initialize_Rest_State_To_Equilateral
//#####################################################################
namespace{
template<class T> UPPER_TRIANGULAR_MATRIX<T,2> Equilateral_Dm(const VECTOR<T,2>&)
{
    return UPPER_TRIANGULAR_MATRIX<T,2>(1,.5,sqrt((T)3)/2);
}
template<class T> UPPER_TRIANGULAR_MATRIX<T,3> Equilateral_Dm(const VECTOR<T,3>&)
{
    T x=sqrt((T)3)/3,d=(T).5*x,h=sqrt((T)6)/3;
    return MATRIX<T,3>(x,0,-h,-d,.5,-h,-d,-.5,-h).R_From_QR_Factorization();
}}
template<class TV,int d> void STRAIN_MEASURE<TV,d>::
Initialize_Rest_State_To_Equilateral(const T side_length)
{
    Dm_inverse.Resize(mesh.elements.m,false,false);
    UPPER_TRIANGULAR_MATRIX<T,d> Dm=side_length*Equilateral_Dm(VECTOR<T,d>());
    Dm_inverse.Fill(Dm.Inverse());
}
//#####################################################################
// Function Print_Altitude_Statistics
//#####################################################################
template<class TV,int d> void STRAIN_MEASURE<TV,d>::
Print_Altitude_Statistics()
{   
    if(!Dm_inverse.m) return;
    ARRAY<T> altitude(Dm_inverse.m,false);
    for(int t=0;t<altitude.m;t++) altitude(t)=Rest_Altitude(t);
    altitude.Sort();
    LOG::cout<<"strain measure - total elements = "<<altitude.m<<std::endl;
    LOG::cout<<"strain measure - smallest altitude = "<<altitude(0)<<std::endl;
    LOG::cout<<"strain measure - one percent altitude = "<<altitude((int)(.01*altitude.m))<<std::endl;
    LOG::cout<<"strain measure - ten percent altitude = "<<altitude((int)(.1*altitude.m))<<std::endl;
    LOG::cout<<"strain measure - median altitude = "<<altitude((int)(.5*altitude.m))<<std::endl;
}
//#####################################################################
#define INSTANTIATION_HELPER_T(T,m,d) \
    template STRAIN_MEASURE<VECTOR<T,m>,d>::STRAIN_MEASURE(T_MESH_OBJECT&); \
    template STRAIN_MEASURE<VECTOR<T,m>,d>::~STRAIN_MEASURE(); \
    template void STRAIN_MEASURE<VECTOR<T,m>,d>::Initialize_Dm_Inverse(ARRAY_VIEW<const VECTOR<T,m> >); \
    template void STRAIN_MEASURE<VECTOR<T,m>,d>::Initialize_Dm_Inverse_Save(); \
    template void STRAIN_MEASURE<VECTOR<T,m>,d>::Copy_Dm_Inverse_Save_Into_Dm_Inverse(const ARRAY<int>&); \
    template void STRAIN_MEASURE<VECTOR<T,m>,d>::Initialize_Rest_State_To_Equilateral(const T); \
    template void STRAIN_MEASURE<VECTOR<T,m>,d>::Print_Altitude_Statistics();
#define INSTANTIATION_HELPER(T) INSTANTIATION_HELPER_T(T,2,2) INSTANTIATION_HELPER_T(T,3,2) INSTANTIATION_HELPER_T(T,3,3)
INSTANTIATION_HELPER(float)
INSTANTIATION_HELPER(double)

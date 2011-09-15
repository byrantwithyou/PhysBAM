//#####################################################################
// Copyright 2003, Zhaosheng Bao.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ALTITUDE_SPRING_TEST
//##################################################################### 
//
// Hanging and stretching test with Altitude Springs.
//
//#####################################################################
// Bao - August 12, 2003
//#####################################################################
#ifndef __ALTITUDE_SPRING_TEST__
#define __ALTITUDE_SPRING_TEST__

#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_S3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include <fstream>
#include "../CLOTH_EXAMPLE.h"
namespace PhysBAM{

template<class T>
class ALTITUDE_SPRING_TEST:public CLOTH_EXAMPLE<T> 
{
public:
    int number_side_panels;
    bool use_altitude_springs;
    T aspect_ratio,side_length;
private:
    LINEAR_SPRINGS<T>* ls;
    LINEAR_ALTITUDE_SPRINGS_S3D<T>* las;
    TRIANGLE_BENDING_ELEMENTS<T>* bend;
    T cfl_number;

public:
    ALTITUDE_SPRING_TEST()
                                    :number_side_panels(1),aspect_ratio((T)1),side_length(1),                      
                                     ls(0),bend(0),las(0),use_altitude_springs(true),cfl_number((T).9)
    {
        final_time=7;
        restart_step_number=0;
        output_directory="Altitude_Spring_Test/output";
        check_initial_mesh_for_self_intersection=false;
    }

    ~ALTITUDE_SPRING_TEST()
    {}

//#####################################################################
// Function Initialize_Cloth_State
//#####################################################################
void Initialize_Cloth(DEFORMABLE_TRIANGULATED_SURFACE<T>& cloth)
{
    if(restart_step_number) Read_Deformable_Triangulated_Surface(cloth,restart_step_number);
    else{
        int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
        cloth.triangulated_surface.triangle_mesh.Initialize_Square_Mesh(m,n);
        for(int k=1;k<=cloth.triangulated_surface.triangle_mesh.number_nodes;k++) cloth.triangulated_surface.particles.array_collection->Add_Element();
        T mass_node=aspect_ratio*sqr(side_length)/(m*n);ARRAY<T>::copy(mass_node,cloth.triangulated_surface.particles.mass.array);
        T dx=aspect_ratio*side_length/(m-1),dy=side_length/(n-1);
        for(int i=1;i<=m;i++) for(int j=1;j<=n;j++){int node=i+m*(j-1);
            cloth.triangulated_surface.particles.X(node)=VECTOR_3D<T>((i-1)*dx,.5,(j-1)*dy);
            cloth.triangulated_surface.particles.V(node)=VECTOR_3D<T>(0,0,0);}}

    // cloth.body_forces.Set_Gravity();
    cloth.Set_CFL_Number((T)cfl_number);
    cloth.Output_Artificial_Damping_Results();
    
    // edge springs 
    delete cloth.triangulated_surface.triangle_mesh.segment_mesh;cloth.triangulated_surface.triangle_mesh.Initialize_Segment_Mesh();
    ls=new LINEAR_SPRINGS<T>(*cloth.triangulated_surface.triangle_mesh.segment_mesh,cloth.triangulated_surface.particles);
    cloth.mesh_based_forces.Append(ls);
    ls->Set_Restlength_From_Particles();
    ls->Set_Stiffness((2*1/(1+(T)root_two)));
    ls->Set_Overdamping_Fraction(2);
    ls->Limit_Time_Step_By_Strain_Rate(true,(T).1);
    ls->Use_Rest_State_For_Strain_Rate();
    ls->Artificially_Damp_Max_Strain_Per_Time_Step();
    ls->Artificially_Damp_Min_And_Max_Strain(true,(T)-.01,(T).1);

    // linear altitude springs
    if(use_altitude_springs){
        las=new LINEAR_ALTITUDE_SPRINGS_S3D<T>(cloth.triangulated_surface.triangle_mesh,cloth.triangulated_surface.particles);
        cloth.mesh_based_forces.Append(las);
        las->Set_Restlength_From_Particles();
        las->Set_Stiffness(10*(4*1/(1+(T)root_two)));
        las->Set_Overdamping_Fraction(2);
        las->Limit_Time_Step_By_Strain_Rate(true,(T).1);
        las->Use_Rest_State_For_Strain_Rate();
        las->Artificially_Damp_Max_Strain_Per_Time_Step();
        las->Artificially_Damp_Min_And_Max_Strain(true,(T)-0.1,(T)0.1);}

    // bending
    bend=new TRIANGLE_BENDING_ELEMENTS<T>(cloth.triangulated_surface.particles);
    cloth.mesh_based_forces.Append(bend);
    bend->Set_Quadruples_From_Triangle_Mesh(cloth.triangulated_surface.triangle_mesh);
    bend->Set_Constants_From_Particles((T).001,(T).001);
    //bend->Set_Stiffness(.001);bend->Set_Damping(.001);

    // set up repulsion springs for collisions!!!!!!
    repulsion_springs_initialized=true;
    T average_restlength=ARRAY<T>::sum(ls->restlength)/ls->restlength.m; 
    T mass_node=cloth.triangulated_surface.particles.mass.Total_Mass(cloth.triangulated_surface.particles)/cloth.triangulated_surface.particles.array_collection->Size();
    collisions_repulsion_spring_constant_over_mass_times_length=2*ls->constant_youngs_modulus/(mass_node*average_restlength);
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY<VECTOR_3D<T> >& V,const T time)
{    
    int i,j;
    int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
    i=1;j=1;V(i+m*(j-1))=VECTOR_3D<T>(-1.5,0,0);
    i=1;j=n;V(i+m*(j-1))=VECTOR_3D<T>(1.5,0,0);
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> >& V,const T time)
{
    int i,j;
    int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
    i=1;j=1;V(i+m*(j-1))=VECTOR_3D<T>(0,0,0);
    i=1;j=n;V(i+m*(j-1))=VECTOR_3D<T>(0,0,0);
}
//#####################################################################
};
}
#endif


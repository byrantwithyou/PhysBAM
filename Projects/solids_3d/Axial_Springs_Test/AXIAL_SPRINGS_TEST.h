//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class AXIAL_SPRINGS_TEST
//##################################################################### 
// Hanging curtain with a moving ball
#ifndef __AXIAL_SPRINGS_TEST__
#define __AXIAL_SPRINGS_TEST__

#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>
namespace PhysBAM{

template<class T,class RW>
class AXIAL_SPRINGS_TEST:public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
public:
    int number_side_panels;
    T aspect_ratio,side_length;

    AXIAL_SPRINGS_TEST()
        :SOLIDS_FLUIDS_EXAMPLE_3D<RW>(fluids_parameters.NONE),number_side_panels(4),aspect_ratio((T)1.7),side_length(1)
    {
        //allow_intersections=true;allow_intersections_tolerance=(T)1e-3;
        last_frame=7*24;
        restart=false;restart_frame=0;
        solids_parameters.cfl=(T)0.5;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
        output_directory="Axial_Springs_Test/output";
        //min_collision_loops=max_collision_loops=1;
    }

    ~AXIAL_SPRINGS_TEST()
    {}

//#####################################################################
// Function Initialize_Body
//#####################################################################
void Initialize_Body(const bool use_axial_springs,const VECTOR_3D<T>& offset)
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Triangulated_Surface();
    TRIANGULATED_SURFACE<T>& triangulated_surface=*solids_parameters.deformable_body_parameters.list(index).triangulated_surface;
    TRIANGLE_MESH& triangle_mesh=triangulated_surface.triangle_mesh;
    PARTICLES<T,VECTOR_3D<T> >& particles=triangulated_surface.particles;

#if 1
    triangle_mesh.number_nodes=4;
    triangle_mesh.triangles.Exact_Resize(3,2);
    for(int i=0;i<4;i++){particles.array_collection->Add_Element();particles.V(i)=VECTOR_3D<T>();}
    particles.X(1)=offset+VECTOR_3D<T>(1,10,0);
    particles.X(2)=offset+VECTOR_3D<T>(0,10,-1);
    particles.X(3)=offset+VECTOR_3D<T>(0,10,1);
    particles.X(4)=offset+VECTOR_3D<T>(-1,10,0);
    triangle_mesh.triangles.Set(1,1,2,3);
    triangle_mesh.triangles.Set(2,3,2,4);
    ARRAY<T>::copy(1,particles.mass.array);
#else
    int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
    triangle_mesh.Initialize_Herring_Bone_Mesh(m,n);
    for(int k=0;k<triangle_mesh.number_nodes;k++) particles.array_collection->Add_Element();
    T mass_node=aspect_ratio*sqr(side_length)/(m*n);ARRAY<T>::copy(mass_node,particles.mass.array); 
    T dx=aspect_ratio*side_length/(m-1),dy=side_length/(n-1);
    for(int i=0;i<m;i++) for(int j=0;j<n;j++){int node=i+m*(j-1);
        particles.X(node)=offset+VECTOR_3D<T>((i-1)*dx,.5,(j-1)*dy);
        particles.V(node)=VECTOR_3D<T>(0,0,0);}
#endif

    solids_parameters.deformable_body_parameters.list(index).Add_Body_Forces(*solids_parameters.deformable_body_parameters.list(index).triangulated_surface);
//    solids_parameters.deformable_body_parameters.list(index).Add_Edge_Springs(solids_parameters.deformable_body_parameters.list(index).triangulated_surface->triangle_mesh,2/(index+sqrt((T)2)),2);
//    solids_parameters.deformable_body_parameters.list(index).Add_Altitude_Springs(solids_parameters.deformable_body_parameters.list(index).triangulated_surface->triangle_mesh,2*4/(index+sqrt((T)2)),4);
//    solids_parameters.deformable_body_parameters.list(index).Add_Bending_Elements(solids_parameters.deformable_body_parameters.list(index).triangulated_surface->triangle_mesh);
    //solids_parameters.deformable_body_parameters.list(index).Add_Bending_Springs(solids_parameters.deformable_body_parameters.list(index).triangulated_surface->triangle_mesh);
    if(use_axial_springs) solids_parameters.deformable_body_parameters.list(index).Add_Axial_Bending_Springs(solids_parameters.deformable_body_parameters.list(index).triangulated_surface->triangle_mesh,.01,1);
}
//#####################################################################
// Function Initialize_Deformable_Objects
//#####################################################################
void Initialize_Bodies()
{
    Initialize_Body(true,VECTOR_3D<T>());
//    Initialize_Body(false,VECTOR_3D<T>(4,0,0));
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Initialize_Bodies();
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
// for external forces and velocities
void Set_External_Velocities(ARRAY<VECTOR_3D<T> >& V,const T time)
{
    Zero_Out_Enslaved_Velocity_Nodes(V,time,id_number);
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
// for external forces and velocities
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> >& V,const T time)
{
#if 0
    int i,j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
    for(int i=1;i<=m/2+1;i++) for(int j=0;j<n;j++)V(i+m*(j-1))=VECTOR_3D<T>(0,0,0);i=1;j=n;
#endif
    V(2)=V(3)=VECTOR_3D<T>();
}
//#####################################################################
};
}
#endif


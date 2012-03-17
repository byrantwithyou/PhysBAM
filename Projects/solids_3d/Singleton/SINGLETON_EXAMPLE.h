//#####################################################################
// Copyright 2003-2004, Geoffrey Irving, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SINGLETON_EXAMPLE
//#####################################################################
#ifndef __SINGLETON_EXAMPLE__
#define __SINGLETON_EXAMPLE__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <Solids_And_Fluids/SOLIDS_3D_EXAMPLE.h>
namespace PhysBAM{
template<class T,class RW>
class SINGLETON_EXAMPLE:public SOLIDS_3D_EXAMPLE<T,RW>
{
public:
    T initial_height;
    T edge_length;//simplex from cube with this edge length
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    bool constrain_face;
    
    SINGLETON_EXAMPLE()
        :initial_height(4),edge_length(1),initial_orientation(0,VECTOR_3D<T>(1,0,0)),
        initial_velocity(0,0,0),initial_angular_velocity(5,2,0),constrain_face(true)
    {
        last_frame=10*24;
        restart=false;restart_frame=0;   
        cg_tolerance=(T)1e-2;
        output_directory="Singleton/output";
    }

    ~SINGLETON_EXAMPLE()
    {}

//#####################################################################
// Function Initialize_Deformable_Object
//#####################################################################
void Initialize_Deformable_Objects()
{
    SOLIDS_3D_EXAMPLE<T,RW>::Initialize_Deformable_Objects();

    solids_parameters.deformable_body_parameters.list(1).Add_Body_Forces(*solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume);
    solids_parameters.deformable_body_parameters.list(1).Add_Neo_Hookean_Elasticity(*solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume,(T)1e3);
    //solids_parameters.deformable_body_parameters.list(1).Add_Linear_Finite_Volume(*solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume,3e4);
    //solids_parameters.deformable_body_parameters.list(1).Add_Diagonalized_Spline_Model(*solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume,(T)1e5,(T).475,0,1,(T).05);
    //solids_parameters.deformable_body_parameters.list(1).Add_Edge_Springs(solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume->tetrahedron_mesh,3000,(T).9);
    //solids_parameters.deformable_body_parameters.list(1).Add_Altitude_Springs(solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume->tetrahedron_mesh,300);
    //solids_parameters.deformable_body_parameters.list(1).Add_Linear_Elasticity(tetrahedralized_volume,2e5,.45,.01);
    //solids_parameters.deformable_body_parameters.list(1).Add_Diagonalized_Neo_Hookean_Elasticity(tetrahedralized_volume,2e5,.45,.01);
}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    solids_parameters.deformable_body_parameters.list.Add_Deformable_Object();
    solids_parameters.deformable_body_parameters.list(1).Allocate_Tetrahedralized_Volume();
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume;
    TETRAHEDRON_MESH& tetrahedron_mesh=tetrahedralized_volume.tetrahedron_mesh;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles=solids_parameters.deformable_body_parameters.list(1).particles;

    particles.Increase_Array_Size(4);
    particles.X(particles.Add_Element())=VECTOR_3D<T>(0,0,0);
    particles.X(particles.Add_Element())=edge_length*VECTOR_3D<T>(1,0,0);
    particles.X(particles.Add_Element())=edge_length*VECTOR_3D<T>(0,1,0);
    particles.X(particles.Add_Element())=edge_length*VECTOR_3D<T>(0,0,1);
    tetrahedron_mesh.Clean_Memory();tetrahedron_mesh.number_nodes=4;tetrahedron_mesh.tetrahedrons.Exact_Resize(4,1);
    tetrahedron_mesh.tetrahedrons(1,1)=1;tetrahedron_mesh.tetrahedrons(2,1)=2;tetrahedron_mesh.tetrahedrons(3,1)=3;tetrahedron_mesh.tetrahedrons(4,1)=4;

    tetrahedralized_volume.Set_Density(1000);
    tetrahedralized_volume.Set_Mass_Of_Particles(false);

    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=0;i<tetrahedralized_volume.particles.array_size;i++){
        tetrahedralized_volume.particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i)=center+initial_orientation.Rotate(tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i).y+=initial_height-bottom;}
}
//#####################################################################
// Function Initialize_Collision_Bodies
//#####################################################################
void Initialize_Collision_Bodies()
{  
    if(!constrain_face){
    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>("../../Public_Data/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(solids_parameters.rigid_body_parameters.list.rigid_bodies.m)->coefficient_of_friction=(T).3;}
    solids_parameters.collision_body_list.Append_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY<VECTOR_3D<T> >& V,const T time){
    switch(id_number){
    case 1:
        if(constrain_face){V(1)=V(2)=V(3)=VECTOR_3D<T>(0,0,0);}
        break;
    default:std::cout<<"Unrecognized deformable object id number"<<std::endl;exit(1);}
}
//#####################################################################
// Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> >& V,const T time){
    switch(id_number){
    case 1:
        if(constrain_face){V(1)=V(2)=V(3)=VECTOR_3D<T>(0,0,0);}
        break;
    default:std::cout<<"Unrecognized deformable object id number"<<std::endl;exit(1);}
}
//#####################################################################
};
}
#endif

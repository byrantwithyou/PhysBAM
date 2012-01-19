//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BICEPS_EXAMPLE
//#####################################################################
#ifndef __BICEPS_EXAMPLE__
#define __BICEPS_EXAMPLE__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <Solids_And_Fluids/SOLIDS_3D_EXAMPLE.h>
namespace PhysBAM{

template<class T,class RW>
class BICEPS_EXAMPLE:public SOLIDS_3D_EXAMPLE<T,RW>
{
public:
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    ARRAY<bool> tendon_tets;
    ARRAY<ARRAY<int> > bone_attachments;
    std::string muscle_name;
    bool constrain_attachments;

    BICEPS_EXAMPLE()
        :initial_height(1),initial_orientation(0,VECTOR_3D<T>(1,0,0)),initial_velocity(0,0,0),initial_angular_velocity(0,0,0),constrain_attachments(true)
    {   
        collisions_repulsion_thickness=(T)1e-2;
        collisions_repulsion_clamp_fraction=(T).9;
        collision_repulsion_spring_multiplier=100;

        last_frame=10*24;
        restart=false;restart_frame=16;   
        cfl=(T).75;
        cg_tolerance=(T)1e-2;
        perform_self_collision=false;
        solids_parameters.collide_with_interior=false;
        output_directory="Biceps/output";
        muscle_name="anconeus_tendon_right";
    }

    ~BICEPS_EXAMPLE()
    {}

    void Postprocess_Substep(const T time,const int substep)
    {std::cout << "minimum volume = " << solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume->Minimum_Signed_Volume() << std::endl;}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    solids_parameters.deformable_body_parameters.list.Add_Deformable_Object();
    DEFORMABLE_OBJECT_3D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(1);
    deformable_object.Allocate_Tetrahedralized_Volume();
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*deformable_object.tetrahedralized_volume;
    TETRAHEDRON_MESH& tetrahedron_mesh=tetrahedralized_volume.tetrahedron_mesh;
    PARTICLES<T,VECTOR_3D<T> >& particles=deformable_object.particles;

    std::string input_file=data_directory+"/VH_Muscles/"+muscle_name+".tet";
    std::ifstream input;FILE_UTILITIES::Safe_Open(input,input_file);tetrahedralized_volume.template Read<RW>(input);input.close();
    std::cout << "total vertices = " << particles.array_collection->Size() << std::endl;std::cout << "total tetrahedra = " << tetrahedron_mesh.tetrahedrons.m << std::endl;
    particles.Store_Velocity(false);particles.Store_Velocity(true);particles.Update_Velocity();
    particles.Store_Mass(false);particles.Store_Mass(true);

    tetrahedralized_volume.Set_Density(1000);tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);
    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=0;i<particles.array_collection->Size();i++){
        particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,particles.X(i)-center);
        particles.X(i)=center+initial_orientation.Rotate(particles.X(i)-center);
        particles.X(i).y+=initial_height-bottom;}

    Initialize_Bone_Attachments();
}
//#####################################################################
// Function Initialize_Deformable_Object
//#####################################################################
void Initialize_Deformable_Objects()
{
    SOLIDS_3D_EXAMPLE<T,RW>::Initialize_Deformable_Objects();
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume;

    solids_parameters.deformable_body_parameters.list(1).Add_Body_Forces(tetrahedralized_volume);
    solids_parameters.deformable_body_parameters.list.deformable_objects(1)->Add_Linear_Finite_Volume(tetrahedralized_volume,(T)1e3,(T).275);
    //solids_parameters.deformable_body_parameters.list(1).Add_Edge_Springs(tetrahedralized_volume.tetrahedron_mesh,3000,2);
    //solids_parameters.deformable_body_parameters.list(1).Add_Altitude_Springs(tetrahedralized_volume.tetrahedron_mesh,300);
    //solids_parameters.deformable_body_parameters.list(1).Add_Linear_Elasticity(tetrahedralized_volume,2e5,.45,.01);
    //solids_parameters.deformable_body_parameters.list(1).Add_Neo_Hookean_Elasticity(tetrahedralized_volume,(T)2e5,(T).45,(T).01);
    //solids_parameters.deformable_body_parameters.list(1).Add_Diagonalized_Neo_Hookean_Elasticity(tetrahedralized_volume,(T)2e5,(T).45,(T).01);
    //solids_parameters.deformable_body_parameters.list(1).Add_Diagonalized_Linear_Finite_Volume(tetrahedralized_volume,(T)3e5,(T).3,(T).01);
    //solids_parameters.deformable_body_parameters.list(1).Disable_Finite_Volume_Damping();
    //solids_parameters.deformable_body_parameters.list(1).Add_Edge_Springs(tetrahedralized_volume.tetrahedron_mesh,3000,(T)1);
    //solids_parameters.deformable_body_parameters.list(1).Add_Altitude_Springs(tetrahedralized_volume.tetrahedron_mesh,300,(T)1);
    //solids_parameters.deformable_body_parameters.list(1).Disable_Spring_Elasticity();
}
//#####################################################################
// Function Initialize_Bone_Attachments
//#####################################################################
void Initialize_Bone_Attachments()
{
    ARRAY<std::string *> bone_input_files(6); 
    bone_input_files(1)=new std::string(data_directory+"/VH_Muscles/Attachments/"+muscle_name+"__sternum.attach");bone_input_files(2)=new std::string(data_directory+"/VH_Muscles/Attachments/"+muscle_name+"__clavicle_right.attach");
    bone_input_files(3)=new std::string(data_directory+"/VH_Muscles/Attachments/"+muscle_name+"__scapula_right.attach");bone_input_files(4)=new std::string(data_directory+"/VH_Muscles/Attachments/"+muscle_name+"__humerus_right.attach");
    bone_input_files(5)=new std::string(data_directory+"/VH_Muscles/Attachments/"+muscle_name+"__ulna_right.attach");bone_input_files(6)= new std::string(data_directory+"/VH_Muscles/Attachments/"+muscle_name+"__radius_right.attach");
    bone_attachments.Resize(6);
    for(int b=0;b<bone_input_files.m;b++){
        std::fstream input;input.open(bone_input_files(b)->c_str(),std::ios::in|std::ios::binary);
        if(input.is_open())bone_attachments(b).template Read<RW>(input);
        input.close();}
}
//#####################################################################
// Function Initialize_Collision_Bodies
//#####################################################################
void Initialize_Collision_Bodies()
{  
    if(constrain_attachments) return;
    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(solids_parameters.rigid_body_parameters.list.rigid_bodies.m)->coefficient_of_friction=(T).3;
    solids_parameters.collision_body_list.Append_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY<VECTOR_3D<T> >& V,const T time){
    if(!constrain_attachments) return;
    switch(id_number){
    case 1:
        for(int b=0;b<bone_attachments.m;b++)for(int t=0;t<bone_attachments(b).m;t++){
            int i,j,k,l;solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume->tetrahedron_mesh.tetrahedrons.Get(bone_attachments(b)(t),i,j,k,l);
            V(i)=V(j)=V(k)=V(l)=VECTOR_3D<T>(0,0,0);}
        break;
    default:std::cout<<"Unrecognized deformable object id number"<<std::endl;exit(1);}
}
//#####################################################################
// Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> >& V,const T time){
    if(!constrain_attachments) return;
    switch(id_number){
    case 1:
        for(int b=0;b<bone_attachments.m;b++)for(int t=0;t<bone_attachments(b).m;t++){
            int i,j,k,l;solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume->tetrahedron_mesh.tetrahedrons.Get(bone_attachments(b)(t),i,j,k,l);
            V(i)=V(j)=V(k)=V(l)=VECTOR_3D<T>(0,0,0);}
        break;
    default:std::cout<<"Unrecognized deformable object id number"<<std::endl;exit(1);}
}
//#####################################################################
};
}
#endif

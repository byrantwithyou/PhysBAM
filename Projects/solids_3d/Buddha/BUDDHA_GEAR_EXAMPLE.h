//#####################################################################
// Copyright 2003-2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BUDDHA_GEAR_EXAMPLE
//#####################################################################
#ifndef __BUDDHA_GEAR_EXAMPLE__
#define __BUDDHA_GEAR_EXAMPLE__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <Solids_And_Fluids/SOLIDS_3D_EXAMPLE.h>
namespace PhysBAM{

template<class T,class RW>
class BUDDHA_GEAR_EXAMPLE:public SOLIDS_3D_EXAMPLE<T,RW>
{
public:
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    std::string input_file;
    int gear_1,gear_2,cylinder;
    T cylinder_time;
    VECTOR_3D<T> cylinder_start,cylinder_velocity;
    QUATERNION<T> roller_orientation;
    T roller_speed,roller_friction;
    
    BUDDHA_GEAR_EXAMPLE()
        :initial_height((T)2.234),initial_orientation(T(pi/2),VECTOR_3D<T>(0,0,1)),initial_velocity(),initial_angular_velocity()
    {
        frame_rate=120;
        first_frame=0;last_frame=10*120;
        
        solids_parameters.use_constant_mass=true; // gah
        restart=true;restart_frame=222;
        cfl=10;
        cg_iterations=1000;
        cg_tolerance=(T)1e-2;

        roller_speed=(T)2;
        roller_friction=1;
        roller_orientation=QUATERNION<T>(0,VECTOR_3D<T>(1,0,0));

        cylinder_time=(T).2;
        cylinder_start=VECTOR_3D<T>(0,2,0);
        cylinder_velocity=VECTOR_3D<T>(0,-5,0);

        input_file=data_directory+"/Tetrahedralized_Volumes/buddha.tet";
        output_directory="Buddha/output_gear";
        
        check_initial_mesh_for_self_intersection=false;
        perform_self_collision=true;
        solids_parameters.collide_with_interior=true;
        max_collision_loops=1;

/*
        collision_repulsion_spring_multiplier=1e10;
        collisions_repulsion_thickness=(T)3e-3;
        collisions_repulsion_clamp_fraction=(T).25;
        collisions_collision_thickness=(T)1e-5;
        collisions_nonrigid_collision_attempts=8;
*/
    }

    ~BUDDHA_GEAR_EXAMPLE()
    {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Tetrahedralized_Volume();
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*solids_parameters.deformable_body_parameters.list(index).tetrahedralized_volume;
    TETRAHEDRON_MESH& tetrahedron_mesh=tetrahedralized_volume.tetrahedron_mesh;
    PARTICLES<T,VECTOR_3D<T> >& particles=tetrahedralized_volume.particles;
 
    std::ifstream input;FILE_UTILITIES::Safe_Open(input,input_file);tetrahedralized_volume.template Read<RW>(input);input.close();
    std::cout << "total vertices = " << particles.array_collection->Size() << std::endl;std::cout << "total tetrahedra = " << tetrahedron_mesh.tetrahedrons.m << std::endl;
    particles.Update_Velocity();particles.Store_Mass(); // in case they're not stored in the file!

    tetrahedralized_volume.Set_Density(1000);tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);
    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());
    for(int i=1;i<=particles.array_size;i++){
        particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,particles.X(i)-center);
        particles.X(i)=center+initial_orientation.Rotate(particles.X(i)-center);
        particles.X(i).y+=initial_height;}
}
//#####################################################################
// Function Initialize_Deformable_Object
//#####################################################################
void Initialize_Deformable_Objects()
{
    SOLIDS_3D_EXAMPLE<T,RW>::Initialize_Deformable_Objects();
    DEFORMABLE_OBJECT_3D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(1);
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*deformable_object.tetrahedralized_volume;

    solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_object.particles));
    solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(tetrahedralized_volume,new DIAGONALIZED_SPLINE_MODEL_3D<T>(tetrahedralized_volume,50000,(T).45,0,10,(T).03),true,(T).15));
}    
//#####################################################################
// Function Initialize_Collision_Bodies
//#####################################################################
void Initialize_Collision_Bodies()
{  
    RIGID_BODY<TV>* rigid_body;
    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/ground");
    rigid_body=solids_parameters.rigid_body_parameters.list.rigid_bodies(solids_parameters.rigid_body_parameters.list.rigid_bodies.m);
    rigid_body->coefficient_of_friction=(T)1;

    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/Rings_Test/cylinder_revolve",.375/2);
    cylinder=solids_parameters.rigid_body_parameters.list.rigid_bodies.m;
    rigid_body=solids_parameters.rigid_body_parameters.list.rigid_bodies(cylinder);
    rigid_body->position=cylinder_start;rigid_body->velocity=cylinder_velocity;
    rigid_body->orientation=QUATERNION<T>(T(pi/2),VECTOR_3D<T>(1,0,0));
    rigid_body->coefficient_of_friction=0;

    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/gear",.375);
    gear_1=solids_parameters.rigid_body_parameters.list.rigid_bodies.m;
    rigid_body=solids_parameters.rigid_body_parameters.list.rigid_bodies(gear_1);
    rigid_body->position=VECTOR_3D<T>(-(T).4,1.5,-.75);
    rigid_body->orientation=roller_orientation;
    rigid_body->angular_velocity=-roller_speed*VECTOR_3D<T>(0,0,1);
    rigid_body->coefficient_of_friction=roller_friction;

    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/gear",.375);
    gear_2=solids_parameters.rigid_body_parameters.list.rigid_bodies.m;
    rigid_body=solids_parameters.rigid_body_parameters.list.rigid_bodies(gear_2);
    rigid_body->position=VECTOR_3D<T>((T).4,1.5,-.75);
    rigid_body->orientation=roller_orientation;
    rigid_body->angular_velocity=roller_speed*VECTOR_3D<T>(0,0,1);
    rigid_body->coefficient_of_friction=roller_friction;

    solids_parameters.collision_body_list.Append_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Update_Collision_Body_Positions_And_Velocities
//#####################################################################
void Update_Collision_Body_Positions_And_Velocities(const T time) PHYSBAM_OVERRIDE
{
    solids_parameters.rigid_body_parameters.list.rigid_bodies(gear_1)->orientation=QUATERNION<T>(-roller_speed*time,VECTOR_3D<T>(0,0,1))*roller_orientation;
    solids_parameters.rigid_body_parameters.list.rigid_bodies(gear_2)->orientation=QUATERNION<T>(roller_speed*time,VECTOR_3D<T>(0,0,1))*roller_orientation;
    if(time>cylinder_time){
        solids_parameters.rigid_body_parameters.list.rigid_bodies(cylinder)->position=cylinder_start;
        solids_parameters.rigid_body_parameters.list.rigid_bodies(cylinder)->velocity=VECTOR_3D<T>();}
    else{
        solids_parameters.rigid_body_parameters.list.rigid_bodies(cylinder)->position=cylinder_start+(time-cylinder_time)*cylinder_velocity;
        solids_parameters.rigid_body_parameters.list.rigid_bodies(cylinder)->velocity=cylinder_velocity;}
}
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE
{
    return;
    ARRAY<bool>& check_collision=solids_parameters.deformable_body_parameters.list(1).collisions.check_collision;
    DIAGONALIZED_FINITE_VOLUME_3D<T>& fvm=*solids_parameters.deformable_body_parameters.list(1).diagonalized_finite_volume_3d(1);
    if(!fvm.Fe_hat.m)return;
    static ARRAY<int> timeout;timeout.Resize(fvm.Fe_hat.m);
    ARRAY<bool>::copy(false,check_collision);
    for(int t=1;t<=fvm.Fe_hat.m;t++){
        if(fvm.Fe_hat(t).x11<3){ // && timeout(t)--<=0){
            int i,j,k,l;fvm.strain_measure.tetrahedron_mesh.tetrahedrons.Get(t,i,j,k,l);
            check_collision(i)=check_collision(j)=check_collision(k)=check_collision(l)=true;}
        else{timeout(t)=10;/*std::cout<<"############ tetrahedron "<<t<<": "<<fvm.Fe_hat(t)<<"\n";*/}}
    for(int p=1;p<=check_collision.m;p++)if(!check_collision(p))std::cout<<"########### disabling collisions for particle "<<p<<"\n";
}
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2003-2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BUDDHA_EXAMPLE
//#####################################################################
#ifndef __BUDDHA_EXAMPLE__
#define __BUDDHA_EXAMPLE__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE.h>
namespace PhysBAM{

template<class T,class RW>
class BUDDHA_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE<VECTOR_3D<T>,RW>
{
public:
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    std::string input_file;
    bool standing_on_ground;
    bool use_gears,use_first_sphere,use_second_sphere,use_cylinder;
    int gear_1,gear_2,sphere_1,sphere_2,cylinder;
    T ball_1_time,ball_2_time,cylinder_time;
    T ball_1_radius,ball_2_radius;
    VECTOR_3D<T> ball_1_start,ball_1_velocity;
    VECTOR_3D<T> ball_2_start,ball_2_velocity;
    VECTOR_3D<T> cylinder_start,cylinder_velocity;
    QUATERNION<T> roller_orientation;
    T roller_speed,roller_friction;
    T hardening_deformation,hardening_strength;
    int example;
    bool joey;
    
    T max_strain_per_time_step;
    T youngs_modulus,poissons_ratio,Rayleigh_coefficient;

    BUDDHA_EXAMPLE()
        :initial_height(0),initial_orientation(0,VECTOR_3D<T>(1,0,0)),initial_velocity(),initial_angular_velocity(),
        use_gears(false),use_first_sphere(false),use_second_sphere(false),use_cylinder(false),standing_on_ground(true),
        ball_1_time(0),ball_2_time(0),cylinder_time(0),ball_1_radius(.25),ball_2_radius(.25)
    {
        example=4;
        
        frame_rate=120;
        first_frame=0;last_frame=10*120;
        max_strain_per_time_step=(T).1;
        
        if(example==1){
            solids_parameters.use_constant_mass=false;
            restart=false;restart_frame=37;
            cfl=10;
            cg_iterations=1000;
            cg_tolerance=(T)1e-2;
            youngs_modulus=25000;
            poissons_ratio=(T).3;
            Rayleigh_coefficient=(T).03;
            hardening_deformation=0;
            hardening_strength=7;
            use_first_sphere=true;
            ball_1_start=VECTOR_3D<T>((T).05,1,-1);
            ball_1_velocity=VECTOR_3D<T>(0,0,15);}
        else if(example==2){
            solids_parameters.use_constant_mass=false;
            restart=false;restart_frame=66;
            cfl=10;
            cg_iterations=1000;
            cg_tolerance=(T)1e-2;
            youngs_modulus=50000;
            poissons_ratio=(T).3;
            Rayleigh_coefficient=(T).03;
            hardening_deformation=0;
            hardening_strength=7;
            use_first_sphere=true;
            ball_1_start=VECTOR_3D<T>((T).05,1,-1);
            ball_1_velocity=VECTOR_3D<T>(0,0,10);}
        else if(example==3){
            solids_parameters.use_constant_mass=true;
            restart=false;restart_frame=0;
            cfl=10;
            cg_iterations=1000;
            cg_tolerance=(T)1e-2;
            youngs_modulus=50000;
            poissons_ratio=(T).45;
            Rayleigh_coefficient=(T).03;
            hardening_deformation=-.5;
            hardening_strength=7;
            use_first_sphere=true;
            ball_1_start=VECTOR_3D<T>((T).05,1,-1);
            ball_1_velocity=VECTOR_3D<T>(0,0,15);}
        else if(example==4){
            solids_parameters.use_constant_mass=true; // gah
            restart=false;restart_frame=108;
            cfl=10;
            cg_iterations=1000;
            cg_tolerance=(T)1e-2;
            youngs_modulus=50000;
            poissons_ratio=(T).45;
            Rayleigh_coefficient=(T).03;
            hardening_deformation=0;
            hardening_strength=10;
            use_gears=true;
            roller_speed=(T)2;
            roller_friction=1;
            use_cylinder=true;
            cylinder_time=(T).2;
            cylinder_start=VECTOR_3D<T>(0,2,0);
            cylinder_velocity=VECTOR_3D<T>(0,-5,0);
            standing_on_ground=false;
            initial_height=(T)2.234;
            initial_orientation=QUATERNION<T>(T(pi/2),VECTOR_3D<T>(0,0,1));}
        else{
            std::cout<<"This has no hope of working.  I might as well give up now.\n";exit(1);}
        
        input_file=data_directory+"/Tetrahedralized_Volumes/buddha.tet";
        output_directory=STRING_UTILITIES::string_sprintf("Buddha/output%d",example);
        
        check_initial_mesh_for_self_intersection=false;
        perform_self_collision=true;
        solids_parameters.collide_with_interior=true;
        max_collision_loops=1;

        collision_repulsion_spring_multiplier=1e10;
        collisions_repulsion_thickness=(T)3e-3;
        collisions_repulsion_clamp_fraction=(T).25;
        collisions_collision_thickness=(T)1e-5;
        collisions_nonrigid_collision_attempts=8;
    }

    ~BUDDHA_EXAMPLE()
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
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=0;i<particles.array_size;i++){
        particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,particles.X(i)-center);
        particles.X(i)=center+initial_orientation.Rotate(particles.X(i)-center);
        T height=standing_on_ground?-bottom:initial_height;
        particles.X(i).y+=height;}
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
    solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(tetrahedralized_volume,new DIAGONALIZED_SPLINE_MODEL_3D<T>(tetrahedralized_volume,
        youngs_modulus,poissons_ratio,hardening_deformation,hardening_strength,Rayleigh_coefficient),true,max_strain_per_time_step));
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

    if(use_first_sphere){
        solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/sphere",ball_1_radius);
        sphere_1=solids_parameters.rigid_body_parameters.list.rigid_bodies.m;
        rigid_body=solids_parameters.rigid_body_parameters.list.rigid_bodies(sphere_1);
        rigid_body->position=ball_1_start;
        rigid_body->velocity=ball_1_velocity;
        rigid_body->coefficient_of_friction=0;}

    if(use_second_sphere){
        solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/sphere",ball_2_radius);
        sphere_2=solids_parameters.rigid_body_parameters.list.rigid_bodies.m;
        rigid_body=solids_parameters.rigid_body_parameters.list.rigid_bodies(sphere_2);
        rigid_body->position=ball_2_start;
        rigid_body->velocity=ball_2_velocity;
        rigid_body->coefficient_of_friction=0;}
    
    if(use_cylinder){
        solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/Rings_Test/cylinder_revolve",.375/2);
        cylinder=solids_parameters.rigid_body_parameters.list.rigid_bodies.m;
        rigid_body=solids_parameters.rigid_body_parameters.list.rigid_bodies(cylinder);
        rigid_body->position=cylinder_start;rigid_body->velocity=cylinder_velocity;
        rigid_body->orientation=QUATERNION<T>(T(pi/2),VECTOR_3D<T>(1,0,0));
        rigid_body->coefficient_of_friction=0;}
    
    if(use_gears){
        roller_orientation=QUATERNION<T>(0,VECTOR_3D<T>(1,0,0));

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
        rigid_body->coefficient_of_friction=roller_friction;}

    solids_parameters.collision_body_list.Append_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Update_Collision_Body_Positions_And_Velocities
//#####################################################################
void Update_Collision_Body_Positions_And_Velocities(const T time) PHYSBAM_OVERRIDE
{
    if(use_first_sphere)solids_parameters.rigid_body_parameters.list.rigid_bodies(sphere_1)->position=ball_1_start+(time-ball_1_time)*ball_1_velocity;
    if(use_second_sphere)solids_parameters.rigid_body_parameters.list.rigid_bodies(sphere_2)->position=ball_2_start+(time-ball_2_time)*ball_2_velocity;
    if(use_gears){
        solids_parameters.rigid_body_parameters.list.rigid_bodies(gear_1)->orientation=QUATERNION<T>(-roller_speed*time,VECTOR_3D<T>(0,0,1))*roller_orientation;
        solids_parameters.rigid_body_parameters.list.rigid_bodies(gear_2)->orientation=QUATERNION<T>(roller_speed*time,VECTOR_3D<T>(0,0,1))*roller_orientation;}
    if(use_cylinder){
        if(time>cylinder_time){
            solids_parameters.rigid_body_parameters.list.rigid_bodies(cylinder)->position=cylinder_start;
            solids_parameters.rigid_body_parameters.list.rigid_bodies(cylinder)->velocity=VECTOR_3D<T>();}
        else{
            solids_parameters.rigid_body_parameters.list.rigid_bodies(cylinder)->position=cylinder_start+(time-cylinder_time)*cylinder_velocity;
            solids_parameters.rigid_body_parameters.list.rigid_bodies(cylinder)->velocity=cylinder_velocity;}}
}
//#####################################################################
};
}
#endif

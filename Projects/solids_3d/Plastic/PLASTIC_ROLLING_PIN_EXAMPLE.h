//#####################################################################
// Copyright 2003, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLASTIC_ROLLING_PIN_EXAMPLE
//#####################################################################
#ifndef __PLASTIC_ROLLING_PIN_EXAMPLE__
#define __PLASTIC_ROLLING_PIN_EXAMPLE__

#include "PLASTICITY_EXAMPLE.h"
namespace PhysBAM{

template<class T>
class PLASTIC_ROLLING_PIN_EXAMPLE:public PLASTICITY_EXAMPLE<T>
{
public:
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    QUATERNION<T> roller_orientation;
    T roller_speed,roller_friction;
    bool use_gears;

    PLASTIC_ROLLING_PIN_EXAMPLE()
        :initial_height((T)1.9),initial_orientation(0,VECTOR_3D<T>(1,0,0)),initial_velocity(0,0,0),initial_angular_velocity(0,0,0)
    {
        final_time=5;
        frame_rate=120;
        restart_step_number=0;
        cfl_number=(T)50;
        max_strain_per_time_step=(T).2;
        cg_tolerance=(T)1e-3;
        cg_iterations=1000;
        use_masses_and_springs=false;use_altitude_springs=false;
        use_diagonalized_fvm=true;
        use_linear_elasticity=false;use_neo_hookean=false;
        diagonalized_neo_hookean_failure_threshold=(T).1;
        youngs_modulus=500000;poissons_ratio=(T).45;Rayleigh_coefficient=(T).01;
        use_plasticity=false;
        use_control=false;
        preserve_volume=false;
        show_goal=false;
        yield_ratio=(T)1.5;
        plastic_clamp_ratio=5;
        use_gears=true;
        roller_speed=(T)4;
        roller_friction=0;//100000;
        strcpy(output_directory,"Plastic/output_elastic");
        strcpy(input_file,"../../Public_Data/Tetrahedralized_Volumes/sphere_finer.tet");
        //strcpy(input_file,"../../Public_Data/Tetrahedralized_Volumes/sphere.tet");
        check_initial_mesh_for_self_intersection=false;
        collisions_repulsion_thickness=(T)1e-3;
        perform_self_collision=true;
        solids_parameters.collide_with_interior=true;
        enforce_tangential_collision_velocities=false;
        max_collision_loops=1;
    }

    ~PLASTIC_ROLLING_PIN_EXAMPLE()
    {}
    
    
    
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    TETRAHEDRON_MESH& tetrahedron_mesh=tetrahedralized_volume.tetrahedron_mesh;
    HEAVY_PARTICLES<T,VECTOR_3D<T> >& particles=tetrahedralized_volume.particles;

    std::fstream input;input.open(input_file,std::ios::in|std::ios::binary);tetrahedralized_volume.Read_Float(input);input.close();
    std::cout << "total vertices = " << particles.array_collection->Size() << std::endl;std::cout << "total tetrhedra = " << tetrahedron_mesh.tetrahedrons.m << std::endl;
    tetrahedralized_volume.particles.Delete_Velocity_And_Acceleration();tetrahedralized_volume.particles.Delete_Mass(); // in case they're accidently stored in the .tet file
    tetrahedralized_volume.particles.Update_Position_And_Velocity();tetrahedralized_volume.particles.Store_Mass(); // add back with the proper sizes

    tetrahedralized_volume.Set_Density(1000);
    tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);
    
    if(use_control){
        plastic_goal.Exact_Resize(1,particles.array_collection->Size());
        T goal_thickness=(T).125,threshold=(T).9;
        MATRIX<T,3> Q=MATRIX<T,3>::Rotation_Matrix(VECTOR_3D<T>(1,0,0),VECTOR_3D<T>(1,0,0));
        for(int p=1;p<=particles.array_collection->Size();p++)plastic_goal(p)=Pancake_Map(particles.X(p),goal_thickness,threshold,Q);
        if(preserve_volume){
            T volume=tetrahedralized_volume.Total_Volume();
            ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >::Exchange_Arrays(particles.X,plastic_goal.array);
            T scale=pow(volume/tetrahedralized_volume.Total_Volume(),(T)one_third);
            for(int p=1;p<=particles.array_collection->Size();p++)particles.X(p)*=scale;
            ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >::Exchange_Arrays(particles.X,plastic_goal.array);}
        if(show_goal){ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >::Exchange_Arrays(particles.X,plastic_goal.array);use_control=false;}}

    for(int p=1;p<=particles.array_collection->Size();p++)particles.X(p)*=(T).5;
    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=1;i<=particles.array_collection->Size();i++){
        particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,particles.X(i)-center);
        particles.X(i)=center+initial_orientation.Rotate(particles.X(i)-center);
        particles.X(i).y+=initial_height-bottom;}
}
//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
void Initialize_Rigid_Bodies(ARRAY<RIGID_BODY<TV>*>& rigid_bodies)
{  
    std::fstream input;char filename[256];
    
    // plane
    rigid_bodies.Resize(rigid_bodies.m+1);rigid_bodies(rigid_bodies.m)=new RIGID_BODY<TV>;
    // triangulated surface
    int index=triangulated_surface_list.Add_Triangulated_Surface();
    sprintf(filename,"../../Public_Data/Rigid_Bodies/ground.tri");input.open(filename,std::ios::in|std::ios::binary);
    triangulated_surface_list.triangulated_surface(index)->Read_Float(input);input.close();
    rigid_bodies(rigid_bodies.m)->Initialize_Triangulated_Surface(*triangulated_surface_list.triangulated_surface(index));
    // implicit surface
    index=implicit_surface_list.Add_Levelset_Implicit_Surface();
    sprintf(filename,"../../Public_Data/Rigid_Bodies/ground.phi");input.open(filename,std::ios::in|std::ios::binary);
    implicit_surface_list.implicit_surface(index)->Read_Float(input);input.close();
    rigid_bodies(rigid_bodies.m)->Initialize_Implicit_Surface(*implicit_surface_list.implicit_surface(index));
    // rigid body
    sprintf(filename,"../../Public_Data/Rigid_Bodies/ground.rgd");input.open(filename,std::ios::in|std::ios::binary);
    rigid_bodies(rigid_bodies.m)->Read_Float(input);input.close();
    rigid_bodies(rigid_bodies.m)->coefficient_of_friction=0;
    
    if(use_gears){
        roller_orientation=QUATERNION<T>(0,VECTOR_3D<T>(1,0,0));
    
        // gear 1
        rigid_bodies.Resize(rigid_bodies.m+1);rigid_bodies(rigid_bodies.m)=new RIGID_BODY<TV>;
        // triangulated surface
        index=triangulated_surface_list.Add_Triangulated_Surface();
        sprintf(filename,"../../Public_Data/Rigid_Bodies/gear.tri");input.open(filename,std::ios::in|std::ios::binary);
        triangulated_surface_list.triangulated_surface(index)->Read_Float(input);input.close();
        triangulated_surface_list.triangulated_surface(index)->Rescale(.375); // resize radius from the default of 1 to .25
        rigid_bodies(rigid_bodies.m)->Initialize_Triangulated_Surface(*triangulated_surface_list.triangulated_surface(index));
        // implicit surface
        index=implicit_surface_list.Add_Levelset_Implicit_Surface();
        sprintf(filename,"../../Public_Data/Rigid_Bodies/gear.phi");input.open(filename,std::ios::in|std::ios::binary);
        implicit_surface_list.implicit_surface(index)->Read_Float(input);input.close();
        implicit_surface_list.implicit_surface(index)->Rescale(.375); // resize radius from the default of 1 to .25
        rigid_bodies(rigid_bodies.m)->Initialize_Implicit_Surface(*implicit_surface_list.implicit_surface(index));
        // rigid body
        rigid_bodies(rigid_bodies.m)->position=VECTOR_3D<T>(-(T).4,1.5,-.75); // reset position
        rigid_bodies(rigid_bodies.m)->orientation=roller_orientation;
        rigid_bodies(rigid_bodies.m)->angular_velocity=-roller_speed*VECTOR_3D<T>(0,0,1);
        rigid_bodies(rigid_bodies.m)->coefficient_of_friction=roller_friction;

        // gear 2
        rigid_bodies.Resize(rigid_bodies.m+1);rigid_bodies(rigid_bodies.m)=new RIGID_BODY<TV>;
        // triangulated surface
        index=triangulated_surface_list.Add_Triangulated_Surface();
        sprintf(filename,"../../Public_Data/Rigid_Bodies/gear.tri");input.open(filename,std::ios::in|std::ios::binary);
        triangulated_surface_list.triangulated_surface(index)->Read_Float(input);input.close();
        triangulated_surface_list.triangulated_surface(index)->Rescale(.375); // resize radius from the default of 1 to .25
        rigid_bodies(rigid_bodies.m)->Initialize_Triangulated_Surface(*triangulated_surface_list.triangulated_surface(index));
        // implicit surface
        index=implicit_surface_list.Add_Levelset_Implicit_Surface();
        sprintf(filename,"../../Public_Data/Rigid_Bodies/gear.phi");input.open(filename,std::ios::in|std::ios::binary);
        implicit_surface_list.implicit_surface(index)->Read_Float(input);input.close();
        implicit_surface_list.implicit_surface(index)->Rescale(.375); // resize radius from the default of 1 to .25
        rigid_bodies(rigid_bodies.m)->Initialize_Implicit_Surface(*implicit_surface_list.implicit_surface(index));
        // rigid body
        rigid_bodies(rigid_bodies.m)->position=VECTOR_3D<T>((T).4,1.5,-.75); // reset position
        rigid_bodies(rigid_bodies.m)->orientation=roller_orientation;
        rigid_bodies(rigid_bodies.m)->angular_velocity=roller_speed*VECTOR_3D<T>(0,0,1);
        rigid_bodies(rigid_bodies.m)->coefficient_of_friction=roller_friction;}
    else{
        roller_orientation=QUATERNION<T>((T)pi/2,VECTOR_3D<T>(1,0,0));
    
        // cylinder 1
        rigid_bodies.Resize(rigid_bodies.m+1);rigid_bodies(rigid_bodies.m)=new RIGID_BODY<TV>;
        // triangulated surface
        index=triangulated_surface_list.Add_Triangulated_Surface();
        sprintf(filename,"../../Public_Data/Rigid_Bodies/Rings_Test/cylinder_revolve.tri");input.open(filename,std::ios::in|std::ios::binary);
        triangulated_surface_list.triangulated_surface(index)->Read_Float(input);input.close();
        triangulated_surface_list.triangulated_surface(index)->Rescale(.5); // resize radius from the default of 1 to .25
        rigid_bodies(rigid_bodies.m)->Initialize_Triangulated_Surface(*triangulated_surface_list.triangulated_surface(index));
        // implicit surface
        index=implicit_surface_list.Add_Levelset_Implicit_Surface();
        sprintf(filename,"../../Public_Data/Rigid_Bodies/Rings_Test/cylinder_revolve.phi");input.open(filename,std::ios::in|std::ios::binary);
        implicit_surface_list.implicit_surface(index)->Read_Float(input);input.close();
        implicit_surface_list.implicit_surface(index)->Rescale(.5); // resize radius from the default of 1 to .25
        rigid_bodies(rigid_bodies.m)->Initialize_Implicit_Surface(*implicit_surface_list.implicit_surface(index));
        // rigid body
        sprintf(filename,"../../Public_Data/Rigid_Bodies/Rings_Test/cylinder_revolve.rgd");input.open(filename,std::ios::in|std::ios::binary);
        rigid_bodies(rigid_bodies.m)->Read_Float(input);input.close();
        rigid_bodies(rigid_bodies.m)->position=VECTOR_3D<T>(-(T).3,1.5,0); // reset position
        rigid_bodies(rigid_bodies.m)->orientation=roller_orientation;
        rigid_bodies(rigid_bodies.m)->angular_velocity=-roller_speed*VECTOR_3D<T>(0,0,1);
        rigid_bodies(rigid_bodies.m)->coefficient_of_friction=roller_friction;

        // cylinder 2
        rigid_bodies.Resize(rigid_bodies.m+1);rigid_bodies(rigid_bodies.m)=new RIGID_BODY<TV>;
        // triangulated surface
        index=triangulated_surface_list.Add_Triangulated_Surface();
        sprintf(filename,"../../Public_Data/Rigid_Bodies/Rings_Test/cylinder_revolve.tri");input.open(filename,std::ios::in|std::ios::binary);
        triangulated_surface_list.triangulated_surface(index)->Read_Float(input);input.close();
        triangulated_surface_list.triangulated_surface(index)->Rescale(.5); // resize radius from the default of 1 to .25
        rigid_bodies(rigid_bodies.m)->Initialize_Triangulated_Surface(*triangulated_surface_list.triangulated_surface(index));
        // implicit surface
        index=implicit_surface_list.Add_Levelset_Implicit_Surface();
        sprintf(filename,"../../Public_Data/Rigid_Bodies/Rings_Test/cylinder_revolve.phi");input.open(filename,std::ios::in|std::ios::binary);
        implicit_surface_list.implicit_surface(index)->Read_Float(input);input.close();
        implicit_surface_list.implicit_surface(index)->Rescale(.5); // resize radius from the default of 1 to .25
        rigid_bodies(rigid_bodies.m)->Initialize_Implicit_Surface(*implicit_surface_list.implicit_surface(index));
        // rigid body
        sprintf(filename,"../../Public_Data/Rigid_Bodies/Rings_Test/cylinder_revolve.rgd");input.open(filename,std::ios::in|std::ios::binary);
        rigid_bodies(rigid_bodies.m)->Read_Float(input);input.close();
        rigid_bodies(rigid_bodies.m)->position=VECTOR_3D<T>((T).3,1.5,0); // reset position
        rigid_bodies(rigid_bodies.m)->orientation=roller_orientation;
        rigid_bodies(rigid_bodies.m)->angular_velocity=roller_speed*VECTOR_3D<T>(0,0,1);
        rigid_bodies(rigid_bodies.m)->coefficient_of_friction=roller_friction;}
}
//#####################################################################
// Function Update_Rigid_Body_Positions
//#####################################################################
void Update_Rigid_Body_Positions(ARRAY<RIGID_BODY<TV>*>& rigid_bodies,const T time)
{
    rigid_bodies(2)->orientation=QUATERNION<T>(-roller_speed*time,VECTOR_3D<T>(0,0,1))*roller_orientation;
    rigid_bodies(3)->orientation=QUATERNION<T>(roller_speed*time,VECTOR_3D<T>(0,0,1))*roller_orientation;
}
//#####################################################################
};
}
#endif

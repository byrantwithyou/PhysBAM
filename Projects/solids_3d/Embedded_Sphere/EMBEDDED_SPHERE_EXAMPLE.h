//#####################################################################
// Copyright 2003-2005, Zhaosheng Bao, Geoffrey Irving, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_SPHERE_EXAMPLE
//#####################################################################
#ifndef __EMBEDDED_SPHERE_EXAMPLE__
#define __EMBEDDED_SPHERE_EXAMPLE__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include "../TET_SIM_EMBEDDED_COLLISIONS_EXAMPLE.h"
namespace PhysBAM{

template <class T>
class EMBEDDED_SPHERE_EXAMPLE:public TET_SIM_EMBEDDED_COLLISIONS_EXAMPLE<T>
{
public:
    SPHERE<T> sphere;
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    int m_input,n_input,mn_input;
    bool cube_mesh;

    EMBEDDED_SPHERE_EXAMPLE()
        :initial_height(1),initial_orientation(0,VECTOR_3D<T>(1,0,0)),initial_velocity(0,0,0),initial_angular_velocity(0,0,0),
        m_input(5),n_input(5),mn_input(5),cube_mesh(false)
    {
        final_time=5;
        restart_step_number=0;   
        cfl_number=(T)10;
        cg_tolerance=(T)1e-2;
        artificially_damp_strain_iterations=1;artificially_damp_strain_rate_iterations=10;
        artificial_min_strain=(T)-.5;artificial_max_strain=(T).5;
        use_masses_and_springs=true;use_altitude_springs=true;
        edge_spring_stiffness=3000;edge_spring_overdamping_fraction=2;
        artificially_damp_edge_spring_strain_rate=false;artificially_damp_edge_spring_strain=false;
        use_shortest_altitude_spring_only=true;
        use_altitude_springs_compressed_by_threshold_only=true;altitude_spring_fraction_compression=(T).1;
        altitude_spring_stiffness=300;altitude_spring_overdamping_fraction=2;
        artificially_damp_altitude_spring_strain_rate=false;artificially_damp_altitude_spring_strain=false;
        use_fvm=false;
        use_diagonalized_fvm=false;
        //use_linear_elasticity=true;
        youngs_modulus=100000;poissons_ratio=(T).45;Rayleigh_coefficient=(T).01;
        use_shortest_fvm_altitude_spring_only=true;
        use_fvm_altitude_springs_compressed_by_threshold_only=true;fvm_altitude_spring_fraction_compression=(T).1;
        artificially_damp_fvm_altitude_spring_strain_rate=false;artificially_damp_fvm_altitude_spring_strain=false;
        strcpy(output_directory,"Sphere/output");
        check_initial_mesh_for_self_intersection=false;
        perform_self_collision=false;
    }

    ~EMBEDDED_SPHERE_EXAMPLE()
    {}

//#####################################################################
// Function Remove_Tetrahedra_Completely_Outside_Level_Set
//#####################################################################
void Remove_Tetrahedra_Completely_Outside_Level_Set(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,const ARRAY<T>&phi)
{
    assert(phi.m == tetrahedralized_volume.particles.Size());
    int t=1;
    while(t <= tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m){
        int i,j,k,l;tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.Get(t,i,j,k,l);
        if(phi(i) > 0 && phi(j) > 0 && phi(k) > 0 && phi(l) > 0)tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.Remove_Index_Lazy(t);
        else t++;}
    tetrahedralized_volume.Discard_Valence_Zero_Particles_And_Renumber();
}
//#####################################################################
// Function Initialize_Embedded_Tetrahedralized_Volume
//#####################################################################
virtual void Initialize_Embedded_Tetrahedralized_Volume(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& embedded_tetrahedralized_volume)
{
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=embedded_tetrahedralized_volume.tetrahedralized_volume;
    tetrahedralized_volume.particles.Update_Position();
    GRID<TV> sphere_grid(m_input,n_input,mn_input,sphere.Bounding_Box());
    if(!cube_mesh) tetrahedralized_volume.Initialize_Octahedron_Mesh_And_Particles(sphere_grid);
    else tetrahedralized_volume.Initialize_Cube_Mesh_And_Particles(sphere_grid);
    tetrahedralized_volume.particles.Delete_Velocity_And_Acceleration();
    tetrahedralized_volume.particles.Delete_Mass(); // in case they're accidently stored in the .tet file
    tetrahedralized_volume.particles.Update_Position_And_Velocity();
    tetrahedralized_volume.particles.Store_Mass();

    ARRAY<T> phi(1,tetrahedralized_volume.particles.Size());
    for(int p=0;p<phi.m;p++) phi(p)=sphere.Signed_Distance(tetrahedralized_volume.particles.X(p));
    Remove_Tetrahedra_Completely_Outside_Level_Set(tetrahedralized_volume,phi);

    embedded_tetrahedralized_volume.embedded_surface.particles.Store_Position_And_Velocity();
    int hash_max=tetrahedralized_volume.particles.Size()*hash_ratio;
    std::cout << "Embedded Particles Hash Table Size: " << hash_max << std::endl;
    embedded_tetrahedralized_volume.Initialize_Parents_To_Embedded_Particles_Hash_Table(hash_max);
    embedded_tetrahedralized_volume.Initialize_Embedded_Triangles_In_Tetrahedron();
    embedded_tetrahedralized_volume.Initialize_Embedded_Children();    
    phi.Resize(1,tetrahedralized_volume.particles.Size());
    for(int p=0;p<phi.m;p++) phi(p)=sphere.Signed_Distance(tetrahedralized_volume.particles.X(p));
    assert(phi.Min() < 0);assert(phi.Max() > 0);
    embedded_tetrahedralized_volume.Calculate_Triangulated_Surface_From_Levelset_On_Tetrahedron_Nodes(phi);
    std::cout << "number of embedded_triangles=" << embedded_tetrahedralized_volume.embedded_surface.triangle_mesh.triangles.m << std::endl;
    std::cout << "total vertices = " << tetrahedralized_volume.particles.Size() << std::endl;
    std::cout << "total tets = " << tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m << std::endl; 
    
    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=0;i<tetrahedralized_volume.particles.array_size;i++){
        tetrahedralized_volume.particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i)=center+initial_orientation.Rotate(tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i).y+=initial_height-bottom;}
    embedded_tetrahedralized_volume.Update_Embedded_Particle_Positions();
    embedded_tetrahedralized_volume.Calculate_Mass_Of_Surface_Nodes();
    tetrahedralized_volume.Set_Density(500);
    tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);                        
}
//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
void Initialize_Rigid_Bodies(ARRAY<RIGID_BODY<TV>*,VECTOR<int,1> >& rigid_bodies)
{  
    std::fstream input;char filename[256];
    
    // plane
    rigid_bodies.Resize(1,rigid_bodies.m+1);rigid_bodies(rigid_bodies.m)=new RIGID_BODY<TV>;
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
    rigid_bodies(rigid_bodies.m)->coefficient_of_friction=(T).3;
}
//#####################################################################
};
}
#endif

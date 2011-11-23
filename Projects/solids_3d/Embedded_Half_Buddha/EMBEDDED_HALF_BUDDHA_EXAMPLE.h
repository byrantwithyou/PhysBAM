//#####################################################################
// Copyright 2003, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_HALF_BUDDHA_EXAMPLE
//#####################################################################
#ifndef __EMBEDDED_HALF_BUDDHA_EXAMPLE__
#define __EMBEDDED_HALF_BUDDHA_EXAMPLE__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include "../TET_SIM_EMBEDDED_COLLISIONS_EXAMPLE.h"
namespace PhysBAM{

template <class T>
class EMBEDDED_HALF_BUDDHA_EXAMPLE:public TET_SIM_EMBEDDED_COLLISIONS_EXAMPLE<T>
{
public:
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    int m_input,n_input,mn_input;
    bool cube_mesh;

    EMBEDDED_HALF_BUDDHA_EXAMPLE()
        :initial_height(1),initial_orientation(0,VECTOR_3D<T>(1,0,0)),initial_velocity(0,0,0),initial_angular_velocity(0,0,0),
        m_input(37),n_input(37),mn_input(13),cube_mesh(false)
    {
        final_time=10;
        restart_step_number=0;//100;   
        cfl_number=(T)10;
        cg_tolerance=(T)1e-2;
        use_masses_and_springs=true;use_altitude_springs=true;
        edge_spring_stiffness=.1*600;edge_spring_overdamping_fraction=2;
        use_shortest_altitude_spring_only=true;
        use_altitude_springs_compressed_by_threshold_only=true;altitude_spring_fraction_compression=(T).1;
        altitude_spring_stiffness=.1*60;altitude_spring_overdamping_fraction=2;
        use_diagonalized_fvm=false;use_linear_elasticity=false;use_neo_hookean=true;
        youngs_modulus=200000;poissons_ratio=(T).45;Rayleigh_coefficient=(T).01;
        strcpy(output_directory,"Half_Buddha/output");
        check_initial_mesh_for_self_intersection=false;
        perform_self_collision=true;
    }

    virtual ~EMBEDDED_HALF_BUDDHA_EXAMPLE()
    {}

//#####################################################################
// Function Remove_Tetrahedra_Completely_Outside_Level_Set
//#####################################################################
void Remove_Tetrahedra_Completely_Outside_Level_Set(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,const ARRAY<T>&phi)
{
    assert(phi.m == tetrahedralized_volume.particles.array_collection->Size());
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
    std::fstream input;char filename[256];
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=embedded_tetrahedralized_volume.tetrahedralized_volume;
    tetrahedralized_volume.particles.Update_Position();
    
    sprintf(filename,"../../Public_Data/Tetrahedralized_Volumes/half_buddha_in_five_cut_brick_low_res.tet");input.open(filename,std::ios::in|std::ios::binary);assert(input.is_open());
    tetrahedralized_volume.Read_Float(input);input.close();
    std::cout << "total vertices = " << tetrahedralized_volume.particles.array_collection->Size() << std::endl;
    std::cout << "total tetrhedra = " << tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m << std::endl; 
    tetrahedralized_volume.particles.Delete_Velocity_And_Acceleration();
    tetrahedralized_volume.particles.Delete_Mass(); // in case they're accidently stored in the .tet file
    tetrahedralized_volume.particles.Update_Position_And_Velocity();
    tetrahedralized_volume.particles.Store_Mass();

      GRID<TV> levelset_grid; // for levelset_implicit_surface
    ARRAY<T,VECTOR<int,3> > phi3d; // for levelset_implicit_surface
    LEVELSET_IMPLICIT_SURFACE<T> levelset_implicit_surface(levelset_grid,phi3d);
    sprintf(filename,"../../Public_Data/Rigid_Bodies/half_buddha.phi");input.open(filename,std::ios::in|std::ios::binary);
    if(!input.is_open()){std::cout << "TROUBLE OPENING " << filename << std::endl;return;}
    levelset_implicit_surface.Read(input);input.close();

    ARRAY<T> phi(1,tetrahedralized_volume.particles.array_collection->Size());
    for(int p=1;p<=phi.m;p++) phi(p)=levelset_implicit_surface(tetrahedralized_volume.particles.X(p));
    Remove_Tetrahedra_Completely_Outside_Level_Set(tetrahedralized_volume,phi);

    embedded_tetrahedralized_volume.embedded_surface.particles.Store_Position_And_Velocity();
    embedded_tetrahedralized_volume.Calculate_Mass_Of_Embedded_Nodes();
    int hash_max=tetrahedralized_volume.particles.array_collection->Size()*hash_ratio;
    std::cout << "Embedded Particles Hash Table Size: " << hash_max << std::endl;
    embedded_tetrahedralized_volume.Initialize_Parents_To_Embedded_Particles_Hash_Table(hash_max);
    embedded_tetrahedralized_volume.Initialize_Embedded_Sub_Elements_In_Parent_Element();
    embedded_tetrahedralized_volume.Initialize_Embedded_Children();
    embedded_tetrahedralized_volume.Set_Interpolation_Fraction_Threshold(interpolation_fraction_threshhold);
    phi.Resize(1,tetrahedralized_volume.particles.array_collection->Size());
    for(int p=1;p<=phi.m;p++) phi(p)=levelset_implicit_surface(tetrahedralized_volume.particles.X(p));
    assert(phi.Min() < 0);assert(phi.Max() > 0);
    std::cout << "about to place embedded triangles" << std::endl;
    embedded_tetrahedralized_volume.embedded_particles.Set_Array_Buffer_Size(1000);
    embedded_tetrahedralized_volume.embedded_triangle_mesh.Initialize_Incident_Triangles();
    embedded_tetrahedralized_volume.Calculate_Triangulated_Surface_From_Levelset_On_Tetrahedron_Nodes(phi);
    std::cout << "number of embedded_triangles = " << embedded_tetrahedralized_volume.embedded_surface.triangle_mesh.triangles.m << std::endl;
    std::cout << "total vertices = " << tetrahedralized_volume.particles.array_collection->Size() << std::endl;
    std::cout << "total tets = " << tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m << std::endl; 
    
    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=1;i<=tetrahedralized_volume.particles.array_size;i++){
        tetrahedralized_volume.particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i)=center+initial_orientation.Rotate(tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i).y+=initial_height-bottom;}
    embedded_tetrahedralized_volume.Update_Embedded_Particle_Positions();
    embedded_tetrahedralized_volume.Calculate_Mass_Of_Embedded_Nodes();
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

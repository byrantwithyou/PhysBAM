//#####################################################################
// Copyright 2003, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COARSE_SEMI_VOID_BRICK_EXAMPLE
//#####################################################################
#ifndef __COARSE_SEMI_VOID_BRICK_EXAMPLE__
#define __COARSE_SEMI_VOID_BRICK_EXAMPLE__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include "../TET_SIM_EMBEDDED_EXAMPLE.h"
#include <Geometry/LEVELSET_IMPLICIT_SURFACE.h>
using namespace PhysBAM;


template <class T>
class COARSE_SEMI_VOID_BRICK_EXAMPLE:public TET_SIM_EMBEDDED_EXAMPLE<T>
{
public:
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    int m_input,n_input,mn_input;
    T xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input; 
    GRID<TV> mattress_grid;
    bool cube_mesh;
    ARRAY<int> constrained_nodes;

    GRID<TV> levelset_grid; // for levelset_implicit_surface
    ARRAY<T,VECTOR<int,3> > phi3d; // for levelset_implicit_surface
    LEVELSET_IMPLICIT_SURFACE<T> levelset_implicit_surface;
    IMPLICIT_SURFACE<T>* implicit_surface;
    ARRAY<T> phi_on_tet_nodes;


    COARSE_SEMI_VOID_BRICK_EXAMPLE()    
                                :initial_height((T).25),initial_orientation(0,VECTOR_3D<T>(1,0,0)),
                                initial_velocity(0,0,0),initial_angular_velocity(0,0,0),
                                m_input(8),n_input(8),mn_input(8),           
                                xmin_input((T)-.5),xmax_input((T)-xmin_input),ymin_input((T)-.5),ymax_input((T)-ymin_input),zmin_input((T)-.5),zmax_input((T)-zmin_input),
                                mattress_grid(m_input,n_input,mn_input,xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input),
                                levelset_implicit_surface(levelset_grid,phi3d),
                                implicit_surface(0)
    
    
    {
       restart_step_number=0;
        final_time=5;                                         
        frame_rate=24;
        cfl_number=(T).4;
        cg_iterations=250;
        youngs_modulus=.01*(7.9e7);poissons_ratio=(T).3;Rayleigh_coefficient=(T).02;
        strcpy(output_directory,"Coarse_Semi_Void_Brick/output");
        
        perform_fracture=true;
        perform_self_collision=false;
        push_out=false;
        T scale=.1;
        fracture_threshold_1=(T)1e4*scale;fracture_threshold_2=(T)5e4*scale;fracture_threshold_3=(T)8e4*scale;
        compressive_threshold_1=(T)-100*1e4*scale;compressive_threshold_2=(T)-5e4*scale;compressive_threshold_3=(T)-8e4*scale;

        interpolation_fraction_threshhold=(T).0002;
        max_number_of_cuts=3;
        number_of_fracture_initiation_points=1000;

        bias_stress=true;

        solids_parameters.use_constant_mass=false;
    }

    virtual ~COARSE_SEMI_VOID_BRICK_EXAMPLE()
    {}



//#####################################################################
// Function Initialize_Tetrahedralized_Volume
//#####################################################################
void Initialize_Tetrahedralized_Volume(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume)
{
    std::cout << "we do nothing in initialize tetrahedralized volume for this example" << std::endl;
}
//#####################################################################
// Function Initialize_Embedded_Tetrahedralized_Volume
//#####################################################################
virtual void Initialize_Embedded_Tetrahedralized_Volume(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& embedded_tetrahedralized_volume)
{
    std::fstream input;char filename[256];
    embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Update_Position();
    
    embedded_tetrahedralized_volume.embedded_surface.particles.Store_Position();
    sprintf(filename,"../../Public_Data/Tetrahedralized_Volumes/coarse_brick_with_quarter_sphere_bite.etv");input.open(filename,std::ios::in|std::ios::binary);assert(input.is_open());
    embedded_tetrahedralized_volume.Read_Float(input);input.close();



    std::cout << "total vertices = " << embedded_tetrahedralized_volume.tetrahedralized_volume.particles.array_collection->Size() << std::endl;
    std::cout << "total tetrhedra = " << embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m << std::endl; 
    std::cout << "total embedded particles = " << embedded_tetrahedralized_volume.embedded_particles.array_collection->Size() << std::endl;
    std::cout << "parents array size = " << embedded_tetrahedralized_volume.parent_particles.m << std::endl;
    embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Delete_Velocity_And_Acceleration();
    embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Delete_Mass(); // in case they're accidently stored in the .tet file
    embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Update_Position_And_Velocity();
    embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Store_Mass();
    





    embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh.Initialize_Incident_Tetrahedrons();
    embedded_tetrahedralized_volume.tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(embedded_tetrahedralized_volume.tetrahedralized_volume.bounding_box->Center());T bottom=embedded_tetrahedralized_volume.tetrahedralized_volume.bounding_box->ymin;
    for(int i=1;i<=embedded_tetrahedralized_volume.tetrahedralized_volume.particles.array_size;i++){
        embedded_tetrahedralized_volume.tetrahedralized_volume.particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i)-center);
        embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i)=center+initial_orientation.Rotate(embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i)-center);
        embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i).y+=initial_height-bottom;}
    std::cout << "total vertices = " << embedded_tetrahedralized_volume.tetrahedralized_volume.particles.array_collection->Size() << std::endl;
    std::cout << "total tets = " << embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m << std::endl;
    embedded_tetrahedralized_volume.tetrahedralized_volume.Set_Density(1000);
    embedded_tetrahedralized_volume.tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);






    int hash_max=embedded_tetrahedralized_volume.tetrahedralized_volume.particles.array_collection->Size()*hash_ratio;
    std::cout << "Embedded Particles Hash Table Size: " << hash_max << std::endl;
    embedded_tetrahedralized_volume.Initialize_Parents_To_Embedded_Particles_Hash_Table(hash_max);
    embedded_tetrahedralized_volume.Initialize_Embedded_Sub_Elements_In_Parent_Element();
    embedded_tetrahedralized_volume.embedded_surface.particles.Set_Array_Buffer_Size(4*particle_buffer_size);
    embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Set_Array_Buffer_Size(particle_buffer_size);
}
//#####################################################################
// Function Initialize_Reference_Fracture_Tetrahedralized_Volume
//#####################################################################
virtual void Initialize_Reference_Fracture_Tetrahedralized_Volume(const EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& etv,FRACTURE_TETRAHEDRALIZED_VOLUME<T>& rftv,
                                                                  ARRAY<int> &reference_orientation_index)
{
    rftv.embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Delete_Velocity_And_Acceleration();
    rftv.embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Delete_Mass(); // in case they're accidently stored in the .tet file
    rftv.embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh.Initialize_Tetrahedron_Mesh(etv.tetrahedralized_volume.tetrahedron_mesh);
    rftv.Initialize();
    rftv.embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Initialize_Particles(etv.tetrahedralized_volume.particles);
    rftv.embedded_tetrahedralized_volume.embedded_particles.Initialize_Particles(etv.embedded_particles);
    rftv.embedded_tetrahedralized_volume.parent_particles.Resize(1,etv.parent_particles.m);
    ARRAY<int>::copy(etv.parent_particles,rftv.embedded_tetrahedralized_volume.parent_particles);
    rftv.embedded_tetrahedralized_volume.embedded_surface.triangle_mesh.Initialize_Triangle_Mesh(etv.embedded_surface.triangle_mesh);
    rftv.embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh.Initialize_Incident_Tetrahedrons();
    rftv.embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh.Initialize_Neighbor_Nodes();
    if(!rftv.embedded_tetrahedralized_volume.embedded_surface.triangle_mesh.number_nodes)rftv.embedded_tetrahedralized_volume.embedded_surface.triangle_mesh.incident_triangles=new ARRAY<ARRAY<int> >(1,0);
    else rftv.embedded_tetrahedralized_volume.embedded_surface.triangle_mesh.Initialize_Incident_Triangles();

    // Initialize hash table
    int hash_max=rftv.embedded_tetrahedralized_volume.tetrahedralized_volume.particles.array_collection->Size()*hash_ratio;
    std::cout << "Reference Embedded Particles Hash Table Size: " << hash_max << std::endl;
    rftv.embedded_tetrahedralized_volume.Initialize_Parents_To_Embedded_Particles_Hash_Table(hash_max);
    

    //new stuff for starting w/ .etv file w/ previous cuts
    rftv.embedded_tetrahedralized_volume.parent_particles.Resize(2,1,etv.parent_particles.m);
    ARRAY<int>::copy(etv.parent_particles,rftv.embedded_tetrahedralized_volume.parent_particles);
    rftv.embedded_tetrahedralized_volume.interpolation_fraction.Resize(1,etv.parent_particles.m);
    ARRAY<T>::copy(etv.interpolation_fraction,rftv.embedded_tetrahedralized_volume.interpolation_fraction);

    rftv.embedded_tetrahedralized_volume.Initialize_Embedded_Sub_Elements_In_Parent_Element();
    rftv.embedded_tetrahedralized_volume.Initialize_Embedded_Children();
    rftv.embedded_tetrahedralized_volume.Set_Interpolation_Fraction_Threshold(interpolation_fraction_threshhold);
    rftv.embedded_tetrahedralized_volume.embedded_surface.particles.Set_Array_Buffer_Size(particle_buffer_size);
    reference_orientation_index.Resize(1,rftv.embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m);
    rftv.embedded_tetrahedralized_volume.node_in_tetrahedron_is_material.Resize(1,etv.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m);
    ARRAY<bool>::copy(etv.node_in_tetrahedron_is_material,rftv.embedded_tetrahedralized_volume.node_in_tetrahedron_is_material);
    Initialize_Bias_Stress_Constants(etv,rftv);
    Calculate_Orientation_Index_Baed_On_Embedded_Triangles(rftv,reference_orientation_index);
}
//#####################################################################
// Function Calculate_Orientation_Index_Baed_On_Embedded_Triangles
//#####################################################################
void Calculate_Orientation_Index_Baed_On_Embedded_Triangles(FRACTURE_TETRAHEDRALIZED_VOLUME<T>& rftv,ARRAY<int> &reference_orientation_index)
{
    ARRAY<int>::copy(1,reference_orientation_index);
    reference_orientation_index.Resize(1,rftv.embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m);
    for(int t=1;t<=rftv.embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m;t++){
        int emb_tri1,emb_tri2,emb_tri3,emb_tri4,emb_tri_count;
        emb_tri_count=rftv.embedded_tetrahedralized_volume.Embedded_Triangles_In_Tetrahedron(t,emb_tri1,emb_tri2,emb_tri3,emb_tri4);
        int i,j,k,l;rftv.embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.Get(t,i,j,k,l);
        int oi=i,oj=j,ok=k,ol=l,count=1;
        if(emb_tri_count == 1){
            int a,b,c;rftv.embedded_tetrahedralized_volume.embedded_triangle_mesh.triangles.Get(emb_tri1,a,b,c);
            while(a != rftv.embedded_tetrahedralized_volume.Embedded_Particle_On_Segment(i,j) || b != rftv.embedded_tetrahedralized_volume.Embedded_Particle_On_Segment(i,k) || c != rftv.embedded_tetrahedralized_volume.Embedded_Particle_On_Segment(i,l)) 
                permute_four_ints(oi,oj,ok,ol,i,j,k,l,++count);
            reference_orientation_index(t)=count;}
        else if (emb_tri_count == 2){
            int a,b,c,d;rftv.embedded_tetrahedralized_volume.embedded_triangle_mesh.triangles.Get(emb_tri1,a,b,c);d=rftv.embedded_tetrahedralized_volume.embedded_triangle_mesh.triangles(2,emb_tri2);
            assert(rftv.embedded_tetrahedralized_volume.embedded_triangle_mesh.triangles(3,emb_tri2) == a);assert(rftv.embedded_tetrahedralized_volume.embedded_triangle_mesh.triangles(1,emb_tri2) == c);            
            while(a != rftv.embedded_tetrahedralized_volume.Embedded_Particle_On_Segment(i,k) || b != rftv.embedded_tetrahedralized_volume.Embedded_Particle_On_Segment(j,k) || c != rftv.embedded_tetrahedralized_volume.Embedded_Particle_On_Segment(j,l) || d != rftv.embedded_tetrahedralized_volume.Embedded_Particle_On_Segment(l,i)) 
                permute_four_ints(oi,oj,ok,ol,i,j,k,l,++count);
            reference_orientation_index(t)=count;}
        else assert(emb_tri_count == 0);}
}


//#####################################################################
// Function Set_External_Velocities
//#####################################################################
virtual void Set_External_Velocities(ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >& V,const T time)
{
    for(int i=1;i<=constrained_nodes.m;i++) V(constrained_nodes(i))=VECTOR_3D<T>(0,0,0);
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
virtual void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >& V,const T time)
{
    for(int i=1;i<=constrained_nodes.m;i++) V(constrained_nodes(i))=VECTOR_3D<T>(0,0,0);
} 
//#####################################################################
// Function Initialize_Bias_Stress_Constants
//#####################################################################
virtual void Initialize_Bias_Stress_Constants(const EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& etv,FRACTURE_TETRAHEDRALIZED_VOLUME<T>& rftv)
{
    TET_SIM_EMBEDDED_EXAMPLE<T>::Initialize_Bias_Stress_Constants(etv,rftv);

    etv.tetrahedralized_volume.Update_Bounding_Box();

}
//#####################################################################
};
#endif







//#####################################################################
// Copyright 2003, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BUDDHA_IN_SHAVED_BRICK_EXAMPLE
//#####################################################################
#ifndef __BUDDHA_IN_SHAVED_BRICK_EXAMPLE__
#define __BUDDHA_IN_SHAVED_BRICK_EXAMPLE__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include "../TET_SIM_EMBEDDED_EXAMPLE.h"
#include <Geometry/LEVELSET_IMPLICIT_SURFACE.h>
using namespace PhysBAM;


template <class T>
class BUDDHA_IN_SHAVED_BRICK_EXAMPLE:public TET_SIM_EMBEDDED_EXAMPLE<T>
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
    BOX_3D<T> box;


    BUDDHA_IN_SHAVED_BRICK_EXAMPLE()    
                                :initial_height((T)0),initial_orientation(0,VECTOR_3D<T>(1,0,0)),
                                initial_velocity(0,0,0),initial_angular_velocity(0,0,0),
                                m_input(8),n_input(8),mn_input(8),           
                                xmin_input((T)-.5),xmax_input((T)-xmin_input),ymin_input((T)-.5),ymax_input((T)-ymin_input),zmin_input((T)-.5),zmax_input((T)-zmin_input),
                                mattress_grid(m_input,n_input,mn_input,xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input),
                                levelset_implicit_surface(levelset_grid,phi3d),
                                implicit_surface(0)
    
    
    {
       restart_step_number=0;
        final_time=50;        
        frame_rate=30;
        cfl_number=(T)20;
        cg_iterations=250;
        youngs_modulus=1000;poissons_ratio=(T).4;Rayleigh_coefficient=(T).02;
        strcpy(output_directory,"F:/Molino_Data/buddha_brick_cut/flush_with_ground");
        
        perform_fracture=false;
        perform_self_collision=true;
        push_out=true;
        perturb_amount_for_collision_freeness=1e-5;
        gravity=1;
        T scale=.1;
        fracture_threshold_1=(T)1e4*scale;fracture_threshold_2=(T)5e4*scale;fracture_threshold_3=(T)8e4*scale;
        compressive_threshold_1=(T)-100*1e4*scale;compressive_threshold_2=(T)-5e4*scale;compressive_threshold_3=(T)-8e4*scale;
        interpolation_fraction_threshhold=(T).0000002;
        max_number_of_cuts=3;
        number_of_fracture_initiation_points=1000;

        bias_stress=true;

        solids_parameters.use_constant_mass=false;
    }

    virtual ~BUDDHA_IN_SHAVED_BRICK_EXAMPLE()
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
    sprintf(filename,"../../Public_Data/Tetrahedralized_Volumes/buddha_in_five_cut_brick_three_levels.tet");input.open(filename,std::ios::in|std::ios::binary);assert(input.is_open());
    embedded_tetrahedralized_volume.tetrahedralized_volume.Read_Float(input);input.close();

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
    for(int i=0;i<embedded_tetrahedralized_volume.tetrahedralized_volume.particles.array_size;i++){
        embedded_tetrahedralized_volume.tetrahedralized_volume.particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i)-center);
        embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i)=center+initial_orientation.Rotate(embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i)-center);
        embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i).y+=initial_height-bottom;}
    std::cout << "total vertices = " << embedded_tetrahedralized_volume.tetrahedralized_volume.particles.array_collection->Size() << std::endl;
    std::cout << "total tets = " << embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m << std::endl;
    embedded_tetrahedralized_volume.tetrahedralized_volume.Set_Density(100);
    embedded_tetrahedralized_volume.tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);


    int hash_max=embedded_tetrahedralized_volume.tetrahedralized_volume.particles.array_collection->Size()*hash_ratio;
    std::cout << "Embedded Particles Hash Table Size: " << hash_max << std::endl;
    embedded_tetrahedralized_volume.Initialize_Parents_To_Embedded_Particles_Hash_Table(hash_max);
    embedded_tetrahedralized_volume.Initialize_Embedded_Sub_Elements_In_Parent_Element();
    embedded_tetrahedralized_volume.embedded_surface.particles.Set_Array_Buffer_Size(4*particle_buffer_size);
    embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Set_Array_Buffer_Size(particle_buffer_size);
    embedded_tetrahedralized_volume.Set_Interpolation_Fraction_Threshold(interpolation_fraction_threshhold);
    embedded_tetrahedralized_volume.node_in_tetrahedron_is_material.Resize(4,1,embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m);
    Constrain_Nodes(embedded_tetrahedralized_volume.tetrahedralized_volume);
}
//#####################################################################
// Function Constrain_Nodes
//#####################################################################
virtual void Constrain_Nodes(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume)
{
    // Recompute Center
    tetrahedralized_volume.Update_Bounding_Box();
    
    VECTOR_3D<T> center=tetrahedralized_volume.bounding_box->Center();
    VECTOR_3D<T> size=tetrahedralized_volume.bounding_box->Size();
    VECTOR_3D<T> minimum_corner=tetrahedralized_volume.bounding_box->Minimum_Corner();
    std::cout << "Center: "<< center << std::endl;    
    std::cout << "Size: "<< size << std::endl;    
        
    box.Reset_Bounds(minimum_corner+VECTOR_3D<T>(1.1*size.x,1.1*size.y,0));
    box.ymin=center.y;
    box.xmin=minimum_corner.x + .5*size.x;
    box.zmax=minimum_corner.z + .5*size.z;
    box.zmin=minimum_corner.z - .001;

    BOX_3D<T> box2=box;
    box2.ymin=minimum_corner.y-.001;

    // Constrain Nodes
    for(int p=0;p<tetrahedralized_volume.particles.array_collection->Size();p++){
        if(box2.Lazy_Outside(tetrahedralized_volume.particles.X(p)))constrained_nodes.Append(p);
    }
    std::cout << "consstrained_nodes.m=" << constrained_nodes.m << std::endl;
}
//#####################################################################
// Function Add_Scripted_Cuts
//#####################################################################
virtual bool Add_Scripted_Cuts(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& reference_embedded_tetrahedralized_volume,
                               ARRAY<int>& reference_orientation_index,const T time)
{
    std::cout << "we're in BUDDHA_IN_SHAVED_BRICK::Add_Scripted_Cuts" << std::endl;
    int new_cuts_added=0;
    TETRAHEDRALIZED_VOLUME<T>& ref_tet_vol=reference_embedded_tetrahedralized_volume.tetrahedralized_volume;
    ref_tet_vol.Update_Bounding_Box();

    VECTOR_3D<T>normal(1,.15,-1);normal.Normalize();
    VECTOR_3D<T>point_on_cut_plane=box.Minimum_Corner()+VECTOR_3D<T>(.1*box.Size().x,box.Size().y,0);
    PLANE<T> plane(normal,point_on_cut_plane);



    T lambda=1;
    T height=LINEAR_INTERPOLATION<T,T>::Linear(box.ymax,((T)1.-lambda)*box.ymax+lambda*box.ymin,2*time/(5. - initial_time));

    int number_fully_in_box=0;
    for(int t=0;t<ref_tet_vol.tetrahedron_mesh.tetrahedrons.m;t++){
        int i,j,k,l;ref_tet_vol.tetrahedron_mesh.tetrahedrons.Get(t,i,j,k,l);
        VECTOR_3D<T> xi=ref_tet_vol.particles.X(i),xj=ref_tet_vol.particles.X(j),xk=ref_tet_vol.particles.X(k),xl=ref_tet_vol.particles.X(l);
        if(!ref_tet_vol.Completely_Inside_Box(t,box))continue;
        if(ref_tet_vol.Centroid(t).y < height)continue;
        number_fully_in_box++;
        int above_plane=(plane.Shortest_Signed_Distance(xi) > 0) + (plane.Shortest_Signed_Distance(xj) > 0) + (plane.Shortest_Signed_Distance(xk) > 0) + (plane.Shortest_Signed_Distance(xl) > 0);
        if(above_plane == 1){
            int oi=i,oj=j,ok=k,ol=l,count=1;
            while(!(plane.Shortest_Signed_Distance(ref_tet_vol.particles.X(i)) > 0))permute_four_ints(oi,oj,ok,ol,i,j,k,l,++count);
            reference_orientation_index(t)=count;
            T phi1=plane.Shortest_Signed_Distance(ref_tet_vol.particles.X(i)),phi2=plane.Shortest_Signed_Distance(ref_tet_vol.particles.X(j)),
            phi3=plane.Shortest_Signed_Distance(ref_tet_vol.particles.X(k)),phi4=plane.Shortest_Signed_Distance(ref_tet_vol.particles.X(l));
            if(!(phi1>0&&phi2<=0&&phi3<=0&&phi4<=0)){std::cout << "this is bad" << std::endl;exit(1);}
            int ri_rj=reference_embedded_tetrahedralized_volume.Add_Embedded_Particle_If_Not_Already_There(i,j,LEVELSET<T>::Theta(phi1,phi2));
            int ri_rk=reference_embedded_tetrahedralized_volume.Add_Embedded_Particle_If_Not_Already_There(i,k,LEVELSET<T>::Theta(phi1,phi3));
            int ri_rl=reference_embedded_tetrahedralized_volume.Add_Embedded_Particle_If_Not_Already_There(i,l,LEVELSET<T>::Theta(phi1,phi4));
            reference_embedded_tetrahedralized_volume.Add_Embedded_Triangle_If_Not_Already_There(ri_rj,ri_rk,ri_rl);
            new_cuts_added++;}    
        else if(above_plane == 2){
            int oi=i,oj=j,ok=k,ol=l,count=1;
            while(!(plane.Shortest_Signed_Distance(ref_tet_vol.particles.X(i)) > 0) || !(plane.Shortest_Signed_Distance(ref_tet_vol.particles.X(j)) > 0))
                permute_four_ints(oi,oj,ok,ol,i,j,k,l,++count);
            reference_orientation_index(t)=count;
            T phi1=plane.Shortest_Signed_Distance(ref_tet_vol.particles.X(i)),phi2=plane.Shortest_Signed_Distance(ref_tet_vol.particles.X(j)),
            phi3=plane.Shortest_Signed_Distance(ref_tet_vol.particles.X(k)),phi4=plane.Shortest_Signed_Distance(ref_tet_vol.particles.X(l));
            if(!(phi1>0&&phi2>0&&phi3<=0&&phi4<=0)){std::cout << "this is bad" << std::endl;exit(1);}
            int ri_rk=reference_embedded_tetrahedralized_volume.Add_Embedded_Particle_If_Not_Already_There(i,k,LEVELSET<T>::Theta(phi1,phi3));
            int ri_rl=reference_embedded_tetrahedralized_volume.Add_Embedded_Particle_If_Not_Already_There(i,l,LEVELSET<T>::Theta(phi1,phi4));
            int rj_rk=reference_embedded_tetrahedralized_volume.Add_Embedded_Particle_If_Not_Already_There(j,k,LEVELSET<T>::Theta(phi2,phi3));
            int rj_rl=reference_embedded_tetrahedralized_volume.Add_Embedded_Particle_If_Not_Already_There(j,l,LEVELSET<T>::Theta(phi2,phi4));
            reference_embedded_tetrahedralized_volume.Add_Embedded_Triangle_If_Not_Already_There(ri_rk,rj_rk,rj_rl);
            reference_embedded_tetrahedralized_volume.Add_Embedded_Triangle_If_Not_Already_There(rj_rl,ri_rl,ri_rk);
            new_cuts_added++;}
        else if(above_plane == 3){
            int oi=i,oj=j,ok=k,ol=l,count=1;
            while(!(plane.Shortest_Signed_Distance(ref_tet_vol.particles.X(i)) <= 0))permute_four_ints(oi,oj,ok,ol,i,j,k,l,++count);
            reference_orientation_index(t)=count;
            T phi1=plane.Shortest_Signed_Distance(ref_tet_vol.particles.X(i)),phi2=plane.Shortest_Signed_Distance(ref_tet_vol.particles.X(j)),
            phi3=plane.Shortest_Signed_Distance(ref_tet_vol.particles.X(k)),phi4=plane.Shortest_Signed_Distance(ref_tet_vol.particles.X(l));
            if(!(phi1<=0&&phi2>0&&phi3>0&&phi4>0)){std::cout << "this is bad" << std::endl;exit(1);}
            int ri_rj=reference_embedded_tetrahedralized_volume.Add_Embedded_Particle_If_Not_Already_There(i,j,LEVELSET<T>::Theta(phi1,phi2));
            int ri_rk=reference_embedded_tetrahedralized_volume.Add_Embedded_Particle_If_Not_Already_There(i,k,LEVELSET<T>::Theta(phi1,phi3));
            int ri_rl=reference_embedded_tetrahedralized_volume.Add_Embedded_Particle_If_Not_Already_There(i,l,LEVELSET<T>::Theta(phi1,phi4));
            reference_embedded_tetrahedralized_volume.Add_Embedded_Triangle_If_Not_Already_There(ri_rj,ri_rk,ri_rl);
            new_cuts_added++;} 
        else {if(!(above_plane == 0 || above_plane == 4)){std::cout << "this is bad" << std::endl;exit(1);}}}
    std::cout << "number_fully_in_box=" << number_fully_in_box << std::endl;
    return (new_cuts_added > 0);
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
    rftv.embedded_tetrahedralized_volume.node_in_tetrahedron_is_material.Resize(4,1,etv.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m);
    ARRAY<bool>::copy(true,rftv.embedded_tetrahedralized_volume.node_in_tetrahedron_is_material);
    Initialize_Bias_Stress_Constants(etv,rftv);
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
virtual void Set_External_Velocities(ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >& V,const T time)
{
    for(int i=0;i<constrained_nodes.m;i++) V(constrained_nodes(i))=VECTOR_3D<T>(0,0,0);
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
virtual void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >& V,const T time)
{
    for(int i=0;i<constrained_nodes.m;i++) V(constrained_nodes(i))=VECTOR_3D<T>(0,0,0);
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







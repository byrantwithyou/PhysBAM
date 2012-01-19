//#####################################################################
// Copyright 2003, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COARSE_SPHERE_EXAMPLE
//#####################################################################
//
//#####################################################################
#ifndef __COARSE_SPHERE_EXAMPLE__
#define __COARSE_SPHERE_EXAMPLE__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include "TET_SIM_EMBEDDED_EXAMPLE.h"
using namespace PhysBAM;

template <class T>
class COARSE_SPHERE_EXAMPLE:public TET_SIM_EMBEDDED_EXAMPLE<T>
{
public:
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    int m_input,n_input,mn_input;
    T xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input; 
    GRID<TV> mattress_grid;
    bool cube_mesh;
    int constraint_type;
    bool constrain_nodes;


    COARSE_SPHERE_EXAMPLE()
                                :initial_height(1),initial_orientation(3.1415/2.0,VECTOR_3D<T>(1,0,0)),
                                initial_velocity(0,0,0),initial_angular_velocity(4,0,0),
                                m_input(4),n_input(4),mn_input(4),
                                ymin_input(-.25),ymax_input(.25),xmin_input(-.25),xmax_input(.25),zmin_input(-.25),zmax_input(.25),
                                mattress_grid(m_input,n_input,mn_input,xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input),
                                cube_mesh(true),constraint_type(1), constrain_nodes(false)
    {
        restart_step_number=0;
        final_time=5;
        frame_rate=100;
        cfl_number=.5;
        cg_iterations=250;
        use_fvm=true;youngs_modulus=900000;poissons_ratio=.4;Rayleigh_coefficient=.005;
        strcpy(output_directory,"Coarse_Sphere/output");
        gravity=10;
        use_linear_elasticity=false;
        substeps_per_fracture_check=1;
        perform_fracture=false;
        estimate_fracture_normal=false;
        use_cauchy_stress=true;
        facture_ratio_before_rebuild=1.;
        limit_time_step_by_strain_rate=true;
        artificially_damp_strain_iterations=0;
        artificially_damp_strain_rate_iterations=0;
        fracture_threshold=2.0e4;
        fracture_bias_constant=1.0e4;
        interpolation_fraction_threshhold_ratio=0.2;
        number_of_cuts=4;
        number_of_fracture_initiation_points=1;
        perform_self_collision=false;
    }

    virtual ~COARSE_SPHERE_EXAMPLE()
    {}


//#####################################################################
// Function Initialize_Tetrahedralized_Volume
//#####################################################################
void Initialize_Tetrahedralized_Volume(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume)
{
    tetrahedralized_volume.particles.Update_Position();
    if(!cube_mesh) tetrahedralized_volume.Initialize_Octahedron_Mesh_And_Particles(mattress_grid);
    else tetrahedralized_volume.Initialize_Cube_Mesh_And_Particles(mattress_grid);

    tetrahedralized_volume.particles.Delete_Velocity_And_Acceleration();
    tetrahedralized_volume.particles.Delete_Mass(); // in case they're accidently stored in the .tet file
    tetrahedralized_volume.particles.Update_Position_And_Velocity();
    tetrahedralized_volume.particles.Store_Mass();

    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=0;i<tetrahedralized_volume.particles.array_size;i++){
        tetrahedralized_volume.particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i)=center+initial_orientation.Rotate(tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i).y+=initial_height-bottom;}
    std::cout << "total vertices = " << tetrahedralized_volume.particles.array_collection->Size() << std::endl;
    std::cout << "total tets = " << tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m << std::endl; 
    tetrahedralized_volume.Set_Density(500);
    tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);
}
//#####################################################################
// Function Initialize_Tetrahedralized_Volume
//#####################################################################

virtual void Initialize_Embedded_Tetrahedralized_Volume(const TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,
                                                        EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& embedded_tetrahedralized_volume)
{
    TET_SIM_EMBEDDED_EXAMPLE<T>::Initialize_Embedded_Tetrahedralized_Volume(tetrahedralized_volume,embedded_tetrahedralized_volume);    
    /*
    * above call does the following
    *
    embedded_tetrahedralized_volume.embedded_surface.particles.Store_Position();
    int hash_max=tetrahedralized_volume.particles.array_collection->Size()*hash_ratio;
    std::cout << "Embedded Particles Hash Table Size: " << hash_max << std::endl;
    embedded_tetrahedralized_volume.Initialize_Parents_To_Embedded_Particles_Hash_Table(hash_max);
    embedded_tetrahedralized_volume.Initialize_Embedded_Triangles_In_Tetrahedron();
    embedded_tetrahedralized_volume.embedded_surface.particles.Set_Array_Buffer_Size(4*particle_buffer_size);
    embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Set_Array_Buffer_Size(particle_buffer_size);    
    */

    // read levelset implicit surface from file
    embedded_tetrahedralized_volume.tetrahedralized_volume.Update_Bounding_Box();
    T radius=(T).8*((tetrahedralized_volume.bounding_box->Size()).Min());
    SPHERE<T> s(tetrahedralized_volume.bounding_box->Center(),radius);
    ARRAY<T> phi(1,tetrahedralized_volume.particles.array_collection->Size());
    for(int p=0;p<phi.m;p++) {
        VECTOR_3D<T> xp=tetrahedralized_volume.particles.X(p);
        std::cout << xp << std::endl;
        phi(p)=((xp - s.center).Magnitude_Squared() - sqr(s.radius));
        if(phi(p) > 0) phi(p)=sqrt(phi(p));
        else phi(p)= -sqrt(-phi(p));
    }

    for(int p=0;p<phi.m;p++) {
        std::cout << phi(p) << std::endl;
    }
    embedded_tetrahedralized_volume.Initialize_Embedded_Triangles_In_Tetrahedron();
    embedded_tetrahedralized_volume.Calculate_Triangulated_Surface_From_Levelset_On_Tetrahedron_Nodes(phi);


                        


    /*    GRID_3D<double> levelset_grid; // for levelset_implicit_surface
    ARRAY<double,VECTOR<int,3> > phi3d; // for levelset_implicit_surface
    LEVELSET_IMPLICIT_SURFACE levelset_implicit_surface(levelset_grid,phi3d);
    IMPLICIT_SURFACE* implicit_surface;    
    char implicit_surface_filename[256];sprintf(implicit_surface_filename,"../../Public_Data/Rigid_Bodies/Sphere_100x100x100.phi"); 
    std::fstream input;input.open(implicit_surface_filename,std::ios::in|std::ios::binary);
    if(!input.is_open()){std::cout << "TROUBLE OPENING" << implicit_surface_filename << std::endl;return;}
    levelset_implicit_surface.Read(input);input.close();implicit_surface=&levelset_implicit_surface;
    implicit_surface->Update_Box();VECTOR_3D<double> size=implicit_surface->box.Size();
    double cell_size=.24*min(size.x,size.y,size.z);
    int m=(int)ceil(size.x/cell_size),n=(int)ceil(size.y/cell_size),mn=(int)ceil(size.z/cell_size);
    GRID_3D<double> meshing_grid(m,n,mn,implicit_surface->box.xmin,implicit_surface->box.xmin+cell_size*m,implicit_surface->box.ymin,implicit_surface->box.ymin+cell_size*n,
                                implicit_surface->box.zmin,implicit_surface->box.zmin+cell_size*mn);
    embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Delete_Velocity_And_Acceleration();// in case they're accidently stored in the .etv file
    embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Delete_Mass(); // in case they're accidently stored in the .etv file
    embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Store_Position();
    embedded_tetrahedralized_volume.tetrahedralized_volume.Initialize_Cube_Mesh_And_Particles(meshing_grid);    
    embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Update_Position_And_Velocity();
    embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Store_Mass();
    embedded_tetrahedralized_volume.Calculate_Phi_On_Tetrahedron_Vertices(*implicit_surface);    
    Remove_Tetrahedra_Outside_Levelset(embedded_tetrahedralized_volume,true); //this must be done after phi on vertices, but before calculating triangulated surface
    embedded_tetrahedralized_volume.Calculate_Phi_On_Tetrahedron_Vertices(*implicit_surface);    
    embedded_tetrahedralized_volume.Calculate_Triangulated_Surface();
    Write_Tetrahedralized_Volume(embedded_tetrahedralized_volume,"initial_mesh.tet");
    Write_Triangulated_Surface(embedded_tetrahedralized_volume,"initial_tri_surf.tri");
    Write_Embedded_Tetrahedralized_Volume(embedded_tetrahedralized_volume,"initial_tet_vol.etv");

    //shift above floor && scale
    int i;
    for(i=0;i<embedded_tetrahedralized_volume.tetrahedralized_volume.particles.array_collection->Size();i++) {embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i).y+=3;embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i)*=.2;}
    embedded_tetrahedralized_volume.tetrahedralized_volume.Update_Bounding_Box();
    
    if(restart_step_number == 0){
        VECTOR_3D<double> center(embedded_tetrahedralized_volume.tetrahedralized_volume.bounding_box->Center());double bottom=embedded_tetrahedralized_volume.tetrahedralized_volume.bounding_box->ymin;
        for(i=0;i<embedded_tetrahedralized_volume.tetrahedralized_volume.particles.array_size;i++){
            embedded_tetrahedralized_volume.tetrahedralized_volume.particles.V(i)=initial_velocity+VECTOR_3D<double>::Cross_Product(initial_angular_velocity,embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i)-center);
            embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i)=center+initial_orientation.Rotate(embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i)-center);
            embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i).y+=initial_height-bottom;}}
    Get_Mass_Of_Each_Node(embedded_tetrahedralized_volume.tetrahedralized_volume);
    if(verbose) {std::cout << "embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m=" << embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m << std::endl;}
    */
}
//#####################################################################
};
#endif

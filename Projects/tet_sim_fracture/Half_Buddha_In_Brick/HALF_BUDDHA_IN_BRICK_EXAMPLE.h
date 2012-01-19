//#####################################################################
// Copyright 2003, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HALF_BUDDHA_IN_BRICK_EXAMPLE
//#####################################################################
#ifndef __HALF_BUDDHA_IN_BRICK_EXAMPLE__
#define __HALF_BUDDHA_IN_BRICK_EXAMPLE__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/GRAIN_BOUNDARIES.h>
#include "../TET_SIM_EMBEDDED_EXAMPLE.h"
#include <Geometry/LEVELSET_IMPLICIT_SURFACE.h>
using namespace PhysBAM;


template <class T>
class HALF_BUDDHA_IN_BRICK_EXAMPLE:public TET_SIM_EMBEDDED_EXAMPLE<T>
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
    GRAIN_BOUNDARIES<T>* grain_boundaries;


    HALF_BUDDHA_IN_BRICK_EXAMPLE()    
                                :initial_height((T).25),initial_orientation(2.1,VECTOR_3D<T>(1,0,0)),
                                initial_velocity(0,0,0),initial_angular_velocity(0,0,0),
                                m_input(8),n_input(8),mn_input(8),           
                                xmin_input((T)-.5),xmax_input((T)-xmin_input),ymin_input((T)-.5),ymax_input((T)-ymin_input),zmin_input((T)-.5),zmax_input((T)-zmin_input),
                                mattress_grid(m_input,n_input,mn_input,xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input),
                                levelset_implicit_surface(levelset_grid,phi3d),
                                implicit_surface(0),grain_boundaries(0)
    
    
    {
       restart_step_number=0;
        final_time=5;                                         
        frame_rate=24;
        cfl_number=(T).4;
        cg_iterations=250;
        youngs_modulus=.01*(7.9e7);poissons_ratio=(T).35;Rayleigh_coefficient=5*(T).03;
        strcpy(output_directory,"Half_Buddha_In_Brick/output");
        
        perform_fracture=true;
        perform_self_collision=false;
        push_out=false;
        T scale=10;
        fracture_threshold_1=(T)1e4*scale;fracture_threshold_2=(T)5e4*scale;fracture_threshold_3=(T)8e4*scale;
        compressive_threshold_1=(T)-1e4*scale;compressive_threshold_2=(T)-5e4*scale;compressive_threshold_3=(T)-8e4*scale;

        interpolation_fraction_threshhold=(T).0002;
        max_number_of_cuts=3;
        number_of_fracture_initiation_points=1000;

        bias_stress=true;

        solids_parameters.use_constant_mass=false;
        Initialize_Levelset();
    }

    virtual ~HALF_BUDDHA_IN_BRICK_EXAMPLE()
    {}


    void Initialize_Levelset()
    {
        char implicit_surface_filename[256];sprintf(implicit_surface_filename,"../../Public_Data/Rigid_Bodies/half_buddha.phi"); 
        std::fstream input;input.open(implicit_surface_filename,std::ios::in|std::ios::binary);
        if(!input.is_open()){std::cout << "TROUBLE OPENING" << implicit_surface_filename << std::endl;return;}
        levelset_implicit_surface.Read_Float(input);input.close();implicit_surface=&levelset_implicit_surface;
        implicit_surface->Update_Box();
    }


//#####################################################################
// Function Initialize_Tetrahedralized_Volume
//#####################################################################
void Initialize_Tetrahedralized_Volume(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume)
{

    std::fstream input;char filename[256];
    tetrahedralized_volume.particles.Update_Position();
    
    sprintf(filename,"../../Public_Data/Tetrahedralized_Volumes/half_buddha_in_five_cut_brick_low_res.tet");input.open(filename,std::ios::in|std::ios::binary);assert(input.is_open());
    tetrahedralized_volume.Read_Float(input);input.close();
    std::cout << "total vertices = " << tetrahedralized_volume.particles.array_collection->Size() << std::endl;
    std::cout << "total tetrhedra = " << tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m << std::endl; 
    tetrahedralized_volume.particles.Delete_Velocity_And_Acceleration();
    tetrahedralized_volume.particles.Delete_Mass(); // in case they're accidently stored in the .tet file
    tetrahedralized_volume.particles.Update_Position_And_Velocity();
    tetrahedralized_volume.particles.Store_Mass();

    tetrahedralized_volume.tetrahedron_mesh.Initialize_Incident_Tetrahedrons();

    phi_on_tet_nodes.Resize(1,tetrahedralized_volume.particles.array_collection->Size());
    for(int t=1;t<=tetrahedralized_volume.particles.array_collection->Size();t++)phi_on_tet_nodes(t)=(*implicit_surface)(tetrahedralized_volume.particles.X(t));

    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=0;i<tetrahedralized_volume.particles.array_size;i++){
        tetrahedralized_volume.particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i)=center+initial_orientation.Rotate(tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i).y+=initial_height-bottom;}
    std::cout << "total vertices = " << tetrahedralized_volume.particles.array_collection->Size() << std::endl;
    std::cout << "total tets = " << tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m << std::endl;
    tetrahedralized_volume.Set_Density(1000);
    tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);

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



    rftv.Set_Fracture_Bias_Stress_Scaling(1);
    rftv.Set_Fracture_Bias_Magnitude(.1*youngs_modulus);
    rftv.Set_Fracture_Bias_Propagation(0);

    rftv.fracture_bias_direction_coefficient=.5;
    rftv.eigenvector_coefficient=.5;
    rftv.fracture_bias_propagation_coefficient=0;
    
    delete grain_boundaries;
    grain_boundaries=new GRAIN_BOUNDARIES<T>(rftv,rftv.embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh);
    int number_of_regions=10;
    grain_boundaries->Fill_Regions_With_Uniform_Vectors(number_of_regions);
    grain_boundaries->Smooth_Fracture_Bias_Directions(19);


    etv.tetrahedralized_volume.Update_Bounding_Box();    
    // set phi
    for(int t=0;t<etv.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m;t++){
        int i,j,k,l;etv.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.Get(t,i,j,k,l);
        T phi1=phi_on_tet_nodes(i);T phi2=phi_on_tet_nodes(j);T phi3=phi_on_tet_nodes(k);T phi4=phi_on_tet_nodes(l);
        int positive_count=(phi1 > 0)+(phi2 > 0)+(phi3 > 0)+(phi4 > 0);
        if(positive_count == 1 || positive_count == 2 || positive_count ==3){
            rftv.Set_Phi(t,phi1,phi2,phi3,phi4);
            rftv.fracture_bias_stress_scaling(t)=4e4;            
            //rftv.fracture_bias_direction(t)=sphere.Grad_Phi(etv.tetrahedralized_volume.Centroid(t));
            //rftv.fracture_bias_magnitude(t)=4e4;
        }
        //if(positive_count == 4)rftv.fracture_bias_stress_scaling(t)=50;
        //if(positive_count == 0)rftv.fracture_bias_stress_scaling(t)=0;
    }


}
////#####################################################################
//// Function Initialize_Tetrahedralized_Volume
////#####################################################################
//void Initialize_Embedded_Tetrahedralized_Volume(EMBEDDED_TETRAHEDRALIZED_VOLUME& embedded_tetrahedralized_volume,bool verbose=true)
//{
//
//    if(restart_step_number == 0) {  
//        // read levelset implicit surface from file
//        embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Store_Position();
//        embedded_tetrahedralized_volume.tetrahedralized_volume.Initialize_Cube_Mesh_And_Particles(meshing_grid);    
//        embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Update_Position_And_Velocity();
//        embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Store_Mass();
//        embedded_tetrahedralized_volume.Calculate_Phi_On_Tetrahedron_Vertices(*implicit_surface);   
//        Remove_Tetrahedra_Outside_Levelset(embedded_tetrahedralized_volume,true); //this must be done after phi on vertices, but before calculating triangulated surface
//        embedded_tetrahedralized_volume.Calculate_Phi_On_Tetrahedron_Vertices(*implicit_surface);   
//        embedded_tetrahedralized_volume.Calculate_Triangulated_Surface();
//        Write_Tetrahedralized_Volume(embedded_tetrahedralized_volume,"initial_mesh.tet");
//        Write_Triangulated_Surface(embedded_tetrahedralized_volume,"initial_tri_surf.tri");
//        Write_Embedded_Tetrahedralized_Volume(embedded_tetrahedralized_volume,"initial_tet_vol.etv");}
//    else {
//        std::fstream input;char embedded_tetrahedralized_volume_filename[256];
//        sprintf(embedded_tetrahedralized_volume_filename,"%s/initial_tet_vol.etv",output_directory);input.open(embedded_tetrahedralized_volume_filename,std::ios::in|std::ios::binary);
//        if(!input.is_open()){std::cout << "TROUBLE OPENING" << embedded_tetrahedralized_volume_filename << std::endl;return;}
//        embedded_tetrahedralized_volume.Read(input);input.close();
//        embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Delete_Velocity_And_Acceleration();// in case they're accidently stored in the .etv file
//        embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Delete_Mass(); // in case they're accidently stored in the .etv file
//        embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Update_Position_And_Velocity();
//        embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Store_Mass();
//    }
//
//    //shift above floor && scale
//    int i;
//    for(i=1;i<=embedded_tetrahedralized_volume.tetrahedralized_volume.particles.array_collection->Size();i++) {embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i).y+=3;embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i)*=.2;}
//    embedded_tetrahedralized_volume.tetrahedralized_volume.Update_Bounding_Box();
//    
//    if(restart_step_number == 0){
//        VECTOR_3D<T> center(embedded_tetrahedralized_volume.tetrahedralized_volume.bounding_box->Center());T bottom=embedded_tetrahedralized_volume.tetrahedralized_volume.bounding_box->ymin;
//        for(i=1;i<=embedded_tetrahedralized_volume.tetrahedralized_volume.particles.array_size;i++){
//            embedded_tetrahedralized_volume.tetrahedralized_volume.particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i)-center);
//            embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i)=center+initial_orientation.Rotate(embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i)-center);
//            embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i).y+=initial_height-bottom;}}
//    Get_Mass_Of_Each_Node(embedded_tetrahedralized_volume.tetrahedralized_volume);
//    if(verbose) {std::cout << "embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m=" << embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m << std::endl;}
//}
//#####################################################################
};
#endif







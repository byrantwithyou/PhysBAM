//#####################################################################
// Copyright 2003, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DUMBBELL_EXAMPLE
//#####################################################################
#ifndef __DUMBBELL_EXAMPLE__
#define __DUMBBELL_EXAMPLE__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include "../TET_SIM_EMBEDDED_EXAMPLE.h"
#include <Geometry/LEVELSET_IMPLICIT_SURFACE.h>
using namespace PhysBAM;


template <class T>
class DUMBBELL_EXAMPLE:public TET_SIM_EMBEDDED_EXAMPLE<T>
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


    DUMBBELL_EXAMPLE()    
                                :initial_height((T)1),initial_orientation(0,VECTOR_3D<T>(1,0,0)),
                                initial_velocity(0,0,0),initial_angular_velocity(0,0,0),
                                m_input(8),n_input(8),mn_input(8),           
                                xmin_input((T)-.5),xmax_input((T)-xmin_input),ymin_input((T)-.5),ymax_input((T)-ymin_input),zmin_input((T)-.5),zmax_input((T)-zmin_input),
                                mattress_grid(m_input,n_input,mn_input,xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input),
                                cube_mesh(true),
                                levelset_implicit_surface(levelset_grid,phi3d),
                                implicit_surface(0)
    
    
    {
       restart_step_number=0;
        final_time=5;                                         
        frame_rate=24;
        cfl_number=(T).5;
        cg_iterations=250;
        youngs_modulus=.01*(7.9e7);poissons_ratio=(T).4;Rayleigh_coefficient=5*(T).02;
        strcpy(output_directory,"Dumbbell/output");
        
        perform_fracture=true;
        perform_self_collision=false;
        T scale=10;
        fracture_threshold_1=(T)1e4*scale;fracture_threshold_2=(T)5e4*scale;fracture_threshold_3=(T)8e4*scale;
        compressive_threshold_1=(T)-1e4*scale;compressive_threshold_2=(T)-5e4*scale;compressive_threshold_3=(T)-8e4*scale;

        interpolation_fraction_threshhold=(T).0002;
        max_number_of_cuts=3;
        number_of_fracture_initiation_points=1000;

        bias_stress=true;

        Initialize_Levelset();
    }

    virtual ~DUMBBELL_EXAMPLE()
    {}


    void Initialize_Levelset()
    {
        char implicit_surface_filename[256];sprintf(implicit_surface_filename,"../../Public_Data/Rigid_Bodies/dumbbell.phi"); 
        std::fstream input;input.open(implicit_surface_filename,std::ios::in|std::ios::binary);
        if(!input.is_open()){std::cout << "TROUBLE OPENING" << implicit_surface_filename << std::endl;return;}
        levelset_implicit_surface.Read(input);input.close();implicit_surface=&levelset_implicit_surface;
        implicit_surface->Update_Box();VECTOR_3D<T> size=implicit_surface->box.Size();
        T cell_size=(T).1*min(size.x,size.y,size.z);
        int m=(int)ceil(size.x/cell_size),n=(int)ceil(size.y/cell_size),mn=(int)ceil(size.z/cell_size);
        std::cout << "box size=" << size << std::endl;
        std::cout << "cell_size=" << cell_size << std::endl;
        std::cout << "(m,n,mn)=" << "(" << m << "," << n << "," << mn << ")" << std::endl;        
        mattress_grid.Initialize(m,n,mn,implicit_surface->box);
    }


//#####################################################################
// Function Initialize_Tetrahedralized_Volume
//#####################################################################
void Initialize_Tetrahedralized_Volume(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume)
{
    tetrahedralized_volume.particles.Update_Position();
    if(!cube_mesh) tetrahedralized_volume.Initialize_Octahedron_Mesh_And_Particles(mattress_grid);
    else tetrahedralized_volume.Initialize_Cube_Mesh_And_Particles(mattress_grid);

    for(int i=0;i<mattress_grid.m;i++)for(int k=0;k<mattress_grid.mn;k++){
        constrained_nodes.Append(tetrahedralized_volume.particles.array_collection->Size()+1-(1+i+k*mattress_grid.m*mattress_grid.n));}


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

    etv.tetrahedralized_volume.Update_Bounding_Box();
    
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
        if(positive_count == 4)rftv.fracture_bias_stress_scaling(t)=50;
        if(positive_count == 0)rftv.fracture_bias_stress_scaling(t)=0;
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







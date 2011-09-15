//#####################################################################
// Copyright 2003, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPHERE_EXAMPLE
//#####################################################################
//
//#####################################################################
// Molino - April 26, 2003
//#####################################################################
#ifndef __SPHERE_EXAMPLE__
#define __SPHERE_EXAMPLE__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Dynamics/Meshing/TETRAHEDRAL_MESHING.h>
#include "TET_SIM_EMBEDDED_EXAMPLE.h"
#include <Geometry/LEVELSET_IMPLICIT_SURFACE.h>
using namespace PhysBAM;

class SPHERE_EXAMPLE:public TET_SIM_EMBEDDED_EXAMPLE 
{
public:
    double initial_height;
    QUATERNION initial_orientation;
    VECTOR_3D<double> initial_velocity,initial_angular_velocity;
    GRID_3D<double> meshing_grid;
    
    SPHERE_EXAMPLE()
                                :initial_height(1),initial_orientation(0/*pi/2*/,VECTOR_3D<double>(1,0,0)),
                                initial_velocity(0,0,0),initial_angular_velocity(0,0,0)
    {
        restart_step_number=0;
        final_time=5;
        use_masses_and_springs=false;
        use_fvm=true;
        fracture_based_on_stress=false;
        youngs_modulus=900000;poissons_ratio=.4;Rayleigh_coefficient=.01;
        strcpy(output_directory,"Sphere/output");
    }

    ~SPHERE_EXAMPLE()
    {}

//#####################################################################
// Function Initialize_Tetrahedralized_Volume
//#####################################################################
void Initialize_Embedded_Tetrahedralized_Volume(EMBEDDED_TETRAHEDRALIZED_VOLUME& embedded_tetrahedralized_volume,bool verbose=true)
{    
    // read levelset implicit surface from file
    GRID_3D<double> levelset_grid; // for levelset_implicit_surface
    ARRAY<double,VECTOR<int,3> > phi3d; // for levelset_implicit_surface
    LEVELSET_IMPLICIT_SURFACE levelset_implicit_surface(levelset_grid,phi3d);
    IMPLICIT_SURFACE* implicit_surface;    
    char implicit_surface_filename[256];sprintf(implicit_surface_filename,"../../Public_Data/Rigid_Bodies/Sphere_100x100x100.phi"); 
    std::fstream input;input.open(implicit_surface_filename,std::ios::in|std::ios::binary);
    if(!input.is_open()){std::cout << "TROUBLE OPENING" << implicit_surface_filename << std::endl;return;}
    levelset_implicit_surface.Read(input);input.close();implicit_surface=&levelset_implicit_surface;
    implicit_surface->Update_Box();VECTOR_3D<double> size=implicit_surface->box.Size();
    double cell_size=.1*min(size.x,size.y,size.z);
    int m=(int)ceil(size.x/cell_size),n=(int)ceil(.5*size.y/cell_size),mn=(int)ceil(size.z/cell_size);
    double contraction_fraction=.7; 
    meshing_grid.Initialize(m,n,mn,implicit_surface->box.Center().x - contraction_fraction*size.x/2,implicit_surface->box.Center().x + contraction_fraction*size.x/2,
                                implicit_surface->box.Center().y,                                implicit_surface->box.Center().y + contraction_fraction*size.y/2,
                                implicit_surface->box.Center().z - contraction_fraction*size.z/2,implicit_surface->box.Center().z + contraction_fraction*size.z/2);
    embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Delete_Velocity_And_Acceleration();// in case they're accidently stored in the .etv file
    embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Delete_Mass(); // in case they're accidently stored in the .etv file
    embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Store_Position();
    embedded_tetrahedralized_volume.tetrahedralized_volume.Initialize_Cube_Mesh_And_Particles(meshing_grid);    
    embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Update_Position_And_Velocity();
    embedded_tetrahedralized_volume.tetrahedralized_volume.particles.Store_Mass();
    embedded_tetrahedralized_volume.Calculate_Phi_On_Tetrahedron_Vertices(*implicit_surface);    
    //Remove_Tetrahedra_Outside_Levelset(embedded_tetrahedralized_volume,true); //this must be done after phi on vertices, but before calculating triangulated surface
    //embedded_tetrahedralized_volume.Calculate_Phi_On_Tetrahedron_Vertices(*implicit_surface);    
    Fracture_Along_Level_Set(embedded_tetrahedralized_volume);    
    embedded_tetrahedralized_volume.Calculate_Triangulated_Surface();
    Write_Tetrahedralized_Volume(embedded_tetrahedralized_volume,"initial_mesh.tet");
    Write_Triangulated_Surface(embedded_tetrahedralized_volume,"initial_tri_surf.tri");
    Write_Embedded_Tetrahedralized_Volume(embedded_tetrahedralized_volume,"initial_tet_vol.etv");


    //shift above floor && scale
    int i;
    for(i=1;i<=embedded_tetrahedralized_volume.tetrahedralized_volume.particles.array_collection->Size();i++) {embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i).y+=3;embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i)*=.7;}
    embedded_tetrahedralized_volume.tetrahedralized_volume.Update_Bounding_Box();
    
    if(restart_step_number == 0){
        VECTOR_3D<double> center(embedded_tetrahedralized_volume.tetrahedralized_volume.bounding_box->Center());double bottom=embedded_tetrahedralized_volume.tetrahedralized_volume.bounding_box->ymin;
        for(i=1;i<=embedded_tetrahedralized_volume.tetrahedralized_volume.particles.array_size;i++){
            embedded_tetrahedralized_volume.tetrahedralized_volume.particles.V(i)=initial_velocity+VECTOR_3D<double>::Cross_Product(initial_angular_velocity,embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i)-center);
            embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i)=center+initial_orientation.Rotate(embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i)-center);
            embedded_tetrahedralized_volume.tetrahedralized_volume.particles.X(i).y+=initial_height-bottom;}}
    Get_Mass_Of_Each_Node(embedded_tetrahedralized_volume.tetrahedralized_volume);
    if(verbose) {std::cout << "embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m=" << embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m << std::endl;}

}
//#####################################################################
// Function Initialize_Attachment_Forces_And_Velocities
//#####################################################################
virtual void Initialize_Attachment_Forces_And_Velocities(COLLISION_FORCES_AND_VELOCITIES& collision_forces_and_velocities)
{
    int constraint_type=3;
    int m=meshing_grid.m,n=meshing_grid.n,mn=meshing_grid.mn;
       if(constraint_type==3){
        for(int i=1;i<=meshing_grid.m;i++)
            for(int k=1;k<=meshing_grid.mn;k++){
                collision_forces_and_velocities.enslaved_nodes.Append((k-1)*m*n + m*(n -1) + i);
            }
    }
    else if(constraint_type==2){
        //for(int i=1;i<=m_input;i++)
        //    attachment_constraints.enslaved_nodes.Append(i);
    }
    else if(constraint_type==1){
        //attachment_constraints.enslaved_nodes.Append(1);
    }
    else return;
}
//#####################################################################
// Function Fracture_Along_Level_Set
//#####################################################################
virtual void Fracture_Along_Level_Set(EMBEDDED_TETRAHEDRALIZED_VOLUME& embedded_tetrahedralized_volume)
{
    PARTICLE_3D<double>& particles=embedded_tetrahedralized_volume.tetrahedralized_volume.particles;
    TETRAHEDRON_MESH& mesh=embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh;
    ARRAY<int>particle_replicated(1,particles.array_collection->Size());
    ARRAY<int>tetrahedron_label(1,mesh.tetrahedrons.m);
    if(particle_replicated.m != particles.array_collection->Size()) particle_replicated.Resize(1,particles.array_collection->Size());
    if(tetrahedron_label.m != mesh.tetrahedrons.m) {std::cout << "error with tetrahedron_label passed" << std::endl;return;}
    int incident_tetrahedrons_defined=(int)mesh.incident_tetrahedrons;if(!incident_tetrahedrons_defined) mesh.Initialize_Incident_Tetrahedrons();
    int number_on_positive_side=0,number_on_negative_side=0,number_intersecting_boundary=0;
    int t,p;
    for(t=1;t<=mesh.tetrahedrons.m;t++){
        int i,j,k,l;mesh.tetrahedrons.Get(t,i,j,k,l);
        double phi1=embedded_tetrahedralized_volume.phi(i),phi2=embedded_tetrahedralized_volume.phi(j),
               phi3=embedded_tetrahedralized_volume.phi(k),phi4=embedded_tetrahedralized_volume.phi(l);
        int positive_count=(phi1>0)+(phi2>0)+(phi3>0)+(phi4>0);
        if(positive_count == 4) tetrahedron_label(t)=1;
        else if(positive_count == 0) tetrahedron_label(t)=-1;
        else tetrahedron_label(t)=0;
    }
    
    
    for(t=1;t<=mesh.tetrahedrons.m;t++){
        if(tetrahedron_label(t) == 1) number_on_positive_side++;
        else if(tetrahedron_label(t) == -1) number_on_negative_side++;
        else if(tetrahedron_label(t) ==  0) number_intersecting_boundary++;
        else{std::cout << "error with tetrahedron_label passed" << std::endl;return;}
    }
    int number_of_new_particles=0;
    for(t=1;t<=mesh.tetrahedrons.m;t++){
        if(tetrahedron_label(t) == 0){
            for(int a=1;a<=4;a++){
                int node=mesh.tetrahedrons(a,t);
                if(!particle_replicated(node)){
                    int new_index=particles.array_collection->Add_Element();assert(new_index == particles.array_collection->Size());particle_replicated(node)=new_index;number_of_new_particles++;
                    particles.X(new_index)=particles.X(node);particles.V(new_index)=particles.V(node);// still need to set mass
                }
            }
        }
    }
    
    for(p=1;p<=particle_replicated.m;p++){
        if(particle_replicated(p)){
            for(t=1;t<=(*mesh.incident_tetrahedrons)(p).m;t++){
                int tetrahedron=(*mesh.incident_tetrahedrons)(p)(t);
                if(tetrahedron_label(tetrahedron) == -1){
                    int i,j,k,l;mesh.tetrahedrons.Get(tetrahedron,i,j,k,l);
                    if(p == i) mesh.tetrahedrons(1,tetrahedron)=particle_replicated(p);
                    else if(p == j) mesh.tetrahedrons(2,tetrahedron)=particle_replicated(p);
                    else if(p == k) mesh.tetrahedrons(3,tetrahedron)=particle_replicated(p);
                    else if(p == l) mesh.tetrahedrons(4,tetrahedron)=particle_replicated(p);
                    else assert(false);
                }
            }
        }
    }

    for(t=1;t<=tetrahedron_label.m;t++){
        if(tetrahedron_label(t) == 0){
            int i,j,k,l;mesh.tetrahedrons.Get(t,i,j,k,l);
            mesh.tetrahedrons.Append(particle_replicated(i),particle_replicated(j),particle_replicated(k),particle_replicated(l));
        }
    }

    //repair phi
    embedded_tetrahedralized_volume.phi.Resize(particles.array_collection->Size());
    for(p=1;p<=particle_replicated.m;p++){
        if(particle_replicated(p)) embedded_tetrahedralized_volume.phi(particle_replicated(p))=embedded_tetrahedralized_volume.phi(p);
    }

    mesh.number_nodes=particles.array_collection->Size();
    if(!incident_tetrahedrons_defined){delete mesh.incident_tetrahedrons;mesh.incident_tetrahedrons=0;}
}
//#####################################################################
};
#endif

//#####################################################################
// Copyright 2004, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARMADILLO_DROP_HEAD_EXAMPLE
//##################################################################### 
#ifndef __ARMADILLO_DROP_HEAD_EXAMPLE__
#define __ARMADILLO_DROP_HEAD_EXAMPLE__

#include "../CLOTH_EXAMPLE.h"
namespace PhysBAM{

template<class T>
class ARMADILLO_DROP_HEAD_EXAMPLE:public CLOTH_EXAMPLE<T>
{
public:
    T height;
    QUATERNION<T> initial_orientation;

    ARMADILLO_DROP_HEAD_EXAMPLE()
        :height(0),initial_orientation((T)3.14,VECTOR_3D<T>(0,0,1))
    {
        use_diagonalized_fvm=false;
        use_masses_and_springs=true;use_altitude_springs=true;
        //edge_spring_restlength_enlargement_fraction=(T).01;
        //altitude_spring_restlength_enlargement_fraction=(T).01;
        use_bending_elements=true;use_bending_springs=false;
        final_time=7;frame_rate=30;
        restart_step_number=0;
        cfl_number=(T)200;
        max_bend_strain_per_time_step=(T).2;
        cg_tolerance=(T)1e-2;cg_iterations=2000;
        youngs_modulus=1e6;poissons_ratio=(T).3;Rayleigh_coefficient=(T).02;
        edge_spring_stiffness=(T)10;edge_spring_overdamping_fraction=4;
        altitude_spring_stiffness=(T)10;altitude_spring_overdamping_fraction=4;
        bending_stiffness=(T)1;bending_damping=(T).02;bending_cutoff_fraction_of_minimum_area=0/*5*/;bending_cutoff_fraction_of_triangles=(T).01;
        bending_spring_stiffness=(T)10/*800*/;bending_spring_overdamping_fraction=4;
        //bending_spring_restlength_enlargement_fraction=0;     
        output_directory="Armadillo_Drop_Head/output";
        check_initial_mesh_for_self_intersection=true;
        perform_self_collision=false;
        allow_intersections=true;allow_intersections_tolerance=(T)1e-8;
        print_residuals=false;
        verbose_dt=true;
    
        gravity=(T)9.8;
    }

    ~ARMADILLO_DROP_HEAD_EXAMPLE()
    {}

//#####################################################################
// Function Initialize_Triangulated_Surface
//#####################################################################
virtual void Get_Initial_Data(TRIANGULATED_SURFACE<T>& triangulated_surface)
{
    std::fstream input;char filename[256];
    TRIANGLE_MESH& triangle_mesh=triangulated_surface.triangle_mesh;PARTICLES<T,VECTOR_3D<T> >& particles=triangulated_surface.particles;        
    sprintf(filename,"../../Public_Data/Triangulated_Surfaces/armadillo_40K.tri");
    input.open(filename,std::ios::in|std::ios::binary);assert(input.is_open());
    triangulated_surface.Read_Float(input);input.close();    
    triangulated_surface.particles.Delete_Velocity_And_Acceleration();triangulated_surface.particles.Delete_Mass(); // in case they're accidently stored in the .tri file
    particles.Store_Position_And_Velocity();particles.Store_Mass();particles.Update_Position_And_Velocity();
    triangulated_surface.Update_Bounding_Box();
    VECTOR_3D<T> center(triangulated_surface.bounding_box->Center());
    //T mass_node=sqr(side_length);ARRAY<T,VECTOR<int,1> >::copy(mass_node,particles.mass);
    T scale=(T).01;
    T velocity=0;//sqrt(2*9.8*height);
    T miny=1e10;
    for(int p=1;p<=particles.array_collection->Size();p++) if(particles.X(p).y<miny) miny=particles.X(p).y;
    for(int p=1;p<=particles.array_collection->Size();p++){
        particles.X(p)=center+initial_orientation.Rotate(particles.X(p)-center);
        particles.X(p)=scale*(particles.X(p)-VECTOR_3D<T>(0,miny,0));
        particles.V(p)=VECTOR_3D<T>(0,-velocity,0);}
    triangulated_surface.Set_Density(8);
    triangulated_surface.Set_Mass_Of_Particles(false);

    LOG::SCOPE scope("mesh statistics");
    triangulated_surface.Print_Statistics(LOG::cout);
}
//#####################################################################
// Function Initialize_Collision_Bodies
//#####################################################################
void Initialize_Collision_Bodies()
{
    int index=solids_parameters.rigid_body_parameters.list.Add_Rigid_Body("../../Public_Data/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.collision_body_list.Append_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2003, Neil Molino, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ONE_TETRAHEDRON_EXAMPLE
//#####################################################################
#ifndef __ONE_TETRAHEDRON_EXAMPLE__
#define __ONE_TETRAHEDRON_EXAMPLE__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include "../TET_SIM_EMBEDDED_EXAMPLE.h"
using namespace PhysBAM;

template <class T>
class ONE_TETRAHEDRON_EXAMPLE:public TET_SIM_EMBEDDED_EXAMPLE<T>
{
public:
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    T xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input; 
    
    ONE_TETRAHEDRON_EXAMPLE()
                                :initial_height((T)1.2),initial_orientation(0,VECTOR_3D<T>(1,0,0)),
                                initial_velocity(0,(T)-sqrt((T)80),0),initial_angular_velocity(10,0,0),     
                                xmin_input((T)-.5),xmax_input((T)-xmin_input),ymin_input((T)-.5),ymax_input((T)-ymin_input),zmin_input((T)-.5),zmax_input((T)-zmin_input)
                                
    {
        //restart_step_number=8;
        final_time=5;                                         
        frame_rate=50;
        cfl_number=.5;
        cg_iterations=250;
        youngs_modulus=(T)1.9e6;poissons_ratio=(T).4;Rayleigh_coefficient=(T).01;        
        strcpy(output_directory,"One_Tet/output");
        gravity=0;
        
        // Fracture Parameters
        perform_fracture=false;
        perform_self_collision=false;
        push_out=false;

        
        T scale =1;
        fracture_threshold_1=(T)1e4*scale;fracture_threshold_2=(T)5e4*scale;fracture_threshold_3=(T)8e4*scale;
        compressive_threshold_1=(T)-100*1e4*scale;compressive_threshold_2=(T)-5e4*scale;compressive_threshold_3=(T)-8e4*scale;
        interpolation_fraction_threshhold=(T).0000002;
        max_number_of_cuts=3;
        number_of_fracture_initiation_points=1000;
    }

    virtual ~ONE_TETRAHEDRON_EXAMPLE()
    {}

//#####################################################################
// Function Initialize_Tetrahedralized_Volume
//#####################################################################
void Initialize_Tetrahedralized_Volume(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume)
{

    GRID<TV> mattress_grid(2,2,2,xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input);
                                
    tetrahedralized_volume.particles.Update_Position();
    tetrahedralized_volume.Initialize_Cube_Mesh_And_Particles(mattress_grid);
    ARRAY<int>deletion_list(1,0);deletion_list.Append(1);deletion_list.Append(2);
    deletion_list.Append(3);deletion_list.Append(4);
    tetrahedralized_volume.tetrahedron_mesh.Delete_Tetrahedrons(deletion_list);
    tetrahedralized_volume.Discard_Valence_Zero_Particles_And_Renumber();

    tetrahedralized_volume.particles.Delete_Velocity_And_Acceleration();
    tetrahedralized_volume.particles.Delete_Mass(); // in case they're accidently stored in the .tet file
    tetrahedralized_volume.particles.Update_Position_And_Velocity();
    tetrahedralized_volume.particles.Store_Mass();

    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=1;i<=tetrahedralized_volume.particles.array_size;i++){
        tetrahedralized_volume.particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i)=center+initial_orientation.Rotate(tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i).y+=initial_height-bottom;}
    std::cout << "total vertices = " << tetrahedralized_volume.particles.array_collection->Size() << std::endl;
    std::cout << "total tets = " << tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m << std::endl;
    tetrahedralized_volume.Set_Density(1000);
    tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);
}
//#####################################################################
// Function Initialize_Reference_Fracture_Tetrahedralized_Volume
//#####################################################################
virtual void Initialize_Reference_Fracture_Tetrahedralized_Volume(const EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& etv,
                                                                  FRACTURE_TETRAHEDRALIZED_VOLUME<T>& rftv,ARRAY<int> &reference_orientation_index)
{
    TET_SIM_EMBEDDED_EXAMPLE<T>::Initialize_Reference_Fracture_Tetrahedralized_Volume(etv,rftv,reference_orientation_index);

    PRINT<T>::Print_Tetrahedron(etv.tetrahedralized_volume.tetrahedron_mesh,1);
    PRINT<T>::Print_Tetrahedron(rftv.embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh,1);
    std::cout << etv.tetrahedralized_volume.particles.array_collection->Size() << std::endl;
    
    reference_orientation_index(1)=1;
    
    int ri,rj,rk,rl;rftv.embedded_tetrahedralized_volume.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.Get(1,ri,rj,rk,rl);
    int ri_rj=rftv.embedded_tetrahedralized_volume.Add_Embedded_Particle_If_Not_Already_There(ri,rj,.5);
    int ri_rk=rftv.embedded_tetrahedralized_volume.Add_Embedded_Particle_If_Not_Already_There(ri,rk,.5);
    int ri_rl=rftv.embedded_tetrahedralized_volume.Add_Embedded_Particle_If_Not_Already_There(ri,rl,.5);
    int rj_rk=rftv.embedded_tetrahedralized_volume.Add_Embedded_Particle_If_Not_Already_There(rj,rk,.5);
    int rj_rl=rftv.embedded_tetrahedralized_volume.Add_Embedded_Particle_If_Not_Already_There(rj,rl,.5);
    int rk_rl=rftv.embedded_tetrahedralized_volume.Add_Embedded_Particle_If_Not_Already_There(rk,rl,.5);
    //rftv.embedded_tetrahedralized_volume.Add_Embedded_Triangle_If_Not_Already_There(ri_rj,ri_rk,ri_rl);
    //rftv.embedded_tetrahedralized_volume.Add_Embedded_Triangle_If_Not_Already_There(rj_rk,rj_rl,ri_rj);
    //rftv.embedded_tetrahedralized_volume.Add_Embedded_Triangle_If_Not_Already_There(rk_rl,ri_rk,rj_rk);
    //rftv.embedded_tetrahedralized_volume.Add_Embedded_Triangle_If_Not_Already_There(rk_rl,rj_rl,ri_rl);

    //first quad added
    rftv.embedded_tetrahedralized_volume.Add_Embedded_Triangle_If_Not_Already_There(ri_rk,rj_rk,rj_rl);
    rftv.embedded_tetrahedralized_volume.Add_Embedded_Triangle_If_Not_Already_There(rj_rl,ri_rl,ri_rk);

    //second quad added
    rftv.embedded_tetrahedralized_volume.Add_Embedded_Triangle_If_Not_Already_There(ri_rk,ri_rj,rj_rl);
    rftv.embedded_tetrahedralized_volume.Add_Embedded_Triangle_If_Not_Already_There(rj_rl,rk_rl,ri_rk);
}
//##################################################################
};
#endif

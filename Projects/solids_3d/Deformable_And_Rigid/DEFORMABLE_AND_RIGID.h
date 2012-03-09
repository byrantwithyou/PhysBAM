//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_AND_RIGID
//#####################################################################
#ifndef __DEFORMABLE_AND_RIGID__
#define __DEFORMABLE_AND_RIGID__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Constitutive_Models/DIAGONALIZED_NEO_HOOKEAN_3D.h>
#include <Forces_And_Torques/BODY_FORCES_3D.h>
namespace PhysBAM{

template<class T,class RW>
class DEFORMABLE_AND_RIGID:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::last_frame;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;

    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    std::string input_file;

    DEFORMABLE_AND_RIGID()
        :BASE(0,fluids_parameters.NONE),initial_height(1),initial_orientation(0,VECTOR_3D<T>(1,0,0)),initial_velocity(0,0,0),initial_angular_velocity(0,0,0)
    {   
        solids_parameters.collisions_repulsion_thickness=(T)1e-2;
        solids_parameters.collisions_repulsion_clamp_fraction=(T).9;
        solids_parameters.collision_repulsion_spring_multiplier=100;

        last_frame=10*24;
        //last_frame=17;//10*24;
        restart=false;restart_frame=16;   
        solids_parameters.cfl=(T)10;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        output_directory="Deformable_And_Rigid/output";
        input_file=data_directory+"/Tetrahedralized_Volumes/adaptive_torus_float.tet";
        solids_parameters.perform_self_collision=false;
        solids_parameters.collide_with_interior=false;

        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    }

    ~DEFORMABLE_AND_RIGID()
    {}

    void Print_Mesh_Information()
    {std::cout << "minimum volume = " << solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume->Minimum_Signed_Volume() << std::endl;}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Tetrahedralized_Volume();
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*solids_parameters.deformable_body_parameters.list(index).tetrahedralized_volume;
    TETRAHEDRON_MESH& tetrahedron_mesh=tetrahedralized_volume.tetrahedron_mesh;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles=tetrahedralized_volume.particles;

    FILE_UTILITIES::Read_From_File<RW>(input_file,tetrahedralized_volume);
    std::cout << "total vertices = " << particles.array_collection->Size() << std::endl;std::cout << "total tetrahedra = " << tetrahedron_mesh.tetrahedrons.m << std::endl;
    particles.Update_Velocity();particles.Store_Mass(); // in case they're not stored in the file!

    tetrahedralized_volume.Set_Density(1000);tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);
    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=0;i<particles.array_collection->Size();i++){
        particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,particles.X(i)-center);
        particles.X(i)=center+initial_orientation.Rotate(particles.X(i)-center);
        particles.X(i).y+=initial_height-bottom;}

    index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/sphere");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->coefficient_of_friction=(T).3;
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->position=VECTOR_3D<T>(-3,2,0);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->velocity=VECTOR_3D<T>(10,0,0);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,
        solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

    index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->coefficient_of_friction=(T).3;
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->is_static=true;

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume;

    solids_parameters.deformable_body_parameters.list(1).Add_Force(Create_Body_Forces<T>(tetrahedralized_volume));
    
    //solids_parameters.deformable_body_parameters.list(1).Add_Edge_Springs(tetrahedralized_volume.tetrahedron_mesh,3000,2);
    //solids_parameters.deformable_body_parameters.list(1).Add_Altitude_Springs(tetrahedralized_volume.tetrahedron_mesh,300);
    //solids_parameters.deformable_body_parameters.list(1).Add_Linear_Elasticity(tetrahedralized_volume,2e5,.45,.01);
    //solids_parameters.deformable_body_parameters.list(1).Add_Neo_Hookean_Elasticity(tetrahedralized_volume,(T)2e5,(T).45,(T).01);
    solids_parameters.deformable_body_parameters.list(1).Add_Force(Create_Diagonalized_Finite_Volume(tetrahedralized_volume,new DIAGONALIZED_NEO_HOOKEAN_3D<T>((T)2e5,(T).45,(T).01)));
    /*    
    solids_parameters.deformable_body_parameters.list(1).Add_Diagonalized_Linear_Finite_Volume(tetrahedralized_volume,(T)3e5,(T).3,(T).01);
    solids_parameters.deformable_body_parameters.list(1).Disable_Finite_Volume_Damping();
    solids_parameters.deformable_body_parameters.list(1).Add_Edge_Springs(tetrahedralized_volume.tetrahedron_mesh,3000,(T)1);
    solids_parameters.deformable_body_parameters.list(1).Add_Altitude_Springs(tetrahedralized_volume.tetrahedron_mesh,300,(T)1);
    solids_parameters.deformable_body_parameters.list(1).Disable_Spring_Elasticity();
    */
}
//#####################################################################
};
}
#endif

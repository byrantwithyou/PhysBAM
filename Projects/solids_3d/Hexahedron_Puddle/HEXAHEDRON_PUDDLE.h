//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HEXAHEDRON_PUDDLE
//#####################################################################
#ifndef __HEXAHEDRON_PUDDLE__
#define __HEXAHEDRON_PUDDLE__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Forces_And_Torques/DIAGONALIZED_FINITE_VOLUME_HEXAHEDRONS.h>
namespace PhysBAM{

template<class T,class RW>
class HEXAHEDRON_PUDDLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR_3D<T>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR_3D<T>,RW> BASE;
    using BASE::last_frame;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::data_directory;using BASE::fluids_parameters;using BASE::solids_parameters;

    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    T mu;

    HEXAHEDRON_PUDDLE()
        :BASE(fluids_parameters.NONE),initial_height(1),initial_orientation(0,VECTOR_3D<T>(1,0,0)),initial_velocity(0,0,0),initial_angular_velocity(.5,.25,0)
    {
        solids_parameters.collisions_repulsion_thickness=(T)1e-2;
        solids_parameters.collisions_repulsion_clamp_fraction=(T).9;
        solids_parameters.collision_repulsion_spring_multiplier=100;

        last_frame=200;
        restart=true;restart_frame=46;
        solids_parameters.cfl=(T)3;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        output_directory="Hexahedron_Puddle/output";
        solids_parameters.perform_self_collision=true;
        solids_parameters.collide_with_interior=true;
    }

    ~HEXAHEDRON_PUDDLE()
    {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Hexahedralized_Volume();
    HEXAHEDRALIZED_VOLUME<T>& hexahedralized_volume=*solids_parameters.deformable_body_parameters.list(index).hexahedralized_volume;
    HEXAHEDRON_MESH& hexahedron_mesh=hexahedralized_volume.hexahedron_mesh;
    PARTICLES<T,VECTOR_3D<T> >& particles=hexahedralized_volume.particles;

    int resolution=15;
    GRID<TV> grid(3*resolution+1,3*resolution+1,resolution+1,(T)-1,(T)1,(T)-1,(T)1,-(T)one_third,(T)one_third);
    hexahedralized_volume.Initialize_Cube_Mesh_And_Particles(grid);
    T radius=(T)(one_third-1e-5);
    for(int h=1;h<=hexahedron_mesh.hexahedrons.m;h++)for(int k=1;k<=8;k++){
        VECTOR_3D<T>& X=particles.X(hexahedron_mesh.hexahedrons(k,h));
        if(fabs(X.x)<radius && fabs(X.y)<radius)hexahedron_mesh.hexahedrons(k,h)=0;}
    hexahedron_mesh.Delete_Hexahedrons_With_Missing_Nodes();
    hexahedralized_volume.Discard_Valence_Zero_Particles_And_Renumber();
    
    hexahedralized_volume.Set_Density(1000);hexahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);
    hexahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(hexahedralized_volume.bounding_box->Center());T bottom=hexahedralized_volume.bounding_box->ymin;
    for(int i=1;i<=particles.array_collection->Size();i++){
        particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,particles.X(i)-center);
        particles.X(i)=center+initial_orientation.Rotate(particles.X(i)-center);
        particles.X(i).y+=initial_height-bottom;}

    index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->is_static=true;

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    Get_Initial_Data();
    DEFORMABLE_OBJECT_3D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(1);
    HEXAHEDRALIZED_VOLUME<T>& hexahedralized_volume=*deformable_object.hexahedralized_volume;
    deformable_object.Add_Body_Forces(hexahedralized_volume);
    deformable_object.Add_Diagonalized_Linear_Finite_Volume(hexahedralized_volume,(T)1e5,(T).45,(T).01);
    mu=deformable_object.diagonalized_finite_volume_hexahedrons(1)->constitutive_model.mu;
}
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
void Update_Time_Varying_Material_Properties(const T time)
{
    T critical=(T)1.6;
    DEFORMABLE_OBJECT_3D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(1);
    DIAGONALIZED_CONSTITUTIVE_MODEL_3D<T>& constitutive_model=deformable_object.diagonalized_finite_volume_hexahedrons(1)->constitutive_model;
    if(time<critical) constitutive_model.mu=0;
    else constitutive_model.mu=mu;
}
//#####################################################################
};
}
#endif

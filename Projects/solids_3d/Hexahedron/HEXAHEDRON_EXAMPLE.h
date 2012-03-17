//#####################################################################
// Copyright 2005, Geoffrey Irving, Tamar Shinar, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HEXAHEDRON_EXAMPLE
//#####################################################################
#ifndef __HEXAHEDRON_EXAMPLE__
#define __HEXAHEDRON_EXAMPLE__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Geometry/Topology/HEXAHEDRON_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Constitutive_Models/DIAGONALIZED_NEO_HOOKEAN_3D.h>
#include <Forces_And_Torques/BODY_FORCES_3D.h>
#include <Forces_And_Torques/DIAGONALIZED_FINITE_VOLUME_3D.h>
namespace PhysBAM{

template<class T,class RW>
class HEXAHEDRON_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::last_frame;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::data_directory;using BASE::fluids_parameters;using BASE::solids_parameters;

    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    std::string input_file;
    VECTOR_3D<T> gravity;
    T max_strain;
    int m,n,mn;
    T xmax,ymax,zmax;
    bool compare_with_tets,read_geometry_from_file;
    T youngs_modulus,poisson_ratio,rayleigh_coeffecient;

    HEXAHEDRON_EXAMPLE()
        :BASE(0,fluids_parameters.NONE),initial_height(.2),initial_orientation(0,VECTOR_3D<T>(1,0,0)),initial_velocity(0,0,0),initial_angular_velocity(0,0,0),max_strain((T).1),
        m(30),n(30),mn(30),xmax((T).1),ymax((T).1),zmax((T).1),compare_with_tets(false),read_geometry_from_file(true)
    {   
        solids_parameters.collisions_repulsion_thickness=(T)1e-2;
        solids_parameters.collisions_repulsion_clamp_fraction=(T).9;
        solids_parameters.collision_repulsion_spring_multiplier=100;
        solids_parameters.use_constant_mass=true;

        last_frame=100;
        restart=false;restart_frame=99;
        solids_parameters.cfl=(T)3;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        output_directory="Hexahedron/output";
        solids_parameters.perform_self_collision=false;
        solids_parameters.collide_with_interior=true;
        youngs_modulus=(T)2e3;poisson_ratio=(T).45;rayleigh_coeffecient=(T).003;
    }

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Hexahedralized_Volume();
    HEXAHEDRALIZED_VOLUME<T>& hexahedralized_volume=*solids_parameters.deformable_body_parameters.list(index).hexahedralized_volume;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& hex_particles=hexahedralized_volume.particles;HEXAHEDRON_MESH& hexahedron_mesh=hexahedralized_volume.hexahedron_mesh;

    if(read_geometry_from_file){
        FILE_UTILITIES::Read_From_File<T>(data_directory+"/Delp_Muscles/idealized_biceps.hex.gz",hexahedralized_volume);
        hex_particles.Store_Velocity(true);hex_particles.Store_Mass(false);hex_particles.Store_Mass(true);}
    else{
        GRID<TV> grid(n,m,mn,0,xmax,0,ymax,0,zmax);
        hexahedralized_volume.Initialize_Cube_Mesh_And_Particles(grid);}

    std::cout << "total vertices = " << hex_particles.Size() << std::endl;std::cout << "total hexahedra = " << hexahedron_mesh.hexahedrons.m << std::endl;
    hexahedralized_volume.Set_Density(1000);hexahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);
    hexahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(hexahedralized_volume.bounding_box->Center());T bottom=hexahedralized_volume.bounding_box->ymin;
    for(int i=0;i<hex_particles.Size();i++){
        hex_particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,hex_particles.X(i)-center);
        hex_particles.X(i)=center+initial_orientation.Rotate(hex_particles.X(i)-center);
        hex_particles.X(i).y+=initial_height-bottom;}

    if(compare_with_tets){
        index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Tetrahedralized_Volume();
        TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*solids_parameters.deformable_body_parameters.list(index).tetrahedralized_volume;
        DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& tet_particles=tetrahedralized_volume.particles;
        tet_particles.Initialize_Particles(hex_particles);HEXAHEDRALIZED_VOLUME<T> hex_copy(hexahedralized_volume.hexahedron_mesh,tet_particles);
        hex_copy.Initialize_Tetrahedralized_Volume();
        tetrahedralized_volume.tetrahedron_mesh.Initialize_Tetrahedron_Mesh(hex_copy.tetrahedralized_volume->tetrahedron_mesh.tetrahedrons);
        tetrahedralized_volume.Set_Density(1000);tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);}

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
    solid_body_collection.Add_Force(Create_Body_Forces<T>(hexahedralized_volume));
    solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(hexahedralized_volume,new DIAGONALIZED_NEO_HOOKEAN_3D<T>(youngs_modulus,poisson_ratio,rayleigh_coeffecient)));
    if(compare_with_tets){
        DEFORMABLE_OBJECT_3D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(2);
        TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*deformable_object.tetrahedralized_volume;
        solid_body_collection.Add_Force(Create_Body_Forces<T>(tetrahedralized_volume));
        solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(tetrahedralized_volume,new DIAGONALIZED_NEO_HOOKEAN_3D<T>(youngs_modulus,poisson_ratio,rayleigh_coeffecient)));}
}
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TORUS_EXAMPLE
//#####################################################################
#ifndef __TORUS_EXAMPLE__
#define __TORUS_EXAMPLE__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TETRAHEDRON_COLLISION_BODY.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Constitutive_Models/DIAGONALIZED_NEO_HOOKEAN_3D.h>
#include <Deformable_Objects/DEFORMABLE_OBJECT_ASYNCHRONOUS_3D.h>
#include <Forces_And_Torques/BODY_FORCES_3D.h>
#include <Geometry/LEVELSET_IMPLICIT_SURFACE.h>
namespace PhysBAM{

template<class T,class RW>
class TORUS_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::last_frame;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::data_directory;using BASE::fluids_parameters;
    using BASE::solids_parameters;

    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    std::string input_file;
    VECTOR_3D<T> gravity;
    T max_strain;
    int force,tweak;

    TETRAHEDRON_COLLISION_BODY<T>* torus_collision_body;
    TRIANGULATED_SURFACE<T>* undeformed_triangulated_surface;
    LEVELSET_IMPLICIT_SURFACE<T>* undeformed_levelset;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> > undeformed_positions;
    bool use_tetrahedron_collisions;
    T bounding_box_percentage_for_levelset_grid_size;

    TORUS_EXAMPLE()
        :BASE(0,fluids_parameters.NONE),initial_height(1),initial_orientation(0,VECTOR_3D<T>(1,0,0)),initial_velocity(0,0,0),initial_angular_velocity(0,0,0),max_strain(.1),
        torus_collision_body(0),undeformed_triangulated_surface(0),undeformed_levelset(0),use_tetrahedron_collisions(false),bounding_box_percentage_for_levelset_grid_size((T).01)
    {   
        solids_parameters.collisions_repulsion_thickness=(T)1e-2;
        solids_parameters.collisions_repulsion_clamp_fraction=(T).9;
        solids_parameters.collision_repulsion_spring_multiplier=100;

        force=1;
        tweak=1;

        last_frame=100;
        restart=false;restart_frame=99;
        solids_parameters.cfl=(T)10;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        input_file=data_directory+"/Tetrahedralized_Volumes/adaptive_torus_float.tet";
        output_directory=STRING_UTILITIES::string_sprintf("Torus/output%d_%d",force,tweak);
        solids_parameters.perform_self_collision=true;
        solids_parameters.collide_with_interior=false;

        gravity=solids_parameters.gravity*solids_parameters.gravity_direction;
        use_tetrahedron_collisions=true;
    }

    ~TORUS_EXAMPLE()
    {delete undeformed_triangulated_surface;delete torus_collision_body;delete undeformed_levelset;}

    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE
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

    std::istream* input=FILE_UTILITIES::Safe_Open_Input(input_file);tetrahedralized_volume.template Read<RW>(*input);delete input;
    std::cout << "total vertices = " << particles.Size() << std::endl;std::cout << "total tetrahedra = " << tetrahedron_mesh.tetrahedrons.m << std::endl;
    particles.Update_Velocity();particles.Store_Mass(); // in case they're not stored in the file!

    tetrahedralized_volume.Set_Density(1000);tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);
    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=0;i<particles.Size();i++){
        particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,particles.X(i)-center);
        particles.X(i)=center+initial_orientation.Rotate(particles.X(i)-center);
        particles.X(i).y+=initial_height-bottom;}

    index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->is_static=true;

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    if(use_tetrahedron_collisions){
        solids_parameters.perform_self_collision=false;
        undeformed_positions.Add_Elements(particles.Size());
        ARRAY<VECTOR_3D<T> >::copy(particles.X.array,undeformed_positions.X.array);
        this->Initialize_Tetrahedron_Collisions();}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    Get_Initial_Data();
    DEFORMABLE_OBJECT_3D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(1);
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*deformable_object.tetrahedralized_volume;
    solid_body_collection.Add_Force(Create_Body_Forces<T>(tetrahedralized_volume));
    
    //deformable_object.Add_Edge_Springs(tetrahedralized_volume.tetrahedron_mesh,3000,2);
    //deformable_object.Add_Altitude_Springs(tetrahedralized_volume.tetrahedron_mesh,300);
    //deformable_object.Add_Linear_Elasticity(tetrahedralized_volume,2e5,.45,.01);
    //deformable_object.Add_Neo_Hookean_Elasticity(tetrahedralized_volume,(T)2e5,(T).45,(T).01);
    //deformable_object.Add_Diagonalized_Neo_Hookean_Elasticity(tetrahedralized_volume,(T)2e5,(T).45,(T).01);

    if(force==1){ // diagonalized finite_volume
        if(tweak==1){solids_parameters.cfl=10;max_strain=(T).1;} // 176.42
        solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(tetrahedralized_volume,new DIAGONALIZED_NEO_HOOKEAN_3D<T>((T)2e5,(T).45,(T).01,(T).25),true,max_strain));}
    if(force==2){ // semi-implicit diagonalized finite volume
        solids_parameters.semi_implicit=true;
        if(tweak==1){solids_parameters.cfl=10;max_strain=(T).1;} // 99.70
        if(tweak==2){solids_parameters.cfl=20;max_strain=(T).2;} // 83.77
        if(tweak==3){solids_parameters.cfl=5;max_strain=(T).05;} // 143.99
        solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(tetrahedralized_volume,new DIAGONALIZED_NEO_HOOKEAN_3D<T>((T)2e5,(T).45,(T).01,(T).25),true,max_strain));}
    if(force==3){ // asynchronous semi-implicit diagonalized finite volume
        solids_parameters.asynchronous=true;
        if(tweak==1){solids_parameters.cfl=1;max_strain=(T).05;} //
        solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(tetrahedralized_volume,new DIAGONALIZED_NEO_HOOKEAN_3D<T>((T)2e5,(T).45,(T).01,(T).25),true,max_strain));
        Setup_Asynchronous_Time_Stepping(deformable_object);
        deformable_object.asynchronous->status_report_average_steps=50;}
}
//#####################################################################
// Function Add_External_Impulses
//#####################################################################
void Add_External_Impulses(ARRAY<VECTOR_3D<T> >& V,const T time,const T dt) PHYSBAM_OVERRIDE
{
    V+=dt*gravity; // since semi-implicit doesn't notice body forces yet 
}
//#####################################################################
// Function Add_External_Impulse
//#####################################################################
void Add_External_Impulse(ARRAY<VECTOR_3D<T> >& V,const int node,const T time,const T dt) PHYSBAM_OVERRIDE
{
    V(node)+=dt*gravity; // since semi-implicit doesn't notice body forces yet 
}
//#####################################################################
// Function Initialize_Tetrahedron_Collisions
//#####################################################################
void Initialize_Tetrahedron_Collisions()
{
    TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume;
    tetrahedralized_volume->Initialize_Triangulated_Surface();
    torus_collision_body=new TETRAHEDRON_COLLISION_BODY<T>(*tetrahedralized_volume);
    undeformed_triangulated_surface=new TRIANGULATED_SURFACE<T>(tetrahedralized_volume->triangulated_surface->triangle_mesh,undeformed_positions);
    undeformed_triangulated_surface->Update_Triangle_List();undeformed_triangulated_surface->Initialize_Triangle_Hierarchy();

    // undeformed levelset
    std::string levelset_file="Torus/undeformed_levelset.phi";
    if(FILE_UTILITIES::File_Exists(levelset_file)) FILE_UTILITIES::Create_From_File<RW>(levelset_file,undeformed_levelset);
    else{
        undeformed_triangulated_surface->Update_Bounding_Box();BOX_3D<T>& box=*undeformed_triangulated_surface->bounding_box;
        undeformed_levelset=LEVELSET_IMPLICIT_SURFACE<T>::Create();
        GRID<TV>& grid=undeformed_levelset->levelset.grid;ARRAY<T,VECTOR<int,3> >& phi=undeformed_levelset->levelset.phi;
        grid=GRID_3D<T>::Create_Grid_Given_Cell_Size(box,bounding_box_percentage_for_levelset_grid_size*box.Edge_Lengths().Max(),false,5);
        phi.Resize(grid);
        LEVELSET_MAKER_UNIFORM<T> levelset_maker;
        levelset_maker.Verbose_Mode();
        levelset_maker.Set_Surface_Padding_For_Flood_Fill((T)1e-3);
        levelset_maker.Use_Fast_Marching_Method(true,0);
        levelset_maker.Compute_Level_Set(*undeformed_triangulated_surface,grid,phi);
        FILE_UTILITIES::Write_To_File<RW>(levelset_file,*undeformed_levelset);} 
    undeformed_levelset->Update_Box();
    torus_collision_body->Set_Implicit_Geometry(undeformed_levelset);
    torus_collision_body->Set_Undeformed_Triangulated_Surface(undeformed_triangulated_surface);
    solids_parameters.collision_body_list.Add_Body(torus_collision_body);
}
//#####################################################################
};
}
#endif

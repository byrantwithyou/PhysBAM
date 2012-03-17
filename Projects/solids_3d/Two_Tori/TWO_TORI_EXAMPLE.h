//#####################################################################
// Copyright 2004-2006, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TWO_TORI_EXAMPLE
//#####################################################################
#ifndef __TWO_TORI_EXAMPLE__
#define __TWO_TORI_EXAMPLE__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/TRIANGULATED_SURFACE_SIGNED_DISTANCE_UNIFORM.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TETRAHEDRON_COLLISION_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Forces_And_Torques/BODY_FORCES_3D.h>
namespace PhysBAM{

template<class T,class RW>
class TWO_TORI_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    int number_of_objects;
    T initial_height1,initial_height2;
    QUATERNION<T> initial_orientation1,initial_orientation2;
    VECTOR_3D<T> initial_velocity1,initial_angular_velocity1,initial_velocity2,initial_angular_velocity2;
    ARRAY<TETRAHEDRALIZED_VOLUME<T>*> ring;
    ARRAY<TETRAHEDRON_COLLISION_BODY<T>*> ring_collision_body;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> > undeformed_positions;
    LEVELSET_IMPLICIT_SURFACE<T>* undeformed_levelset;
    TRIANGULATED_SURFACE<T>* undeformed_triangulated_surface;
    ARRAY<int> intersecting_points;
    ARRAY<ARRAY<int> > tetrahedrons_near_points;
    ARRAY<ARRAY<int> > triangles_near_points;
    ARRAY<ARRAY<VECTOR_3D<T> > > outward_direction;
    ARRAY<ARRAY<T> > distances_to_surface;
    bool use_tetrahedron_collisions,read_implicit_surfaces;
    T bounding_box_percentage_for_levelset_grid_size;
    std::string input_file;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::last_frame;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;
    using BASE::data_directory;using BASE::write_output_files;

    TWO_TORI_EXAMPLE()
        :BASE(0,fluids_parameters.NONE),number_of_objects(2),
         initial_height1(1.5),initial_height2(5),initial_orientation1(0/*T(.15*pi)*/,VECTOR_3D<T>(1,0,0)),initial_orientation2(0/*T(.35*pi)*/,VECTOR_3D<T>(0,1,0)),
         initial_velocity1(0,0,0),initial_angular_velocity1(0,0/*1*/,0),initial_velocity2(0,0,0),initial_angular_velocity2(0,0,0/*-2*/),
         use_tetrahedron_collisions(true),read_implicit_surfaces(true),bounding_box_percentage_for_levelset_grid_size((T).01)
    {   
        solids_parameters.collisions_repulsion_thickness=(T)1e-2;
        solids_parameters.collisions_repulsion_clamp_fraction=(T).9;
        solids_parameters.collision_repulsion_spring_multiplier=100;

        last_frame=10*24;
        restart=false;restart_frame=19;   
        solids_parameters.cfl=(T)10;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        output_directory="Two_Tori/output";
        input_file=data_directory+"/Tetrahedralized_Volumes/adaptive_torus_float.tet";
        //input_file=data_directory+"/Tetrahedralized_Volumes/sphere.tet";
        solids_parameters.perform_self_collision=true;
        if(use_tetrahedron_collisions)solids_parameters.perform_self_collision=false;
        solids_parameters.collide_with_interior=false;
        solids_parameters.synchronize_multiple_objects=true;
        write_output_files=true;
        number_of_objects=2;
    }

    ~TWO_TORI_EXAMPLE()
    {delete undeformed_triangulated_surface;ring_collision_body.Delete_Pointers_And_Clean_Memory();}
//#####################################################################
// Function Initialize_Tetrahedron_Collisions
//#####################################################################
void Initialize_Tetrahedron_Collisions()
{
    for(int object=0;object<number_of_objects;object++){
        ring_collision_body.Append(new TETRAHEDRON_COLLISION_BODY<T>(*ring(object)));
        ring(object)->Initialize_Triangulated_Surface();}
    undeformed_triangulated_surface=new TRIANGULATED_SURFACE<T>(ring(1)->triangulated_surface->triangle_mesh,undeformed_positions);
    undeformed_triangulated_surface->Update_Triangle_List();undeformed_triangulated_surface->Initialize_Triangle_Hierarchy();
    
    //implicit surfaces
    std::string input_file=output_directory+"/torus_implicit.phi";std::istream* input=FILE_UTILITIES::Safe_Open_Input(input_file,true,false);
    if(input){
        undeformed_levelset=LEVELSET_IMPLICIT_SURFACE<T>::Create();
        undeformed_levelset->template Read<RW>(*input);}
    else{
        ring(1)->triangulated_surface->Update_Bounding_Box();
        BOX_3D<T>& box=*ring(1)->triangulated_surface->bounding_box;
        T largest_dimension=max(box.xmax-box.xmin,box.ymax-box.ymin,box.zmax-box.zmin);
        T dx=bounding_box_percentage_for_levelset_grid_size*largest_dimension;T extra_padding=(T)10*dx;
        int m=int((box.xmax-box.xmin)/dx),n=int((box.ymax-box.ymin)/dx),mn=int((box.zmax-box.zmin)/dx);
        undeformed_levelset=LEVELSET_IMPLICIT_SURFACE<T>::Create();
        undeformed_levelset->levelset.phi.Resize(1,m,1,n,1,mn);
        undeformed_levelset->levelset.grid.Initialize(m,n,mn,box.xmin-extra_padding,box.xmax+extra_padding,
            box.ymin-extra_padding,box.ymax+extra_padding,box.zmin-extra_padding,box.zmax+extra_padding);
        SIGNED_DISTANCE::Calculate(*ring(1)->triangulated_surface,undeformed_levelset->levelset.grid,undeformed_levelset->levelset.phi,true);
        FILE_UTILITIES::Create_Directory(output_directory);FILE_UTILITIES::Write_To_File<RW>(output_directory+"/torus_implicit.phi",*undeformed_levelset);}
    undeformed_levelset->Update_Box();
    for(int object=0;object<number_of_objects;object++){
        ring_collision_body(object)->Set_Implicit_Geometry(undeformed_levelset);
        ring_collision_body(object)->Set_Undeformed_Triangulated_Surface(undeformed_triangulated_surface);
        solids_parameters.collision_body_list.Add_Body(ring_collision_body(object));}
}
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    T initial_height[2]={initial_height1,initial_height2};
    VECTOR_3D<T> initial_velocity[2]={initial_velocity1,initial_velocity2};
    VECTOR_3D<T> initial_angular_velocity[2]={initial_angular_velocity1,initial_angular_velocity2};
    QUATERNION<T> initial_orientation[2]={initial_orientation1,initial_orientation2};

    for(int object=0;object<number_of_objects;object++){
        int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Tetrahedralized_Volume();
        solids_parameters.deformable_body_parameters.list.deformable_objects(index)->id_number=index;
        TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*solids_parameters.deformable_body_parameters.list(index).tetrahedralized_volume;
        ring.Append(solids_parameters.deformable_body_parameters.list(index).tetrahedralized_volume);
        TETRAHEDRON_MESH& tetrahedron_mesh=tetrahedralized_volume.tetrahedron_mesh;
        DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles=tetrahedralized_volume.particles;
        
        std::istream* input=FILE_UTILITIES::Safe_Open_Input(input_file);tetrahedralized_volume.template Read<RW>(*input);delete input;
        std::cout << "total vertices = " << particles.Size() << std::endl;std::cout << "total tetrahedra = " << tetrahedron_mesh.tetrahedrons.m << std::endl;
        particles.Update_Velocity();particles.Store_Mass(); // in case they're not stored in the file!

        if(number_of_objects==1){
            int old_number_of_particles=particles.Size();int old_number_of_tets=tetrahedron_mesh.tetrahedrons.m;
            for(int p=0;p<old_number_of_particles;p++){int p2=particles.Add_Element();particles.X(p2)=particles.X(p);particles.X(p2).y+=initial_height2-initial_height1;}
            for(int t=0;t<old_number_of_tets;t++){
                int i,j,k,l;tetrahedron_mesh.tetrahedrons.Get(t,i,j,k,l);
                tetrahedron_mesh.tetrahedrons.Append(i+old_number_of_particles,j+old_number_of_particles,k+old_number_of_particles,l+old_number_of_particles);}
            tetrahedron_mesh.number_nodes=particles.Size();}

        tetrahedralized_volume.Set_Density(1000);tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);
        tetrahedralized_volume.Update_Bounding_Box();
        VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
        for(int i=0;i<particles.Size();i++){
            particles.V(i)=initial_velocity[object-1]+VECTOR_3D<T>::Cross_Product(initial_angular_velocity[object-1],particles.X(i)-center);
            particles.X(i)=center+initial_orientation[object-1].Rotate(particles.X(i)-center);
            particles.X(i).y+=initial_height[object-1]-bottom;
            if(object==1 && use_tetrahedron_collisions){undeformed_positions.Add_Element();undeformed_positions.X(i)=particles.X(i);}}}
        
    int index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->is_static=true;

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    if(use_tetrahedron_collisions) Initialize_Tetrahedron_Collisions();
}
//#####################################################################
// Function Create_Singleton_Geometry
//#####################################################################
void Create_Singleton_Geometry(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,const T edge_length=(T).25)
{
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles=tetrahedralized_volume.particles;TETRAHEDRON_MESH& tetrahedron_mesh=tetrahedralized_volume.tetrahedron_mesh;
    particles.Increase_Array_Size(4);
    particles.X(particles.Add_Element())=VECTOR_3D<T>(0,0,0);
    particles.X(particles.Add_Element())=edge_length*VECTOR_3D<T>(1,0,0);
    particles.X(particles.Add_Element())=edge_length*VECTOR_3D<T>(0,1,0);
    particles.X(particles.Add_Element())=edge_length*VECTOR_3D<T>(0,0,1);
    tetrahedron_mesh.Clean_Memory();tetrahedron_mesh.number_nodes=4;tetrahedron_mesh.tetrahedrons.Exact_Resize(4,1);
    tetrahedron_mesh.tetrahedrons(1,1)=1;tetrahedron_mesh.tetrahedrons(2,1)=2;tetrahedron_mesh.tetrahedrons(3,1)=3;tetrahedron_mesh.tetrahedrons(4,1)=4;
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    Get_Initial_Data();
    for(int object=0;object<number_of_objects;object++){
        if(use_tetrahedron_collisions)solids_parameters.deformable_body_parameters.list(object).collisions.collision_body_list_id=object+1;
        TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*solids_parameters.deformable_body_parameters.list(object).tetrahedralized_volume;
        solids_parameters.deformable_body_parameters.list(object).Add_Force(Create_Body_Forces<T>(tetrahedralized_volume));
        solids_parameters.deformable_body_parameters.list(object).Add_Force(Create_Edge_Springs<T>(tetrahedralized_volume,5000,(T)1));
        solids_parameters.deformable_body_parameters.list(object).Add_Force(Create_Altitude_Springs<T>(tetrahedralized_volume,500,(T)1));}

        //solids_parameters.deformable_body_parameters.list(1).Add_Linear_Elasticity(tetrahedralized_volume,2e5,.45,.01);
        //solids_parameters.deformable_body_parameters.list(1).Add_Neo_Hookean_Elasticity(tetrahedralized_volume,(T)2e5,(T).45,(T).01);
        //solids_parameters.deformable_body_parameters.list(1).Add_Diagonalized_Neo_Hookean_Elasticity(tetrahedralized_volume1,(T)2e6,(T).45,(T).01);
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

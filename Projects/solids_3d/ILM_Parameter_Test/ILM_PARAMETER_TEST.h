//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ILM_PARAMETER_TEST
//#####################################################################
#ifndef __ILM_PARAMETER_TEST__
#define __ILM_PARAMETER_TEST__

#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>
namespace PhysBAM{

template <class T,class RW>
class ILM_PARAMETER_TEST:public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart_frame;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::data_directory;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::solids_parameters;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::verbose_dt;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::initial_time;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::frame_rate;

    int example_number;
    T cloth_density;
    GRID<TV> cloth_grid;
    T enslaved_square_size;
    ARRAY<int> enslaved_nodes;
    T cloth_height;

    bool use_edge_springs;
    T edge_springs_stiffness;
    T edge_springs_overdamping_fraction;
    bool edge_springs_limit_time_step_by_strain_rate;
    T edge_springs_max_strain_per_time_step;
    bool edge_springs_use_rest_state_for_strain_rate;
    T edge_springs_restlength_enlargement_fraction;

    bool use_altitude_springs;
    T altitude_springs_stiffness;
    T altitude_springs_overdamping_fraction;
    bool altitude_springs_use_shortest_only;
    bool altitude_springs_use_compressed_by_threshold_only;
    T altitude_springs_fraction_compression;
    bool altitude_springs_limit_time_step_by_strain_rate;
    T altitude_springs_max_strain_per_time_step;
    bool altitude_springs_use_rest_state_for_strain_rate;
    T altitude_springs_restlength_enlargement_fraction;

    ILM_PARAMETER_TEST(const int example_number_input=0)
        :SOLIDS_FLUIDS_EXAMPLE_3D<RW>(fluids_parameters.NONE),example_number(example_number_input)
    {
        last_frame=500;
        restart=false;restart_frame=0;   
        output_directory="ILM_Parameter_Test/output";
        verbose_dt=true;

        // dynSolver - Dynamics - solver properties
        solids_parameters.implicit_solve_parameters.cg_iterations=200;
        solids_parameters.implicit_solve_parameters.cg_tolerance=0.001;
        solids_parameters.cfl=(T)0.9;

        // dynSolver - Dynamics - collisions
        solids_parameters.collision_tolerance=1e-6;
        solids_parameters.enforce_tangential_collision_velocity=false;

        // dynSolver - Dynamics - self collisions
        solids_parameters.perform_self_collision=true;
        solids_parameters.turn_off_all_collisions=true;
        solids_parameters.collisions_collision_thickness=1e-6;

        solids_parameters.turn_off_all_repulsions=false;
        solids_parameters.collisions_repulsion_thickness=0.1;
        solids_parameters.collisions_repulsion_clamp_fraction=1;
        solids_parameters.collision_repulsion_spring_multiplier=0.1;
//        solids_parameters.collision_repulsion_spring_multiplier=1e3;
        solids_parameters.spring_limiter_fraction=0.1;
        solids_parameters.collisions_disable_repulsions_based_on_proximity_factor=1.1;
        solids_parameters.self_collision_friction_coefficient=0;
        solids_parameters.allow_intersections=false;
        solids_parameters.allow_intersections_tolerance=1e-8;
        solids_parameters.collisions_nonrigid_collision_attempts=20;
        solids_parameters.min_collision_loops=1;
        solids_parameters.max_collision_loops=64;

        // cloth - Dynamics - material properties
        cloth_density=0.4;

        // cloth - Dynamics - body forces
        solids_parameters.gravity=32.1522;

        // cloth - Dynamics - edge springs
        use_edge_springs=true;
        edge_springs_stiffness=40;
        edge_springs_overdamping_fraction=1;
        edge_springs_limit_time_step_by_strain_rate=true;
        edge_springs_max_strain_per_time_step=.1;
        edge_springs_use_rest_state_for_strain_rate=true;
        edge_springs_restlength_enlargement_fraction=0;
      
        // cloth - Dynamics - altitude springs
        use_altitude_springs=false;
        altitude_springs_stiffness=4;
        altitude_springs_overdamping_fraction=1;
        altitude_springs_use_shortest_only=true;
        altitude_springs_use_compressed_by_threshold_only=true;
        altitude_springs_fraction_compression=0.1;
        altitude_springs_limit_time_step_by_strain_rate=true;
        altitude_springs_max_strain_per_time_step=0.1;
        altitude_springs_use_rest_state_for_strain_rate=true;
        altitude_springs_restlength_enlargement_fraction=0;

        // cloth - Dynamics - bending elements
        
        // cloth - Dynamics - bending springs

        // cloth - Dynamics - self collision

//        VECTOR_2D<int> cloth_resolution(5,5);
        VECTOR_2D<int> cloth_resolution(21,21);
        VECTOR_2D<T> cloth_size(1,1);
        enslaved_square_size=.401;
        cloth_height=0.2;

        cloth_grid=GRID_2D<T>(cloth_resolution.x,cloth_resolution.y,
                              -cloth_size.x/2,cloth_size.x/2,
                              -cloth_size.y/2,cloth_size.y/2);
    }

    ~ILM_PARAMETER_TEST()
    {}

int Node_Index(const int i,const int j)
{
    return i+cloth_grid.m*(j-1);
}

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    int id;

    id=Add_Deformable_Object(cloth_grid,MATRIX<T,4>::Rotation_Matrix_X_Axis(pi/2));
    DEFORMABLE_OBJECT_3D<T>* deformable_object=solids_parameters.deformable_body_parameters.list.deformable_objects(id);
    deformable_object->triangulated_surface->Set_Density(cloth_density);
    deformable_object->triangulated_surface->Set_Mass_Of_Particles();
    deformable_object->Add_Body_Forces(*deformable_object->triangulated_surface,solids_parameters.gravity,solids_parameters.gravity_direction);

    if(use_edge_springs)
        deformable_object->Add_Edge_Springs(deformable_object->triangulated_surface->triangle_mesh,
            edge_springs_stiffness,edge_springs_overdamping_fraction,edge_springs_limit_time_step_by_strain_rate,
            edge_springs_max_strain_per_time_step,edge_springs_use_rest_state_for_strain_rate,
            edge_springs_restlength_enlargement_fraction);

    if(use_altitude_springs)
        deformable_object->Add_Altitude_Springs(deformable_object->triangulated_surface->triangle_mesh,
            altitude_springs_stiffness,altitude_springs_overdamping_fraction,
            altitude_springs_use_shortest_only,altitude_springs_use_compressed_by_threshold_only,
            altitude_springs_fraction_compression,altitude_springs_limit_time_step_by_strain_rate,
            altitude_springs_max_strain_per_time_step,altitude_springs_use_rest_state_for_strain_rate,
            altitude_springs_restlength_enlargement_fraction);

    // Adjust particles
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles=deformable_object->triangulated_surface->particles;
#if 0
    MATRIX<T,4> rotation=MATRIX<T,4>::Rotation_Matrix_Z_Axis(-2*pi/3);
    for(int i=0;i<cloth_grid.m/2;i++) for(int j=0;j<cloth_grid.n;j++){
        int node=i+cloth_grid.m*(j-1);particles.X(node)=rotation*particles.X(node);}
#endif
    particles.X.array+=VECTOR_3D<T>(0,cloth_height,0);

    enslaved_nodes.Append(Node_Index(1,1));
//    enslaved_nodes.Append(Node_Index(1,cloth_grid.n));
#if 0
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<T>("../../Public_Data/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list(id)->is_static=true;
    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    BOX_2D<T> enslaved_square(-enslaved_square_size/2,enslaved_square_size/2,-enslaved_square_size/2,enslaved_square_size/2);
    for(int i=0;i<cloth_grid.m;i++) for(int j=0;j<cloth_grid.n;j++){
        if(enslaved_square.Lazy_Inside(cloth_grid.X(i,j))){
            int node=i+cloth_grid.m*(j-1);
            enslaved_nodes.Append(node);
            //particles.mass(node)=FLT_MAX;
        }}
#endif

    // static cloth for ground
    id=Add_Deformable_Object(GRID_2D<T>(2,2,-1,1,-1,1),MATRIX<T,4>::Rotation_Matrix_Y_Axis(pi/2)*MATRIX<T,4>::Rotation_Matrix_X_Axis(-pi/2));
    deformable_object=solids_parameters.deformable_body_parameters.list.deformable_objects(id);
    //deformable_object->triangulated_surface->Set_Density(1e10);
    deformable_object->triangulated_surface->Set_Density(100000);
    deformable_object->triangulated_surface->Set_Mass_Of_Particles();

    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Initialize_Bodies();
}
//#####################################################################
// Function Add_Circle_Deformable_Object
//#####################################################################
int Add_Deformable_Object(const GRID<TV>& cloth_grid,const MATRIX<T,4>& transform)
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Triangulated_Surface();
    TRIANGULATED_SURFACE<T>& triangulated_surface=*solids_parameters.deformable_body_parameters.list(index).triangulated_surface;
    TRIANGLE_MESH& triangle_mesh=triangulated_surface.triangle_mesh;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles=triangulated_surface.particles;

    triangle_mesh.Initialize_Herring_Bone_Mesh(cloth_grid.m,cloth_grid.n);
    particles.Add_Elements(triangle_mesh.number_nodes);
    for(int i=0;i<cloth_grid.m;i++) for(int j=0;j<cloth_grid.n;j++){
        int node=i+cloth_grid.m*(j-1);particles.X(node)=transform*VECTOR_3D<T>(cloth_grid.X(i,j));particles.V(node)=VECTOR_3D<T>();}

    return index;
}
//#####################################################################
// Function Solids_Example_Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY<VECTOR_3D<T> >& V,const T time)
{
    Zero_Out_Enslaved_Velocity_Nodes(V,time,id_number);
}
//#####################################################################
// Function Solids_Example_Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> >& V,const T time)
{
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles=solids_parameters.deformable_body_parameters.list(1).triangulated_surface->particles;
//    if(time<0)
    {for(int i=0;i<enslaved_nodes.m;i++) V(enslaved_nodes(i))=VECTOR_3D<T>();}
//    else{for(int i=0;i<enslaved_nodes.m;i++) V(enslaved_nodes(i))=VECTOR_3D<T>::Cross_Product(particles.X(enslaved_nodes(i)),VECTOR_3D<T>(0,1,0));}
}
//#####################################################################
};
}
#endif

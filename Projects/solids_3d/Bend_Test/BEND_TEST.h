//#####################################################################
// Copyright 2003, 2004, Robert Bridson, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BEND_TEST
//##################################################################### 
//
// Triangular cloth draped over a ball (for convergence tests)
//
//#####################################################################
#ifndef __BEND_TEST__
#define __BEND_TEST__

#include "../CLOTH_EXAMPLE.h"
namespace PhysBAM{

template<class T>
class BEND_TEST:public CLOTH_EXAMPLE<T>
{
public:
    int number_side_panels;
    T aspect_ratio,side_length;
    T density;

    BEND_TEST()
    {
        number_side_panels=20;
        aspect_ratio=1;
        side_length=10;

        final_time=(T)4;
        frame_rate=30;
        restart_step_number=0;
        output_directory="Bend_Test/output";
        output_prefix="cloth_";
        check_initial_mesh_for_self_intersection=true;
        print_residuals=false;
        verbose_dt=true;

        // --- begin auto generated portion

        // Not using fvm, so Rayleigh_coefficient, poissons_ratio, youngs_modulus are irrelevant
        use_masses_and_springs=true;
        use_fvm=false;
        use_diagonalized_fvm=false;

        allow_intersections=true;   // always on (no knob)
        allow_intersections_tolerance=(T)1e-8;  // no knob
        bending_cutoff_fraction_of_minimum_area=0;  // not used at ILM
        limit_time_step_by_strain_rate=true;    // always on (no knob)
        edge_spring_restlength_enlargement_fraction=0;
        altitude_spring_restlength_enlargement_fraction=(T).01;
        limit_time_step_by_strain_rate=true;    // always on!

        // --- begin values set in sim file

        cg_tolerance=(T)0.01; // cgTolerance
        cg_iterations=900000; // cgIterations
        cfl_number=(T)1; // cflMult

        gravity=(T)9.8; // gravityMult
        density=(T)0.4; // skinDensity
        dynamic_ether_viscosity=(T)0; // dynamicEtherViscosity

        edge_spring_stiffness=(T)40; // edgeSpringConstant
        edge_spring_overdamping_fraction=(T)1; // edgeSpringOverdampingFraction
        max_strain_per_time_step=(T)0.1; // tensileSpringCfl, triangleAltitudeSpringCfl, triangleBendingSpringCfl

        perform_self_collision=false; // doSelfCollision
        collisions_small_number=(T)1e-08; // collisionsSmallNumber
        collisions_repulsion_thickness=(T)0.001; // collisionsRepulsionThickness
        collisions_collision_thickness=(T)1e-06; // collisionsCollisionThickness
        collisions_repulsion_clamp_fraction=(T)0.4; // collisionsRepulsionClampFraction
        spring_limiter_fraction=(T)0.1; // springLimiterFraction
        self_collision_friction_coefficient=(T)0.4; // selfCollisionFrictionCoefficient
        collisions_nonrigid_collision_attempts=10; // collisionsNonrigidCollisionAttempts
        collision_repulsion_spring_multiplier=100; // collisionRepulsionSpringMultiplier

        use_altitude_springs=true; // useTriangleAltitudeSprings
        use_shortest_altitude_spring_only=true; // triangleAltitudeShortestSpringOnly
        use_altitude_springs_compressed_by_threshold_only=true; // triangleAltitudeCompressByThresholdOnly
        altitude_spring_fraction_compression=(T)0.1; // triangleAltitudeFractionCompression
        altitude_spring_stiffness=(T)4; // triangleAltitudeStiffness
        altitude_spring_overdamping_fraction=2; // triangleAltitudeOverdampingFraction
        // triangleAltitudeSpringCfl=0.1 (see max_strain_per_time_step above)

        use_bending_elements=true; // useTriangleBendingElements
        bending_stiffness=(T)10; // triangleBendingElementsStiffness
        bending_damping=(T)1; // triangleBendingElementsDamping
        max_bend_strain_per_time_step=(T)0.01; // triangleBendingElementsCfl
        bending_cutoff_fraction_of_triangles=(T)0; // triangleBendingElementsFractionOfElementsIgnored

        use_bending_springs=false; // useTriangleBendingSprings
        bending_spring_stiffness=40; // triangleBendingSpringStiffness
        bending_spring_overdamping_fraction=2; // triangleBendingSpringDamping
        // triangleBendingSpringCfl=0.5 (see max_strain_per_time_step above)
        bending_spring_restlength_enlargement_fraction=(T)0.01; // triangleBendingSpringRestlengthEnlargementFraction

        // --- end auto generated portion
    }

    ~BEND_TEST()
    {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data(TRIANGULATED_SURFACE<T>& triangulated_surface)
{
    TRIANGLE_MESH& triangle_mesh=triangulated_surface.triangle_mesh;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles=triangulated_surface.particles;

    particles.Update_Velocity();particles.Store_Mass();

    int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
    triangle_mesh.Initialize_Herring_Bone_Mesh(m,n);
    for(int k=0;k<triangle_mesh.number_nodes;k++) particles.array_collection->Add_Element();
    T dx=aspect_ratio*side_length/(m-1),dy=side_length/(n-1);
    for(int i=0;i<m;i++) for(int j=0;j<n;j++){int node=i+m*(j-1);
        particles.X(node)=VECTOR_3D<T>((i-1)*dx,.5,(j-1)*dy);
        particles.V(node)=VECTOR_3D<T>(0,0,0);}

    triangulated_surface.Set_Density(density);
    triangulated_surface.Set_Mass_Of_Particles(false);

    LOG::SCOPE scope("mesh statistics");
    triangulated_surface.Print_Statistics(LOG::cout);
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
// for external forces and velocities
void Set_External_Velocities(ARRAY<VECTOR_3D<T> >& V,const T time)
{
    Zero_Out_Enslaved_Velocity_Nodes(V,time);
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
// for external forces and velocities
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> >& V,const T time)
{
    int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
    for(int i=0;i<m/2+1;i++) for(int j=0;j<n;j++) V(i+m*(j-1))=VECTOR_3D<T>(0,0,0);
}
//#####################################################################
};
}
#endif


//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RHOMBIC_DODECAHEDRON
//#####################################################################
#ifndef __RHOMBIC_DODECAHEDRON__
#define __RHOMBIC_DODECAHEDRON__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Geometry/Red_Green/RED_GREEN_GRID_3D.h>
#include <Solids_And_Fluids/MELTING_EXAMPLE_3D.h>
namespace PhysBAM{

template <class T,class RW>
class RHOMBIC_DODECAHEDRON:public MELTING_EXAMPLE_3D<T,RW>
{
public:
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity;
    VECTOR_3D<T> initial_angular_velocity;
    T side_length;
    int m,n,mn;
    
    RHOMBIC_DODECAHEDRON()
        :initial_height(0),initial_orientation(),initial_velocity(),initial_angular_velocity(1,1,0),
        side_length(2),m(1),n(1),mn(1)
    {
        frame_rate=24*5;
        last_frame=10*(int)frame_rate;
        restart=false;restart_frame=205;
        solids_parameters.cfl=(T).5;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        melting_parameters.maximum_depth=6;

        melting_parameters.refine_for_high_deformation=true;
        melting_parameters.refine_boundary_only=true;
        melting_parameters.successive_refinement_multiplier=(T)1.5;
        melting_parameters.expansive_refinement_threshold=(T).05;
        melting_parameters.compressive_refinement_threshold=(T).05;

        melting_parameters.write_overlay_levelset=false;
        solids_parameters.perform_self_collision=false;
        
        melting_parameters.number_of_objects=1;
        output_directory="Rhombic_Dodecahedron/output";
        solids_parameters.collide_with_interior=false;
        verbose_dt=true;
    }

    ~RHOMBIC_DODECAHEDRON()
    {}

//#####################################################################
// Function Initialize_Particle_Positions_And_Velocities
//#####################################################################
void Initialize_Grid(const int object,RED_GREEN_GRID_3D<T>& grid)
{
    assert(object==1);T y=2;
    grid.Initialize(GRID_3D<T>(4*m+1,4*n+1,4*mn+1,0,side_length,y,y+side_length*n/m,0,side_length*mn/m),melting_parameters.maximum_depth);
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi(const int object,ARRAY<T>& phi)
{
    ARRAY<T>::copy(-1,phi);
}
//#####################################################################
// Function Initialize_Levelset_Velocity
//#####################################################################
void Initialize_Levelset_Velocity(const int object,ARRAY<VECTOR_3D<T> >& V)
{
    ARRAY<VECTOR_3D<T> >::copy(VECTOR_3D<T>(),V);
}
//#####################################################################
// Function Initialize_Particle_Positions_And_Velocities
//#####################################################################
void Initialize_Particle_Positions_And_Velocities(const int object)
{
    DEFORMABLE_OBJECT_3D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(1);
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_object.embedded_tetrahedralized_volume->tetrahedralized_volume;
    PARTICLES<T,VECTOR_3D<T> >& particles=deformable_object.particles;

    particles.Update_Velocity();
    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());
    for(int i=0;i<particles.array_collection->Size();i++){
        particles.X(i)=center+initial_orientation.Rotate(particles.X(i)-center);
        particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,particles.X(i)-center);
        particles.X(i).y+=initial_height;}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Forces()
{
    DEFORMABLE_OBJECT_3D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(1);
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_object.embedded_tetrahedralized_volume->tetrahedralized_volume;

    deformable_object.particles.Store_Mass();
    tetrahedralized_volume.Set_Density(1000);
    tetrahedralized_volume.Set_Mass_Of_Particles(false);

    deformable_object.Delete_Forces();
    solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_object.particles,solids_parameters.gravity,solids_parameters.gravity_direction));
    solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(tetrahedralized_volume,new DIAGONALIZED_SPLINE_MODEL_3D<T>(tetrahedralized_volume,(T)5e4,(T).3,(T).5,7,(T).01)));
}
//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
void Initialize_Rigid_Bodies()
{
    int index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>("../../Public_Data/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list(index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list(index)->is_static=true;

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
};
}
#endif


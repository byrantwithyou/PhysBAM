//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MELTING_SPHERE
//#####################################################################
#ifndef __MELTING_SPHERE__
#define __MELTING_SPHERE__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <Level_Sets/LEVELSET_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Red_Green/RED_GREEN_GRID_3D.h>
#include <Solids_And_Fluids/MELTING_EXAMPLE_3D.h>
namespace PhysBAM{

template <class T,class RW>
class MELTING_SPHERE:public MELTING_EXAMPLE_3D<T,RW,SOLIDS_FLUIDS_EXAMPLE_3D<RW> >
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::data_directory;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::verbose_dt;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::solids_parameters;
    using MELTING_EXAMPLE_3D<T,RW,SOLIDS_FLUIDS_EXAMPLE_3D<RW> >::melting_parameters;

    SPHERE<T> sphere;
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity;
    VECTOR_3D<T> initial_angular_velocity;
    T side_length;
    int m,n,mn;
    
    MELTING_SPHERE()
        :MELTING_EXAMPLE_3D<T,RW,SOLIDS_FLUIDS_EXAMPLE_3D<RW> >(fluids_parameters.NONE),initial_height(0),initial_orientation(),initial_velocity(),initial_angular_velocity(),
        side_length(2),m(2),n(2),mn(2)
    {
        frame_rate=24*10;
        last_frame=10*(int)frame_rate;
        restart=false;restart_frame=205;
        solids_parameters.cfl=(T).5;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        melting_parameters.maximum_depth=3;

        melting_parameters.refine_near_interface=true;
        //successive_refinement_multiplier=(T)1.5;
        //expansive_refinement_threshold=(T).4;
        //compressive_refinement_threshold=(T).4;

        melting_parameters.write_overlay_levelset=false;
        solids_parameters.perform_self_collision=false;
        
        output_directory="Melting_Sphere/output";
        solids_parameters.collide_with_interior=false;
        verbose_dt=true;
    }

    ~MELTING_SPHERE()
    {}

//#####################################################################
// Function Initialize_Deformable_And_Rigid_Bodies
//#####################################################################
void Initialize_Deformable_And_Rigid_Bodies()
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Embedded_Tetrahedralized_Volume();
    Add_Melting_Object(melting_parameters.DEFORMABLE,index);

    index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>("../../Public_Data/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->is_static=true;

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Initialize_Particle_Positions_And_Velocities
//#####################################################################
void Initialize_Grid(const int object,RED_GREEN_GRID_3D<T>& grid)
{
    assert(object==1);T y=2;
    grid.Initialize(GRID_3D<T>(4*m+1,4*n+1,4*mn+1,0,side_length,y,y+side_length*n/m,0,side_length*mn/m),melting_parameters.maximum_depth);
    sphere.center=grid.uniform_grid.Domain().Center();
    sphere.radius=.9;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi(const int object,ARRAY<T>& phi)
{
    assert(object==1);
    RED_GREEN_GRID_3D<T>& grid=melting_parameters.levelsets(1)->grid;
    ARRAY<VECTOR_3D<T> >& node_locations=grid.Node_Locations();
    
    for(int p=0;p<phi.m;p++) phi(p)=sphere.Signed_Distance(node_locations(p));    
    //ARRAY<T>::copy(-1,phi);
}
//#####################################################################
// Function Initialize_Levelset_Velocity
//#####################################################################
void Initialize_Levelset_Velocity(const int object,ARRAY<VECTOR_3D<T> >& V)
{
    assert(object==1);
    RED_GREEN_GRID_3D<T>& grid=melting_parameters.levelsets(1)->grid;
    ARRAY<VECTOR_3D<T> >& node_locations=grid.Node_Locations();
    
    //for(int p=0;p<V.m;p++) V(p)=VECTOR_3D<T>(node_locations(p).y-sphere.center.y,0,0);
    for(int p=0;p<V.m;p++) V(p)=VECTOR_3D<T>(3*max(1-2*fabs(node_locations(p).y-sphere.center.y),0.),0,0);
}
//#####################################################################
// Function Initialize_Particle_Positions_And_Velocities
//#####################################################################
void Initialize_Particle_Positions_And_Velocities(const int object)
{
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=solids_parameters.deformable_body_parameters.list(object).embedded_tetrahedralized_volume->tetrahedralized_volume;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles=solids_parameters.deformable_body_parameters.list(object).particles;

    particles.Update_Velocity();
    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());
    for(int i=0;i<particles.array_collection->Size();i++){
        particles.X(i)=center+initial_orientation.Rotate(particles.X(i)-center);
        particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,particles.X(i)-center);
        particles.X(i).y+=initial_height;}
}
//#####################################################################
// Function Initialize_Forces
//#####################################################################
void Initialize_Forces()
{
    DEFORMABLE_OBJECT_3D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(1);
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_object.embedded_tetrahedralized_volume->tetrahedralized_volume;

    solids_parameters.deformable_body_parameters.list(1).particles.Store_Mass();
    tetrahedralized_volume.Set_Density(1000);
    tetrahedralized_volume.Set_Mass_Of_Particles(false);

    solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_object.particles));
    solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(tetrahedralized_volume,new DIAGONALIZED_SPLINE_MODEL_3D<T>(tetrahedralized_volume,(T)5e4,(T).3,(T).5,7,(T).01)));
}
//#####################################################################
// Function Postprocess_Solids_Substep
//#####################################################################
void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE
{
    Update_Solids_Topology_For_Melting(1/frame_rate,time); // ridiculous hack to test if things work
}
//#####################################################################
};
}
#endif


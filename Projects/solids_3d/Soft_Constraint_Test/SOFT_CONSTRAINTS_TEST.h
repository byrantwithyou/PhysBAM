//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOFT_CONSTRAINTS_TEST
//#####################################################################
#ifndef __SOFT_CONSTRAINTS_TEST__
#define __SOFT_CONSTRAINTS_TEST__

#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <Forces_And_Torques/CONSTRAINED_POINT_FIXED.h>
#include <Forces_And_Torques/CONSTRAINED_POINT_IN_RIGID_BODY_3D.h>
#include <Forces_And_Torques/CONSTRAINED_POINT_IN_TETRAHEDRON.h>
#include <Forces_And_Torques/CONSTRAINED_POINT_IN_TRIANGLE.h>
#include <Forces_And_Torques/SOFT_CONSTRAINTS_3D.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>
namespace PhysBAM{

template <class T,class RW>
class SOFT_CONSTRAINTS_TEST:public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart_frame;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::data_directory;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::solids_parameters;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::verbose_dt;

    DEFORMABLE_OBJECT_LIST_3D<T> soft_constraints_body_list;
    SOFT_CONSTRAINTS_3D<T> soft_constraints;

    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    std::string input_file;
    ARRAY<ARRAY<int> > rigid_body_constrained_points;
    ARRAY<ARRAY<int> > deformable_object_constrained_points;

    SOFT_CONSTRAINTS_TEST()
        :SOLIDS_FLUIDS_EXAMPLE_3D<RW>(fluids_parameters.NONE),initial_height(1),initial_orientation(0,VECTOR_3D<T>(1,0,0)),initial_velocity(0,0,0),initial_angular_velocity(0,0,0)
    {
        solids_parameters.collisions_repulsion_thickness=(T)1e-2;
        solids_parameters.collisions_repulsion_clamp_fraction=(T).9;
        solids_parameters.collision_repulsion_spring_multiplier=100;

        last_frame=50*24;
        //last_frame=17;//10*24;
        restart=false;restart_frame=16;   
        solids_parameters.cfl=(T)10;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        output_directory="Soft_Constraint_Test/output";
        input_file=data_directory+"/Tetrahedralized_Volumes/adaptive_torus_float.tet";
        solids_parameters.perform_self_collision=false;
        solids_parameters.collide_with_interior=false;
        solids_parameters.synchronize_multiple_objects=false;

        verbose_dt=true;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    }

    ~SOFT_CONSTRAINTS_TEST()
    {}

    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return false;}

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    RIGID_BODY_LIST_3D<T>& rigid_body_list=solids_parameters.rigid_body_parameters.list;

    int id;

    GRID<TV> cloth_grid(5,5,0,1,0,1);
    MATRIX<T,4> transform=MATRIX<T,4>::Rotation_Matrix_Z_Axis(-3*pi/4);

#if 0
    id=solids_parameters.deformable_body_parameters.list.Add_Deformable_Triangulated_Surface();
    TRIANGULATED_SURFACE<T>* triangulated_surface=solids_parameters.deformable_body_parameters.list(id).triangulated_surface;
    PARTICLES<T,VECTOR_3D<T> >* particles=&triangulated_surface->particles;
    triangulated_surface->triangle_mesh.Initialize_Herring_Bone_Mesh(cloth_grid.m,cloth_grid.n);
    particles->array_collection->Add_Elements(triangulated_surface->triangle_mesh.number_nodes);
    for(int i=0;i<cloth_grid.m;i++) for(int j=0;j<cloth_grid.n;j++){
        int node=i+cloth_grid.m*(j-1);particles->X(node)=transform*VECTOR_3D<T>(cloth_grid.X(i,j));particles->V(node)=VECTOR_3D<T>();}
    triangulated_surface->Set_Density(1);
    triangulated_surface->Set_Mass_Of_Particles(true);
    solids_parameters.deformable_body_parameters.list(id).Add_Body_Forces(*triangulated_surface);
    solids_parameters.deformable_body_parameters.list(id).Add_Edge_Springs(triangulated_surface->triangle_mesh,200/(1+sqrt((T)2)),2);

    transform=MATRIX<T,4>::Translation_Matrix(VECTOR_3D<T>(0,-2,0))*MATRIX<T,4>::Rotation_Matrix_Z_Axis(-3*pi/4);
    id=solids_parameters.deformable_body_parameters.list.Add_Deformable_Triangulated_Surface();
    triangulated_surface=solids_parameters.deformable_body_parameters.list(id).triangulated_surface;
    particles=&triangulated_surface->particles;
    triangulated_surface->triangle_mesh.Initialize_Herring_Bone_Mesh(cloth_grid.m,cloth_grid.n);
    particles->array_collection->Add_Elements(triangulated_surface->triangle_mesh.number_nodes);
    for(int i=0;i<cloth_grid.m;i++) for(int j=0;j<cloth_grid.n;j++){
        int node=i+cloth_grid.m*(j-1);particles->X(node)=transform*VECTOR_3D<T>(cloth_grid.X(i,j));particles->V(node)=VECTOR_3D<T>();}
    triangulated_surface->Set_Density(1);
    triangulated_surface->Set_Mass_Of_Particles(true);
    solids_parameters.deformable_body_parameters.list(id).Add_Body_Forces(*triangulated_surface);
    solids_parameters.deformable_body_parameters.list(id).Add_Edge_Springs(triangulated_surface->triangle_mesh,200/(1+sqrt((T)2)),2);
#endif

#if 0
    id=solids_parameters.deformable_body_parameters.list.Add_Deformable_Tetrahedralized_Volume();
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*solids_parameters.deformable_body_parameters.list(id).tetrahedralized_volume;
    TETRAHEDRON_MESH& tetrahedron_mesh=tetrahedralized_volume.tetrahedron_mesh;
    PARTICLES<T,VECTOR_3D<T> >& particles=tetrahedralized_volume.particles;

    std::istream* input=FILE_UTILITIES::Safe_Open_Input(input_file);tetrahedralized_volume.template Read<RW>(*input);delete input;
    std::cout << "total vertices = " << particles.array_collection->Size() << std::endl;std::cout << "total tetrahedra = " << tetrahedron_mesh.tetrahedrons.m << std::endl;
    particles.Update_Velocity();particles.Store_Mass(); // in case they're not stored in the file!

    tetrahedralized_volume.tetrahedron_mesh.Initialize_Incident_Tetrahedrons();

    tetrahedralized_volume.Set_Density(100);tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);
    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=1;i<=particles.array_collection->Size();i++){
        particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,particles.X(i)-center);
        particles.X(i)=center+initial_orientation.Rotate(particles.X(i)-center);
        particles.X(i).y+=initial_height-bottom;}

    solids_parameters.deformable_body_parameters.list(id).Add_Body_Forces(tetrahedralized_volume);
    solids_parameters.deformable_body_parameters.list(id).Add_Diagonalized_Neo_Hookean_Elasticity(tetrahedralized_volume,(T)2e5,(T).45,(T).01);

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/subdivided_box");
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list(id)->position=VECTOR_3D<T>(-5,14,0);
    solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction);

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/subdivided_box");
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list(id)->position=VECTOR_3D<T>(-10,15,0);
    solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction);

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/subdivided_box",2);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list(id)->position=VECTOR_3D<T>(-8,10,0);
    solids_parameters.rigid_body_parameters.list(id)->is_static=true;
#endif
#if 0
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/subdivided_box");
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list(id)->position=VECTOR_3D<T>(0,25,0);
    solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction);
    solids_parameters.rigid_body_parameters.list(id)->is_static=true;
#endif
#if 1
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/box");
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list(id)->position=VECTOR_3D<T>(0,6,0);
//    solids_parameters.rigid_body_parameters.list(id)->orientation=QUATERNION<T>(.4,VECTOR_3D<T>(0,1,0));
    solids_parameters.rigid_body_parameters.list(id)->is_kinematic=true;
    solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction);
#endif
#if 1
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/box");
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list(id)->position=VECTOR_3D<T>(3,6,0);
//    solids_parameters.rigid_body_parameters.list(id)->orientation=QUATERNION<T>(-.4,VECTOR_3D<T>(0,1,0));
    solids_parameters.rigid_body_parameters.list(id)->is_kinematic=true;
    solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction);

//    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/box");
//    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction((T).3);
//    solids_parameters.rigid_body_parameters.list(id)->position=VECTOR_3D<T>(6,6,0);
//    solids_parameters.rigid_body_parameters.list(id)->orientation=QUATERNION<T>(-.4,VECTOR_3D<T>(0,1,0));
//    solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction);

    //id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/ground");
    //rigid_body_particles.Rigid_Body(id).is_static=true;
    //rigid_body_particles.Rigid_Body(id).add_to_spatial_partition=false;
#endif

#if 0
    ARRAY<int> segments(2,3);
    segments.Set(1,1,2);
    segments.Set(2,2,3);
    segments.Set(3,3,4);

    soft_constraints.Add_Constrained_Point(new CONSTRAINED_POINT_FIXED<T,VECTOR_3D<T> >(VECTOR_3D<T>(0,15,0)));
    soft_constraints.Add_Constrained_Point(new CONSTRAINED_POINT_IN_RIGID_BODY<TV>(*solids_parameters.rigid_body_parameters.list(1),VECTOR_3D<T>()));
    soft_constraints.Add_Constrained_Point(new CONSTRAINED_POINT_IN_RIGID_BODY<TV>(*solids_parameters.rigid_body_parameters.list(2),VECTOR_3D<T>()));
    soft_constraints.Add_Constrained_Point(new CONSTRAINED_POINT_IN_TETRAHEDRON<T,VECTOR_3D<T> >(tetrahedralized_volume.tetrahedron_mesh,particles,3820,VECTOR_3D<T>(0,0,1)));
#endif
#if 0
    ARRAY<int> segments(2,1);
    segments.Set(1,1,2);

    soft_constraints.Add_Constrained_Point(new CONSTRAINED_POINT_IN_RIGID_BODY<TV>(*solids_parameters.rigid_body_parameters.list(1),VECTOR_3D<T>()));
    soft_constraints.Add_Constrained_Point(new CONSTRAINED_POINT_IN_RIGID_BODY<TV>(*solids_parameters.rigid_body_parameters.list(2),VECTOR_3D<T>()));
#endif

    rigid_body_constrained_points.Resize(solids_parameters.rigid_body_parameters.list.rigid_bodies.m);
    deformable_object_constrained_points.Resize(solids_parameters.deformable_body_parameters.list.deformable_objects.m);

    int point_id;
    ARRAY<int> segments(2,0);

#if 0
    point_id=soft_constraints.Add_Constrained_Point(new CONSTRAINED_POINT_FIXED<T,VECTOR_3D<T> >(VECTOR_3D<T>(0,5,0)));

#if 1
    point_id=soft_constraints.Add_Constrained_Point(new CONSTRAINED_POINT_IN_TRIANGLE<T,VECTOR_3D<T> >(
        solids_parameters.deformable_body_parameters.list(1).triangulated_surface->triangle_mesh,
        solids_parameters.deformable_body_parameters.list(1).triangulated_surface->particles,
        1,VECTOR_3D<T>(1,0,0)));
    deformable_object_constrained_points(1).Append(point_id);
    point_id=soft_constraints.Add_Constrained_Point(new CONSTRAINED_POINT_IN_TRIANGLE<T,VECTOR_3D<T> >(
        solids_parameters.deformable_body_parameters.list(1).triangulated_surface->triangle_mesh,
        solids_parameters.deformable_body_parameters.list(1).triangulated_surface->particles,
        31,VECTOR_3D<T>(0,1,0)));
    deformable_object_constrained_points(1).Append(point_id);
    point_id=soft_constraints.Add_Constrained_Point(new CONSTRAINED_POINT_IN_TRIANGLE<T,VECTOR_3D<T> >(
        solids_parameters.deformable_body_parameters.list(2).triangulated_surface->triangle_mesh,
        solids_parameters.deformable_body_parameters.list(2).triangulated_surface->particles,
        1,VECTOR_3D<T>(1,0,0)));
    deformable_object_constrained_points(2).Append(point_id);
    segments.Append(1,2);
    segments.Append(3,4);
#endif
#endif

#if 1
    VECTOR_3D<T> world_pos(1.5,10,0);
    point_id=soft_constraints.Add_Constrained_Point(new CONSTRAINED_POINT_IN_RIGID_BODY<TV>(*rigid_body_list(1),VECTOR_3D<T>(0,0,0)));
    rigid_body_constrained_points(1).Append(point_id);
    point_id=soft_constraints.Add_Constrained_Point(new CONSTRAINED_POINT_IN_RIGID_BODY<TV>(*rigid_body_list(2),VECTOR_3D<T>(0,0,0)));
    rigid_body_constrained_points(2).Append(point_id);
//    point_id=soft_constraints.Add_Constrained_Point(new CONSTRAINED_POINT_IN_RIGID_BODY<TV>(*rigid_body_list(3),VECTOR_3D<T>(0,0,0)));
//    rigid_body_constrained_points(3).Append(point_id);
    segments.Append(1,2);
//    segments.Append(2,3);
#endif

#if 0
    point_id=soft_constraints.Add_Constrained_Point(new CONSTRAINED_POINT_IN_RIGID_BODY<TV>(*solids_parameters.rigid_body_parameters.list(1),VECTOR_3D<T>()));
    rigid_body_constrained_points(1).Append(point_id);

    segments.Append(1,point_id);
#endif

    soft_constraints.deformable_object.segmented_curve->segment_mesh.Initialize_Segment_Mesh(segments);
    soft_constraints.Synchronize_Particles_With_Constrained_Points();
    soft_constraints.deformable_object.Add_Edge_Springs(soft_constraints.deformable_object.segmented_curve->segment_mesh,2000,1);
//    soft_constraints.deformable_object.Add_Edge_Springs(soft_constraints.deformable_object.segmented_curve->segment_mesh,20,1);
//    soft_constraints.deformable_object.linear_springs(1)->Clamp_Restlength(0.1);

    soft_constraints_body_list.Add_Deformable_Object(&soft_constraints.deformable_object);

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Initialize_Bodies();
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(id==int(1)){
        INTERPOLATION_CURVE<T,VECTOR_3D<T> > motion_curve;
        motion_curve.Add_Control_Point(0,VECTOR_3D<T>(0,6,0));
        motion_curve.Add_Control_Point(2,VECTOR_3D<T>(0,10,0));
        frame.t=motion_curve.Value(time);}
    else if(id==int(2)){
        INTERPOLATION_CURVE<T,VECTOR_3D<T> > motion_curve;
        motion_curve.Add_Control_Point(0,VECTOR_3D<T>(3,6,0));
        motion_curve.Add_Control_Point(2,VECTOR_3D<T>(3.4,10,0));
        frame.t=motion_curve.Value(time);}
}
//#####################################################################
// Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
    if(time>1){
        solids_parameters.rigid_body_parameters.list(1)->is_kinematic=false;
        solids_parameters.rigid_body_parameters.list(2)->is_kinematic=false;
    }
}
//#####################################################################
// Apply_Constraints
//#####################################################################
void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE
{
    std::cout << "Applying constraint " << std::endl;
    soft_constraints.Apply_Constraint_Correcting_Impulses(dt,time);
}
//#####################################################################
// Constraints_CFL
//#####################################################################
T Constraints_CFL()
{
    T cfl=soft_constraints.CFL();
    std::cout << "Got constraints cfl: " << cfl << std::endl;
    return cfl;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE
{
    std::string prefix=output_directory+"/soft_constraints_";
    if(frame==first_frame) soft_constraints_body_list.template Write_Static_Variables<RW>(prefix);
    soft_constraints_body_list.template Write_Dynamic_Variables<RW>(prefix,frame);
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Write_Output_Files(frame);
}
//#####################################################################
};
}
#endif

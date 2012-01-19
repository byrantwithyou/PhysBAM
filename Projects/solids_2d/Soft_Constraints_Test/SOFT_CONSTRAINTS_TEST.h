//#####################################################################
// Copyright 2004-2006, Eran Guendelman, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOFT_CONSTRAINTS_TEST
//#####################################################################
#ifndef __SOFT_CONSTRAINTS_TEST__
#define __SOFT_CONSTRAINTS_TEST__

#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Bindings/RIGID_BODY_BINDING.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Constitutive_Models/DIAGONALIZED_SPLINE_MODEL_2D.h>
#include <Deformable_Objects/KINEMATIC_BINDING.h>
namespace PhysBAM{

template <class T,class RW>
class SOFT_CONSTRAINTS_TEST:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    typedef VECTOR<T,2> TV;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::last_frame;using BASE::output_directory;using BASE::data_directory;using BASE::verbose_dt;

    T initial_height;
    T initial_orientation;
    VECTOR<T,2> initial_velocity;
    T initial_angular_velocity;
    GRID<TV> mattress_grid;
    bool fix_boundary;
    bool use_ground;

    SOFT_CONSTRAINTS_TEST()
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>(0,fluids_parameters.NONE),
        initial_height(4),initial_orientation(0),initial_velocity(0,0),initial_angular_velocity(0),
        mattress_grid(3,5,-.25,.25,-.5,.5),fix_boundary(false),use_ground(false)
    {
        last_frame=100*24;
        solids_parameters.cfl=(T).5;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        output_directory="Soft_Constraints_Test/output";
        solids_parameters.collide_with_interior=false;
        verbose_dt=true;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    }

    ~SOFT_CONSTRAINTS_TEST()
    {}

//#####################################################################
// Function Get_Initial_Data()
//#####################################################################
void Get_Initial_Data()
{
    DEFORMABLE_OBJECT<TV>& deformable_object=solid_body_collection.deformable_object;

    // triangulated area
    TRIANGULATED_AREA<T>& triangulated_area=*TRIANGULATED_AREA<T>::Create();
    PARTICLES<TV>& area_particles=triangulated_area.particles;
    triangulated_area.Initialize_Square_Mesh_And_Particles(mattress_grid);
    area_particles.Update_Velocity();area_particles.Store_Mass();
    triangulated_area.Set_Density(10);triangulated_area.Set_Mass_Of_Particles(false);
    triangulated_area.Update_Bounding_Box();
    VECTOR<T,2> center(triangulated_area.bounding_box->Center());T bottom=triangulated_area.bounding_box->ymin;
    for(int i=0;i<area_particles.array_collection->Size();i++){
        area_particles.X(i)=center+MATRIX<T,2>::Rotation_Matrix(initial_orientation)*(area_particles.X(i)-center);
        VECTOR<T,2> radial=area_particles.X(i)-center;
        T temp=radial.y;radial.y=-radial.x;radial.x=temp;
        area_particles.V(i)=initial_velocity+initial_angular_velocity*(radial);
        area_particles.X(i).y+=initial_height-bottom;}
    triangulated_area.Check_Signed_Areas_And_Make_Consistent(true);
    deformable_object.Add_Structure(triangulated_area.Append_Particles_And_Create_Copy(deformable_object.particles));
    deformable_object.collisions.collision_structures.Append(deformable_object.structures.Last());

    // segmented curve
    SEGMENTED_CURVE_2D<T>& segmented_curve=*SEGMENTED_CURVE_2D<T>::Create();
    PARTICLES<TV>& curve_particles=segmented_curve.particles;
    segmented_curve.mesh.Initialize_Straight_Mesh(4);curve_particles.array_collection->Add_Elements(segmented_curve.mesh.number_nodes);
    curve_particles.Update_Velocity();curve_particles.Store_Mass();
    for(int i=0;i<4;i++) segmented_curve.particles.X(i)=VECTOR<T,2>(i,1); // TODO: remove
    deformable_object.Add_Structure(segmented_curve.Append_Particles_And_Create_Copy(deformable_object.particles));

    // correct number nodes
    for(int i=0;i<deformable_object.structures.m;i++) deformable_object.structures(i)->Update_Number_Nodes();

    // rigid bodies
    int index;
    index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/square");
    solids_parameters.rigid_body_parameters.list(index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list(index)->frame.t=VECTOR<T,2>(-5,14);
    //solids_parameters.rigid_body_parameters.list(index)->velocity=VECTOR<T,2>(10,0);
    solids_parameters.rigid_body_parameters.list(index)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction);

    index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/square");
    solids_parameters.rigid_body_parameters.list(index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list(index)->frame.t=VECTOR<T,2>(-10,15);
    //solids_parameters.rigid_body_parameters.list(index)->velocity=VECTOR<T,2>(10,0);
    solids_parameters.rigid_body_parameters.list(index)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction);

    if(use_ground){
        index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>("../../Public_Data/Rigid_Bodies_2D/ground");
        solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T).3);
        solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->is_static=true;}

    // constraints and bindings
    int curve_offset=triangulated_area.mesh.number_nodes;
    solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new KINEMATIC_BINDING<TV>(deformable_object.particles,curve_offset+1,VECTOR<T,2>(0,15)));
    solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new RIGID_BODY_BINDING<TV>(deformable_object.particles,curve_offset+2,*solids_parameters.rigid_body_parameters.list(1),VECTOR<T,2>()));
    solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new RIGID_BODY_BINDING<TV>(deformable_object.particles,curve_offset+3,*solids_parameters.rigid_body_parameters.list(2),VECTOR<T,2>()));
    solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new LINEAR_BINDING<TV,3>(deformable_object.particles,curve_offset+4,triangulated_area.mesh.elements(1),VECTOR<T,3>(1,0,0)));
    for(int i=0;i<segmented_curve.mesh.number_nodes;i++){
        int p=curve_offset+i;
        deformable_object.particles.X(p)=solid_body_collection.deformable_body_collection.binding_list.Average_Target_Position(p);
        deformable_object.particles.V(p)=solid_body_collection.deformable_body_collection.binding_list.Average_Target_Velocity(p);
        deformable_object.particles.mass(p)=solid_body_collection.deformable_body_collection.binding_list.Effective_Mass(p,VECTOR<T,2>(1,0));} // TODO: determine direction dynamically
    deformable_object.particles.mass.effective_mass=deformable_object.particles.mass.array;

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    Get_Initial_Data();

    DEFORMABLE_OBJECT<TV>& deformable_object=solid_body_collection.deformable_object;
    TRIANGULATED_AREA<T>& triangulated_area=deformable_object.template Find_Structure<TRIANGULATED_AREA<T>&>();
    PARTICLES<TV>& particles=deformable_object.particles;
    PARTICLE_SUBSET<PARTICLES<TV> >& area_particles=*new PARTICLE_SUBSET<PARTICLES<TV> >(particles);
    triangulated_area.mesh.elements.Flattened().Get_Unique(area_particles.active_indices);
    area_particles.Update_Subset_Index_From_Element_Index();
    solid_body_collection.Add_Force(new GRAVITY<TV>(area_particles));
    solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(triangulated_area,new DIAGONALIZED_SPLINE_MODEL_2D<T>((T)5e4,(T).3,(T).5,7,(T).01)));
    SEGMENTED_CURVE_2D<T>& segmented_curve=deformable_object.template Find_Structure<SEGMENTED_CURVE_2D<T>&>();
    solid_body_collection.Add_Force(Create_Edge_Springs(segmented_curve,(T)2e2,(T)1));

    solid_body_collection.Update_Fragments();

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Initialize_Bodies();
}
//#####################################################################
// Apply_Constraints
//#####################################################################
void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE
{
    DEFORMABLE_OBJECT<TV>& deformable_object=solid_body_collection.deformable_object;
    ARRAY<TV> F(deformable_object.particles.array_collection->Size());
    for(int f=0;f<deformable_object.particles_of_fragment.m;f++) solid_body_collection.deformable_object.Add_All_Forces(F,time,f);
    for(int b=0;b<solid_body_collection.deformable_body_collection.binding_list.bindings.m;b++){
        RIGID_BODY_BINDING<TV>* rigid_body_binding=dynamic_cast<RIGID_BODY_BINDING<TV>*>(solid_body_collection.deformable_body_collection.binding_list.bindings(b));
        if(rigid_body_binding) rigid_body_binding->Apply_Impulse_To_Parents_Based_On_Embedding(F(rigid_body_binding->particle_index)*dt);}
}
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_FLUID_COUPLING_TEST 
//##################################################################### 
#ifndef __SOLID_FLUID_COUPLING_TEST__
#define __SOLID_FLUID_COUPLING_TEST__

#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Math_Tools/choice.h>
#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include <Level_Sets/EXTRAPOLATION_2D.h>
#include <Rigid_Bodies/RIGID_BODY_FLUID_FORCES_2D.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_2D.h>
namespace PhysBAM{

template<class T,class RW=T>
class SOLID_FLUID_COUPLING_TEST:public SOLIDS_FLUIDS_EXAMPLE_2D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::solids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::verbose_dt;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::write_time;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::write_frame_title;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::data_directory;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::abort_when_dt_below;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Use_Thin_Shells_Fluid_Coupling_Defaults;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::set_fluid_boundary_conditions_for_pressure_jump_solve;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::need_to_compute_pressure_coupling;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::collision_body_affected_by_fluid;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::use_object_pseudo_velocity_for_boundary_conditions;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::rigid_body_fluid_forces;

    int example_number;
    bool use_source;
    BOX_2D<T> source;
    MATRIX<T,3> world_to_source;
    VECTOR_2D<T> source_velocity;
    T initial_water_level;

    PARAMETER_LIST parameter_list;
    bool use_source_override;
    T domain_scale;

    ARRAY<int,VECTOR<int,2> > id_objects;
    ARRAY<T,VECTOR<int,2> > phi_objects;
    ARRAY<VECTOR_2D<T> ,VECTOR<int,2> > V_objects;
    bool add_gravity_to_object;
    bool use_next_dt_for_velocity_adjustment;
    bool subtract_last_pressure_force_from_object;
    bool extrapolate_combined_solid_fluid_velocity;
    bool set_faces_to_invalid;
    T velocity_extrapolation_band_width;
    bool use_extrapolated_velocity_for_laplace;
    bool only_use_object_velocity_in_mixed_face;
    bool average_velocities_in_mixed_face;
    bool mixed_face_is_solid;
    bool use_effective_velocity_for_poisson;

    ARRAY<RIGID_BODY<TV>*> bodies_in_fluid;
    ARRAY<T> body_densities;

    SOLID_FLUID_COUPLING_TEST(const int example_number_input)
        :SOLIDS_FLUIDS_EXAMPLE_2D<RW>(fluids_parameters.WATER),example_number(example_number_input),domain_scale(1)
    {
        std::cout << "Example number " << example_number << std::endl;
        {std::string filename=STRING_UTILITIES::string_sprintf("Solid_Fluid_Coupling_Test/example_%d.param",example_number);
        if(FILE_UTILITIES::File_Exists(filename)){std::cout << "Reading parameter file '" << filename << "'" << std::endl;parameter_list.Read(filename);}}

        // Two way coupling
        fluids_parameters.use_solid_fluid_coupling=false;
        Set_Parameter(fluids_parameters.use_solid_fluid_coupling,"fluids_parameters.use_solid_fluid_coupling");
        if(fluids_parameters.use_solid_fluid_coupling) Use_Thin_Shells_Fluid_Coupling_Defaults();

        // TODO: is lockstep necessary?
        fluids_parameters.use_volumetric_solid_fluid_coupling=true;
        Set_Parameter(fluids_parameters.use_volumetric_solid_fluid_coupling,"fluids_parameters.use_volumetric_solid_fluid_coupling");
        add_gravity_to_object=false;
        Set_Parameter(add_gravity_to_object,"add_gravity_to_object");
        use_next_dt_for_velocity_adjustment=false;
        Set_Parameter(use_next_dt_for_velocity_adjustment,"use_next_dt_for_velocity_adjustment");
        subtract_last_pressure_force_from_object=false;
        Set_Parameter(subtract_last_pressure_force_from_object,"subtract_last_pressure_force_from_object");
        extrapolate_combined_solid_fluid_velocity=false;
        Set_Parameter(extrapolate_combined_solid_fluid_velocity,"extrapolate_combined_solid_fluid_velocity");
        fluids_parameters.scalar_substeps=(T).95;
        set_faces_to_invalid=false;
        Set_Parameter(set_faces_to_invalid,"set_faces_to_invalid");
        velocity_extrapolation_band_width=0;
        use_extrapolated_velocity_for_laplace=true;
        Set_Parameter(use_extrapolated_velocity_for_laplace,"use_extrapolated_velocity_for_laplace");
        only_use_object_velocity_in_mixed_face=false;
        Set_Parameter(only_use_object_velocity_in_mixed_face,"only_use_object_velocity_in_mixed_face");
        average_velocities_in_mixed_face=false;
        Set_Parameter(average_velocities_in_mixed_face,"average_velocities_in_mixed_face");
        mixed_face_is_solid=false;
        Set_Parameter(mixed_face_is_solid,"mixed_face_is_solid");
        use_effective_velocity_for_poisson=false;
        Set_Parameter(use_effective_velocity_for_poisson,"use_effective_velocity_for_poisson");

        Set_Parameter(domain_scale,"domain_scale");

        // Common parameters
        output_directory="Solid_Fluid_Coupling_Test/output";
        first_frame=0;last_frame=1000;frame_rate=60;
        restart=false;restart_frame=0;

        // Fluids parameters
        fluids_parameters.grid.Initialize(101,101,0,domain_scale,0,domain_scale);
        fluids_parameters.number_particles_per_cell=0;
        fluids_parameters.reseeding_frame_rate=10;
        fluids_parameters.incompressible_iterations=1000;
        fluids_parameters.cfl=0.9;
        fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[1][2]=true;fluids_parameters.domain_walls[2][1]=true;fluids_parameters.domain_walls[2][2]=false;
        fluids_parameters.bias_towards_negative_particles=true;
        fluids_parameters.use_removed_positive_particles=false;fluids_parameters.use_removed_negative_particles=true;
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;fluids_parameters.variable_viscosity=false;fluids_parameters.second_order_cut_cell_method=false;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;
        fluids_parameters.write_particles=true;fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=true;
        fluids_parameters.store_particle_ids=true;
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.enforce_divergence_free_extrapolation=false;
        fluids_parameters.use_old_velocities_for_boundary_conditions=false;

        // Solids parameters
        solids_parameters.gravity=9.8;
        solids_parameters.cfl=(T).9;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies=true;
        solids_parameters.deformable_body_parameters.write=true;
        verbose_dt=true;

        initial_water_level=domain_scale*0.6123;
        Set_Parameter(initial_water_level,"initial_water_level");
        fluids_parameters.simulate=true;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        use_source=false;

        // Debugging
        fluids_parameters.write_debug_data=true;
        abort_when_dt_below=1e-7;
        write_time=true;
        write_frame_title=true;

        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Parameters_From_Parameter_List(*this,parameter_list);
        use_source_override=true;Set_Parameter(use_source_override,"use_source");
        Set_Parameter(source_velocity,"source_velocity");
    }

    ~SOLID_FLUID_COUPLING_TEST() 
    {}

    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return false;}

//#####################################################################
// Function Set_Parameter
//#####################################################################
template<class T2> void Set_Parameter(T2& parameter,const std::string& name)
{
    if(parameter_list.Is_Defined(name)){
        parameter=parameter_list.template Get_Parameter<T2>(name);
        std::cout << "[param] set " << name << " = " << parameter << std::endl;}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    int id;
    T mu=0.8,epsilon=0.2;
    Set_Parameter(mu,"mu");
    Set_Parameter(epsilon,"epsilon");

    T density_factor,pressure_force_scale;
    RIGID_BODY<TV>* rigid_body=0;

    bool two_way_coupling=true;
    Set_Parameter(two_way_coupling,"two_way_coupling");

    int number_of_bodies=3;
    Set_Parameter(number_of_bodies,"number_of_bodies");

    if(example_number==8){
        std::string prefix=STRING_UTILITIES::string_sprintf("rigid_body_%d.",1);
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(true,data_directory+"/Rigid_Bodies_2D/cross",.2);
        density_factor=1;
        pressure_force_scale=1;
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Rigid_Body_Parameters_2D(id,density_factor,pressure_force_scale,parameter_list);
        RIGID_BODY<TV>* rigid_body=solids_parameters.rigid_body_parameters.list(id);
        rigid_body->Set_Mass(rigid_body->mass*density_factor*fluids_parameters.density/rigid_body->Volumetric_Density());
        rigid_body->Set_Coefficient_Of_Friction(mu);
        rigid_body->Set_Coefficient_Of_Restitution(epsilon);
        rigid_body->position=VECTOR_2D<T>(domain_scale*.50333,domain_scale*.50333);
        Set_Parameter(rigid_body->position,prefix+"position");
        Set_Parameter(rigid_body->orientation,prefix+"orientation");
        Set_Parameter(rigid_body->velocity,prefix+"velocity");
        Set_Parameter(rigid_body->angular_velocity,prefix+"angular_velocity");
        rigid_body->Update_Angular_Momentum();
        bool is_static=false;
        bool basic_forces=true;
        Set_Parameter(basic_forces,prefix+"basic_forces");
        if(basic_forces) rigid_body->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction.Vector_2D(),solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
        if(!is_static){Add_Volumetric_Body_To_Fluid_Simulation(*rigid_body,true,two_way_coupling,pressure_force_scale);
            bodies_in_fluid.Append(rigid_body);body_densities.Append(rigid_body->Volumetric_Density());}
    }
    else if(example_number==10){
        for(int i=0;i<2;i++){
            std::string prefix=STRING_UTILITIES::string_sprintf("rigid_body_%d.",i);
            id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(true,data_directory+"/Rigid_Bodies_2D/square_extra_refined",.1);
            density_factor=1;pressure_force_scale=1;
            THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Rigid_Body_Parameters_2D(id,density_factor,pressure_force_scale,parameter_list);
            RIGID_BODY<TV>* rigid_body=solids_parameters.rigid_body_parameters.list(id);
            rigid_body->Set_Mass(rigid_body->mass*density_factor*fluids_parameters.density/rigid_body->Volumetric_Density());
            rigid_body->Update_Angular_Momentum();
            rigid_body->rigid_body_collection.rigid_body_particle.kinematic(rigid_body->particle_index)=true;
            Add_Volumetric_Body_To_Fluid_Simulation(*rigid_body,true,false);
            bodies_in_fluid.Append(rigid_body);body_densities.Append(rigid_body->Volumetric_Density());}
    }
    else for(int i=0;i<number_of_bodies;i++){
        std::string prefix=STRING_UTILITIES::string_sprintf("rigid_body_%d.",i);
        T scale=(T).08123*domain_scale;Set_Parameter(scale,prefix+"scale");
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>(true,data_directory+"/Rigid_Bodies_2D/square_extra_refined",scale);
        if(i==1) density_factor=0.5; else if(i==2) density_factor=1; else if(i==3) density_factor=2; else density_factor=0.5;
        pressure_force_scale=1;
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Rigid_Body_Parameters_2D(id,density_factor,pressure_force_scale,parameter_list);
        RIGID_BODY<TV>* rigid_body=solids_parameters.rigid_body_parameters.list(id);
        rigid_body->Set_Mass(rigid_body->mass*density_factor*fluids_parameters.density/rigid_body->Volumetric_Density());
        rigid_body->Set_Coefficient_Of_Friction(mu);
        rigid_body->Set_Coefficient_Of_Restitution(epsilon);
        if(i==1) rigid_body->position=VECTOR_2D<T>(domain_scale*.20333,domain_scale*.750333);
        else if(i==2) rigid_body->position=VECTOR_2D<T>(domain_scale*.50333,domain_scale*.750333);
        else if(i==3) rigid_body->position=VECTOR_2D<T>(domain_scale*.80333,domain_scale*.750333);
        Set_Parameter(rigid_body->position,prefix+"position");
        Set_Parameter(rigid_body->orientation,prefix+"orientation");
        Set_Parameter(rigid_body->velocity,prefix+"velocity");
        Set_Parameter(rigid_body->angular_velocity,prefix+"angular_velocity");
        rigid_body->Update_Angular_Momentum();
        bool is_static=false;
        bool basic_forces=true;
        Set_Parameter(is_static,prefix+"is_static");
        if(is_static){rigid_body->is_static=true;basic_forces=false;}
        Set_Parameter(basic_forces,prefix+"basic_forces");
        if(basic_forces) rigid_body->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction.Vector_2D(),solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
        if(!is_static){Add_Volumetric_Body_To_Fluid_Simulation(*rigid_body,true,two_way_coupling,pressure_force_scale);
            bodies_in_fluid.Append(rigid_body);body_densities.Append(rigid_body->Volumetric_Density());}}

    THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Add_Rigid_Body_Walls(*this,epsilon,mu);

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    solids_parameters.rigid_body_parameters.Reset_Kinematic_Rigid_Bodies(0);
    SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Initialize_Bodies();
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    INTERPOLATION_CURVE<T,VECTOR_2D<T> > position;
    if(example_number==10){
        if(id==int(1)){
            position.Add_Control_Point(0,VECTOR_2D<T>((T).2,(T).5));
            position.Add_Control_Point(1,VECTOR_2D<T>((T).3993,(T).5));
            position.Add_Control_Point(9,VECTOR_2D<T>((T).3993,(T).5));
            position.Add_Control_Point(10,VECTOR_2D<T>((T).2,(T).5));}
        else if(id==int(2)){
            position.Add_Control_Point(0,VECTOR_2D<T>((T).8,(T).5));
            position.Add_Control_Point(1,VECTOR_2D<T>((T).6007,(T).5));
            position.Add_Control_Point(9,VECTOR_2D<T>((T).6007,(T).5));
            position.Add_Control_Point(10,VECTOR_2D<T>((T).8,(T).5));}
        frame.t=position.Value(time);}
}
//#####################################################################
// Function Construct_Levelsets_For_Objects
//#####################################################################
void Construct_Levelsets_For_Objects(const T time)
{
    GRID<TV>& grid=fluids_parameters.grid;
    phi_objects.Resize(grid,3);ARRAY<T,VECTOR<int,2> >::copy(5*grid.max_dx_dy,phi_objects);
    id_objects.Resize(grid,3);ARRAY<int,VECTOR<int,2> >::copy(0,id_objects);
    for(int object=0;object<bodies_in_fluid.m;object++){
        RIGID_BODY<TV>& rigid_body=*bodies_in_fluid(object);
        for(int i=-2;i<=grid.m+3;i++)for(int j=-2;j<=grid.n+3;j++){
            T phi_object=rigid_body.Implicit_Curve_Extended_Value(grid.X(i,j));
            if(phi_object<phi_objects(i,j)){
                id_objects(i,j)=object;phi_objects(i,j)=phi_object;}}}

    // Make signed distance
    FAST_LEVELSET_2D<T> levelset_objects(grid,phi_objects);levelset_objects.Set_Band_Width(3);levelset_objects.Reinitialize();

    // use effective velocity whenever possible
    bool use_effective_velocity=false; if(bodies_in_fluid.m&&bodies_in_fluid(1)->saved_states.m>=2) use_effective_velocity=true;
    if(!use_effective_velocity) std::cout << "Construct_Levelsets_For_Objects: not using effectively velocity (possibly start of simulation)" << std::endl;
    Construct_Velocities_For_Objects(use_effective_velocity,false,false,0,time);
}
//#####################################################################
// Function Construct_Velocities_For_Objects
//#####################################################################
void Construct_Velocities_For_Objects(const bool use_effective_velocity,const bool subtract_last_pressure_force,const bool add_gravity,const T dt,const T time)
{
    T dt_for_velocity_adjustment=use_next_dt_for_velocity_adjustment?fluids_parameters.next_dt:dt;

    VECTOR_2D<T> extra_velocity;
    if(add_gravity){
        if(use_next_dt_for_velocity_adjustment) std::cout << "Using next dt " << fluids_parameters.next_dt << " for gravity" << std::endl;
        extra_velocity=dt_for_velocity_adjustment*solids_parameters.gravity*solids_parameters.gravity_direction.Vector_2D();
    }

    ARRAY<RIGID_BODY_STATE_2D<T> > state_for_object_velocity(bodies_in_fluid.m);
    if(subtract_last_pressure_force){
        if(use_next_dt_for_velocity_adjustment) std::cout << "Using next dt " << fluids_parameters.next_dt << " to subtract last pressure force" << std::endl;
        for(int i=0;i<bodies_in_fluid.m;i++){
            assert(&rigid_body_fluid_forces(i)->rigid_body==bodies_in_fluid(i));
            state_for_object_velocity(i)=*bodies_in_fluid(i)->saved_states(COLLISION_BODY_2D<T>::THIN_SHELLS_NEW_STATE);
            std::cout << i << ":  VELOCITY BEFORE " << state_for_object_velocity(i).velocity << std::endl;
            state_for_object_velocity(i).velocity-=dt_for_velocity_adjustment*rigid_body_fluid_forces(i)->last_computed_force/bodies_in_fluid(i)->mass;
            state_for_object_velocity(i).angular_momentum-=dt_for_velocity_adjustment*rigid_body_fluid_forces(i)->last_computed_torque;
            bodies_in_fluid(i)->Update_Angular_Velocity(state_for_object_velocity(i));
            std::cout << i << ":  VELOCITY AFTER " << state_for_object_velocity(i).velocity << std::endl;}}

    GRID<TV>& grid=fluids_parameters.grid;
    V_objects.Resize(grid,3);ARRAY<VECTOR_2D<T> ,VECTOR<int,2> >::copy(VECTOR_2D<T>(0,0),V_objects);
    for(int i=-2;i<=grid.m+3;i++)for(int j=-2;j<=grid.n+3;j++)if(phi_objects(i,j)<0){
        int body_id=id_objects(i,j);
        if(use_effective_velocity) V_objects(i,j)=bodies_in_fluid(body_id)->Pointwise_Object_Pseudo_Velocity(0,grid.X(i,j),
                                                    COLLISION_BODY_2D<T>::THIN_SHELLS_OLD_STATE,COLLISION_BODY_2D<T>::THIN_SHELLS_NEW_STATE);
        else if(subtract_last_pressure_force) V_objects(i,j)=bodies_in_fluid(body_id)->Pointwise_Object_Velocity(grid.X(i,j),state_for_object_velocity(body_id));
        else V_objects(i,j)=extra_velocity+bodies_in_fluid(body_id)->Pointwise_Object_Velocity(grid.X(i,j));}
}
//#####################################################################
// Function Extrapolate_Phi_Into_Objects
//#####################################################################
void Extrapolate_Phi_Into_Objects(const T time)
{
    fluids_parameters.Extrapolate_Phi_Into_Object(phi_objects);
}
//#####################################################################
// Function Adjust_Phi_with_Objects
//#####################################################################
void Adjust_Phi_With_Objects(const T time)
{
    fluids_parameters.Adjust_Phi_With_Object(phi_objects,V_objects,time);
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    // Not so good to set up a heaviside function here because then the interface will
    // be exactly between the two nodes which can lead to roundoff issues when setting dirichlet cells, etc.
    GRID<TV>& grid=fluids_parameters.grid;
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++){
        fluids_parameters.particle_levelset_evolution.phi(i,j)=grid.y(j)-grid.ymin-initial_water_level;
        if(phi_objects(i,j)<0) fluids_parameters.particle_levelset_evolution.phi(i,j)=-phi_objects(i,j);
    }
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
void Adjust_Phi_With_Sources(const T time)
{
    if(use_source) SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Adjust_Phi_With_Source(source,world_to_source,time);
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,2> >*& cell_centered_mask,const T time)
{
    if(use_source) SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Get_Source_Reseed_Mask(source,world_to_source,cell_centered_mask,true,time);
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time)
{
    if(use_source) SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Get_Source_Velocities(source,world_to_source,source_velocity,time);
}
//#####################################################################
// Utilities
//#####################################################################
bool Mixed_Face(const int axis,const VECTOR_2D<T>& face_index){return phi_objects(GRID_2D<T>::Face_Node_Index(axis,face_index,1))<0 || phi_objects(GRID_2D<T>::Face_Node_Index(axis,face_index,2))<0;}
//#####################################################################
// Function Average_Node_Velocities_To_Faces
//#####################################################################
void Average_Node_Velocities_To_Faces()
{
    if(set_faces_to_invalid){
        std::cout << "Setting psi_D_old to all true" << std::endl;
        ARRAY<bool,VECTOR<int,2> >::copy(true,fluids_parameters.incompressible.psi_D_old);
    }

    if(!set_fluid_boundary_conditions_for_pressure_jump_solve){
        // Average from nodes to faces but ignore nodes inside objects
        // Note: if both nodes are in object the face must be psi_N so we don't process it

        PROJECTION_2D<T>& projection=fluids_parameters.incompressible.projection;
        ARRAY<VECTOR_2D<T> ,VECTOR<int,2> > &V=fluids_parameters.incompressible.V,&V_ghost=fluids_parameters.incompressible.V_ghost;

        for(int axis=0;axis<2;axis++){
            GRID<TV>& face_grid=choice(axis,projection.u_grid,projection.v_grid);
            ARRAY<bool,VECTOR<int,2> >& psi_N=choice(axis,projection.elliptic_solver->psi_N_u,projection.elliptic_solver->psi_N_v);
            ARRAY<T,VECTOR<int,2> >& face_velocity=choice(axis,projection.u,projection.v);

            for(int i=0;i<face_grid.m;i++) for(int j=0;j<face_grid.n;j++) if(!psi_N(i,j)){
                VECTOR_2D<int> face_index(i,j),first_cell_index=face_index,second_cell_index=face_index;first_cell_index[axis]--;
                int count=0;T velocity=0;bool mixed_face=Mixed_Face(axis,face_index);
                if(only_use_object_velocity_in_mixed_face && mixed_face){
                    for(int c=0;c<2;c++){VECTOR_2D<int> index=GRID_2D<T>::Face_Node_Index(axis,face_index,c);
                        if(phi_objects(index)<0){velocity+=V(index)[axis];count++;}}
                    assert(count);face_velocity(face_index)=velocity/(T)count;}
                else if(average_velocities_in_mixed_face && mixed_face){
                    for(int c=0;c<2;c++){velocity+=V(GRID_2D<T>::Face_Node_Index(axis,face_index,c))[axis];}
                    face_velocity(face_index)=(T).5*velocity;}
                else if(!fluids_parameters.incompressible.psi_D_old(first_cell_index)||!fluids_parameters.incompressible.psi_D_old(second_cell_index)){
                    for(int c=0;c<2;c++){VECTOR_2D<int> index=GRID_2D<T>::Face_Node_Index(axis,face_index,c);
                        if(phi_objects(index)>=0){velocity+=V(index)[axis]-V_ghost(index)[axis];count++;}}
                    assert(count);face_velocity(face_index)+=velocity/(T)count;}
                else{
                    for(int c=0;c<2;c++){VECTOR_2D<int> index=GRID_2D<T>::Face_Node_Index(axis,face_index,c);
                        if(phi_objects(index)>=0){velocity+=V(index)[axis];count++;}}
                    assert(count);face_velocity(face_index)=velocity/(T)count;}}
            }
    }
    else SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Average_Node_Velocities_To_Faces();
}
//#####################################################################
// Function Adjust_Velocity_With_Objects
//#####################################################################
void Adjust_Velocity_With_Objects(const T time)
{
    // V_objects should be set with the effective velocity at this point (and phi_objects is in object's initial position for this time step)
    fluids_parameters.Extrapolate_Velocity_Into_Object(phi_objects,V_objects,velocity_extrapolation_band_width,false,time);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time)
{
    if(!set_fluid_boundary_conditions_for_pressure_jump_solve){
        GRID<TV> &grid=fluids_parameters.grid,&u_grid=fluids_parameters.u_grid,&v_grid=fluids_parameters.v_grid;
        LAPLACE_2D<T>& elliptic_solver=*fluids_parameters.incompressible.projection.elliptic_solver;

        Construct_Velocities_For_Objects(use_object_pseudo_velocity_for_boundary_conditions,false,false,dt,time);

        ARRAY<VECTOR_2D<T> ,VECTOR<int,2> >& V=fluids_parameters.incompressible.V;
        if(use_extrapolated_velocity_for_laplace) fluids_parameters.Extrapolate_Velocity_Into_Object(phi_objects,V_objects,velocity_extrapolation_band_width,false,time);
        else for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) if(phi_objects(i,j)<0) V(i,j)=V_objects(i,j);

        // Get psi_N velocities by selectively averaging only from object when at least one node is inside object
        for(int i=0;i<u_grid.m;i++) for(int j=0;j<u_grid.n;j++) 
            if((!mixed_face_is_solid && phi_objects(i,j)+phi_objects(i,j+1)<0) || (mixed_face_is_solid && (phi_objects(i,j)<0 || phi_objects(i,j+1)<0))){
                int count=0;T velocity=0;
                if(phi_objects(i,j)<0 /*&& collision_body_affected_by_fluid(id_objects(i,j))*/){velocity+=V(i,j).x;count++;}
                if(phi_objects(i,j+1)<0 /*&& collision_body_affected_by_fluid(id_objects(i,j+1))*/){velocity+=V(i,j+1).x;count++;}
                assert(count);if(count==2){velocity*=(T).5;}
                elliptic_solver.psi_N_u(i,j)=true;fluids_parameters.incompressible.projection.u(i,j)=velocity;}

        for(int i=0;i<v_grid.m;i++) for(int j=0;j<v_grid.n;j++) 
            if((!mixed_face_is_solid && phi_objects(i,j)+phi_objects(i+1,j)<0) || (mixed_face_is_solid && (phi_objects(i,j)<0 || phi_objects(i+1,j)<0))){
                int count=0;T velocity=0;
                if(phi_objects(i,j)<0 /*&& collision_body_affected_by_fluid(id_objects(i,j))*/){velocity+=V(i,j).y;count++;}
                if(phi_objects(i+1,j)<0 /*&& collision_body_affected_by_fluid(id_objects(i+1,j))*/){velocity+=V(i+1,j).y;count++;}
                assert(count);if(count==2){velocity*=(T).5;}
                elliptic_solver.psi_N_v(i,j)=true;fluids_parameters.incompressible.projection.v(i,j)=velocity;}
    }
    else{
        GRID<TV> &grid=fluids_parameters.grid,&u_grid=fluids_parameters.u_grid,&v_grid=fluids_parameters.v_grid;
        POISSON_2D<T>* poisson=fluids_parameters.incompressible.projection.poisson;
        assert(poisson && need_to_compute_pressure_coupling && (fluids_parameters.fire || poisson->use_variable_beta));
        ARRAY<T,VECTOR<int,2> >::copy((T)1/fluids_parameters.density,poisson->beta_right);ARRAY<T,VECTOR<int,2> >::copy((T)1/fluids_parameters.density,poisson->beta_top);
        poisson->Use_Body_Psi_N(true);
        ARRAY<bool,VECTOR<int,2> >::copy(false,poisson->body_psi_N_u);ARRAY<bool,VECTOR<int,2> >::copy(false,poisson->body_psi_N_v);

        Construct_Velocities_For_Objects(use_effective_velocity_for_poisson,subtract_last_pressure_force_from_object,add_gravity_to_object,dt,time);

        if(extrapolate_combined_solid_fluid_velocity){ // extrapolate combined fluid and solid velocity into air
            ARRAY<T,VECTOR<int,2> > combined_solid_fluid_phi(grid);
            FAST_LEVELSET_2D<T> combined_solid_fluid_levelset(grid,combined_solid_fluid_phi);
            combined_solid_fluid_levelset.Set_Band_Width(3);combined_solid_fluid_levelset.Reinitialize();
            for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++){
                combined_solid_fluid_phi(i,j)=min(phi_objects(i,j),fluids_parameters.particle_levelset_evolution.phi(i,j));
                if(phi_objects(i,j)<0) fluids_parameters.incompressible.V(i,j)=V_objects(i,j);}
            EXTRAPOLATION_2D<T,VECTOR_2D<T> > extrapolate(grid,combined_solid_fluid_phi,fluids_parameters.incompressible.V);extrapolate.Set_Band_Width(3);extrapolate.Extrapolate();}

        for(int i=0;i<u_grid.m;i++) for(int j=0;j<u_grid.n;j++)
            if((!mixed_face_is_solid && phi_objects(i,j)+phi_objects(i,j+1)<0) || (mixed_face_is_solid && (phi_objects(i,j)<0 || phi_objects(i,j+1)<0))){
                int count=0;T density=0,velocity=0;
                if(phi_objects(i,j)<0 && collision_body_affected_by_fluid(id_objects(i,j))){int body_id=id_objects(i,j);
                    density+=body_densities(body_id);velocity+=V_objects(i,j).x;count++;}
                if(phi_objects(i,j+1)<0 && collision_body_affected_by_fluid(id_objects(i,j+1))){int body_id=id_objects(i,j+1);
                    density+=body_densities(body_id);velocity+=V_objects(i,j+1).x;count++;}
                if(count==2){density*=(T).5;velocity*=(T).5;}
                poisson->psi_N_u(i,j)=true;fluids_parameters.incompressible.projection.u(i,j)=velocity;
                poisson->body_psi_N_u(i,j)=true;poisson->beta_right(i-1,j)=1/density;}

        for(int i=0;i<v_grid.m;i++) for(int j=0;j<v_grid.n;j++)
            if((!mixed_face_is_solid && phi_objects(i,j)+phi_objects(i+1,j)<0) || (mixed_face_is_solid && (phi_objects(i,j)<0 || phi_objects(i+1,j)<0))){
                int count=0;T density=0,velocity=0;
                if(phi_objects(i,j)<0 && collision_body_affected_by_fluid(id_objects(i,j))){int body_id=id_objects(i,j);
                    density+=body_densities(body_id);velocity+=V_objects(i,j).y;count++;}
                if(phi_objects(i+1,j)<0 && collision_body_affected_by_fluid(id_objects(i+1,j))){int body_id=id_objects(i+1,j);
                    density+=body_densities(body_id);velocity+=V_objects(i+1,j).y;count++;}
                if(count==2){density*=(T).5;velocity*=(T).5;}
                poisson->psi_N_v(i,j)=true;fluids_parameters.incompressible.projection.v(i,j)=velocity;
                poisson->body_psi_N_v(i,j)=true;poisson->beta_top(i,j-1)=1/density;}}
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
void Set_Dirichlet_Boundary_Conditions(const T time)
{
    if(!set_fluid_boundary_conditions_for_pressure_jump_solve){
        SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Set_Dirichlet_Boundary_Conditions(time);
        return;}

    GRID<TV>& p_grid=fluids_parameters.incompressible.projection.p_grid;
    ARRAY<T,VECTOR<int,2> >& phi=fluids_parameters.particle_levelset_evolution.phi;
    for(int i=0;i<p_grid.m;i++) for(int j=0;j<p_grid.n;j++) 
        if((phi(i,j)+phi(i+1,j)+phi(i,j+1)+phi(i+1,j+1) > 0) && (phi_objects(i,j)+phi_objects(i+1,j)+phi_objects(i,j+1)+phi_objects(i+1,j+1)>0)){
            fluids_parameters.incompressible.projection.elliptic_solver->psi_D(i,j)=true;fluids_parameters.incompressible.projection.p(i,j)=0;}
}
//#####################################################################
// Function Adjust_Particle_For_Objects
//#####################################################################
bool Adjust_Particle_For_Objects(PARTICLE_LEVELSET_PARTICLES<T,VECTOR_2D<T> >& particles,const int index,VECTOR_2D<T>& V,const typename PARTICLE_LEVELSET<T,VECTOR_2D<T> >::PARTICLE_TYPE particle_type,const T dt,const T time)
{
    LEVELSET_2D<T> levelset(fluids_parameters.grid,phi_objects);
    fluids_parameters.Adjust_Particle_For_Object(levelset,V_objects,particles,index,V,particle_type,dt,time);
    return true;
} 
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
void Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<T,VECTOR_2D<T> >& particles,const typename PARTICLE_LEVELSET<T,VECTOR_2D<T> >::PARTICLE_TYPE particle_type,const T time)
{
    LEVELSET_2D<T> levelset_objects(fluids_parameters.grid,phi_objects);
    if(particle_type == PARTICLE_LEVELSET<T,VECTOR_2D<T> >::NEGATIVE || particle_type == PARTICLE_LEVELSET<T,VECTOR_2D<T> >::REMOVED_NEGATIVE){
        for(int k=particles.array_collection->Size();k>=1;k--) if(levelset_objects.Extended_Phi(particles.X(k))<-fluids_parameters.grid.min_dx_dy) particles.Delete_Particle(k);}
   else for(int k=particles.array_collection->Size();k>=1;k--) if(levelset_objects.Extended_Phi(particles.X(k))<0) particles.Delete_Particle(k);
}
//#####################################################################
// Function Adjust_Particle_For_Objects
//#####################################################################
void Write_Output_Files(const int frame) const
{
    const LEVELSET_2D<T> levelset_objects((GRID_2D<T>&)fluids_parameters.grid,(ARRAY<T,VECTOR<int,2> >&)phi_objects);
    FILE_UTILITIES::Write_To_File<RW>(STRING_UTILITIES::string_sprintf("%s/phi_objects.%d",output_directory.c_str(),frame),levelset_objects);
    ARRAY<PAIR<VECTOR_2D<T>,T> > rigid_body_forces_and_torques(solids_parameters.rigid_body_parameters.list.Number_Of_Elements());
    for(int i=0;i<solids_parameters.rigid_body_parameters.list.Number_Of_Elements();i++){const RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list(i);
        if(rigid_body.fluid_forces.m){rigid_body_forces_and_torques(i).x=rigid_body.fluid_forces(1)->last_computed_force;rigid_body_forces_and_torques(i).y=rigid_body.fluid_forces(1)->last_computed_torque;}}
    FILE_UTILITIES::Write_To_File<RW>(STRING_UTILITIES::string_sprintf("%s/rigid_body_forces_and_torques.%d",output_directory.c_str(),frame),rigid_body_forces_and_torques);
    SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Write_Output_Files(frame);
}
//#####################################################################
};    
}
#endif

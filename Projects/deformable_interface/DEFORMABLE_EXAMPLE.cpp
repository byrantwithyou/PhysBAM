#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Collision_Detection/COLLISION_GEOMETRY_SPATIAL_PARTITION.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/COROTATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLE_GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "DATA_EXCHANGE_CONVERSION.h"
#include "DEFORMABLE_EXAMPLE.h"
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> DEFORMABLE_EXAMPLE<T>::
DEFORMABLE_EXAMPLE(const STREAM_TYPE stream_type)
    :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),fully_implicit(false),driver(*new SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> >(*this)),added_body(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> DEFORMABLE_EXAMPLE<T>::
~DEFORMABLE_EXAMPLE()
{
    object_wrappers.Delete_Pointers_And_Clean_Memory();
    force_wrappers.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
template<class T> void DEFORMABLE_EXAMPLE<T>::
Initialize_Bodies() PHYSBAM_OVERRIDE
{
    fluids_parameters.simulate=false;
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.cfl=1;
    solids_parameters.triangle_collision_parameters.perform_self_collision=true;
    solids_parameters.triangle_collision_parameters.output_interaction_pairs=true;
    solids_parameters.rigid_body_collision_parameters.use_push_out=true;
    solids_parameters.use_rigid_deformable_contact=true;
    solids_parameters.rigid_body_collision_parameters.collision_bounding_box_thickness=(T)1e-3;
    solids_parameters.triangle_collision_parameters.collisions_output_number_checked=false;
    solid_body_collection.deformable_body_collection.soft_bindings.use_gauss_seidel_for_impulse_based_collisions=true;
    solids_parameters.verbose_dt=true;
    solids_parameters.triangle_collision_parameters.total_collision_loops=1;
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=20;
    solids_parameters.enforce_repulsions_in_cg=false;
    solids_parameters.enforce_poststabilization_in_cg=true;
    solids_parameters.triangle_collision_parameters.collisions_nonrigid_collision_attempts=30;
    solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness*=(T)5;
    solid_body_collection.collision_body_list.Set_Collision_Body_Thickness((T).1);
    solids_parameters.triangle_collision_parameters.collisions_collision_thickness=(T)1e-3;
    solids_parameters.implicit_solve_parameters.project_nullspace_frequency=1;
}
//#####################################################################
// Function Add_New_Forces
//#####################################################################
template<class T> void DEFORMABLE_EXAMPLE<T>::
Add_New_Forces()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    for(int i=1;i<=soft_bindings.use_impulses_for_collisions.m;i++)
        if(soft_bindings.use_impulses_for_collisions(i))
            soft_bindings.bindings_using_impulses_for_collisions.Append(i);

    binding_list.Clamp_Particles_To_Embedded_Positions();
    binding_list.Clamp_Particles_To_Embedded_Velocities();
    soft_bindings.Clamp_Particles_To_Embedded_Positions(true);
    soft_bindings.Clamp_Particles_To_Embedded_Velocities(true);

    // correct number nodes
    for(int i=1;i<=deformable_body_collection.deformable_geometry.structures.m;i++)
        deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();

    for(int i=1;i<=new_forces_relations.m;i++){
        PHYSBAM_ASSERT(new_forces_relations(i).x->id==DEFORMABLE_BODY_WRAPPER::fixed_id);
        Instantiate_Force(new_forces_relations(i).y,*static_cast<DEFORMABLE_BODY_WRAPPER*>(new_forces_relations(i).x));}
    new_forces_relations.Remove_All();
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T> void DEFORMABLE_EXAMPLE<T>::
Write_Output_Files(const int frame) const
{
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    binding_list.Clamp_Particles_To_Embedded_Velocities();
    BASE::Write_Output_Files(frame);
}
//#####################################################################
// Function Add_Deformable_Body
//#####################################################################
template<class T> DEFORMABLE_BODY_WRAPPER* DEFORMABLE_EXAMPLE<T>::
Add_Deformable_Body(const data_exchange::deformable_body& body)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    BINDING_LIST<TV>& binding_list=deformable_body_collection.binding_list;
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    ARRAY<int> parent_index;
    DEFORMABLE_BODY_WRAPPER* wrap=new DEFORMABLE_BODY_WRAPPER(*this);
    Triangulated_Surface_From_Data_Exchange(*surface,body.mesh,body.position,&parent_index);
    T density=body.mass?body.mass/surface->Volumetric_Volume():1000;

    TETRAHEDRALIZED_VOLUME<T>* tet_volume=0;
    TRIANGULATED_SURFACE<T>* tri_surface=0;
    tests.Create_Regular_Embedded_Surface(binding_list,soft_bindings,*surface,density,1000,(T)1e-5,wrap->particle_map,&tri_surface,&tet_volume,true);

    wrap->structure_index=deformable_body_collection.deformable_geometry.structures.m;
    wrap->free_particles_index=wrap->structure_index-1;
    wrap->enclosing_structure_index=wrap->structure_index-2;

    tests.Substitute_Soft_Bindings_For_Embedded_Nodes(*tri_surface,soft_bindings);

    deformable_body_collection.collisions.collision_structures.Append(deformable_body_collection.deformable_geometry.structures(wrap->free_particles_index));
    if(solids_parameters.triangle_collision_parameters.perform_self_collision)
        solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(tri_surface);

    return wrap;
}
//#####################################################################
// Function Add_Rigid_Body
//#####################################################################
template<class T> SCRIPTED_GEOMETRY_WRAPPER* DEFORMABLE_EXAMPLE<T>::
Add_Rigid_Body(const data_exchange::ground_plane& body)
{
    SCRIPTED_GEOMETRY_WRAPPER* wrap=new SCRIPTED_GEOMETRY_WRAPPER(*this);
    RIGID_BODY<TV>& ground=tests.Add_Ground();
    ground.X()=TV(To_Pb(body.position));
    ground.Rotation()=ROTATION<TV>::From_Rotated_Vector(TV(0,1,0),TV(To_Pb(body.normal)));
    wrap->rigid_index=ground.particle_index;
    return wrap;
}
//#####################################################################
// Function Add_Rigid_Body
//#####################################################################
template<class T> SCRIPTED_GEOMETRY_WRAPPER* DEFORMABLE_EXAMPLE<T>::
Add_Rigid_Body(const data_exchange::scripted_geometry& body)
{
    SCRIPTED_GEOMETRY_WRAPPER* wrap=new SCRIPTED_GEOMETRY_WRAPPER(*this);
    PARTICLES<TV>& particles=*new PARTICLES<TV>;
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create(particles);
    Triangulated_Surface_From_Data_Exchange(*surface,body.mesh,body.position,0);
    RIGID_BODY<TV>& rigid_body=*tests.Create_Rigid_Body_From_Triangulated_Surface(*surface,solid_body_collection.rigid_body_collection,1,20);
    rigid_body.is_static=true;
    solid_body_collection.rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
    wrap->rigid_index=rigid_body.particle_index;
    return wrap;
}
//#####################################################################
// Function Instantiate_Force
//#####################################################################
template<class T> void DEFORMABLE_EXAMPLE<T>::
Instantiate_Force(GRAVITY_WRAPPER* wrapper,DEFORMABLE_BODY_WRAPPER& body)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    ARRAY<int>* influenced_particles=new ARRAY<int>(particles.X.m);

    TETRAHEDRALIZED_VOLUME<T>* volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(deformable_body_collection.deformable_geometry.structures(body.enclosing_structure_index));
    PHYSBAM_ASSERT(volume);
    volume->mesh.elements.Flattened().Get_Unique(*influenced_particles);
    DEFORMABLE_GRAVITY<TV>* gravity=new DEFORMABLE_GRAVITY<TV>(particles,influenced_particles,wrapper->magnitude,TV(wrapper->direction));
    gravity->Own_Influenced_Particles();
    deformable_body_collection.Add_Force(gravity);
    wrapper->force_instances.Append(gravity);

    gravity->limit_time_step_by_strain_rate=false;
    gravity->use_implicit_velocity_independent_forces=fully_implicit;
}
//#####################################################################
// Function Instantiate_Force
//#####################################################################
template<class T> void DEFORMABLE_EXAMPLE<T>::
Instantiate_Force(VOLUMETRIC_FORCE_WRAPPER* wrapper,DEFORMABLE_BODY_WRAPPER& body)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    TETRAHEDRALIZED_VOLUME<T>* volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(deformable_body_collection.deformable_geometry.structures(body.enclosing_structure_index));
    PHYSBAM_ASSERT(volume);
    FINITE_VOLUME<TV,3>* finite_volume=Create_Finite_Volume(*volume,new COROTATED<T,3>(wrapper->stiffness,wrapper->poissons_ratio,wrapper->damping),false);
    solid_body_collection.Add_Force(finite_volume);
    wrapper->force_instances.Append(finite_volume);

    finite_volume->limit_time_step_by_strain_rate=false;
    finite_volume->use_implicit_velocity_independent_forces=fully_implicit;
}
//#####################################################################
// Function Instantiate_Force
//#####################################################################
template<class T> void DEFORMABLE_EXAMPLE<T>::
Instantiate_Force(FORCE_WRAPPER* wrapper,DEFORMABLE_BODY_WRAPPER& body)
{
    switch(wrapper->id){
        case VOLUMETRIC_FORCE_WRAPPER::fixed_id: Instantiate_Force(static_cast<VOLUMETRIC_FORCE_WRAPPER*>(wrapper),body); break;
        case GRAVITY_WRAPPER::fixed_id: Instantiate_Force(static_cast<GRAVITY_WRAPPER*>(wrapper),body); break;
        default: PHYSBAM_FATAL_ERROR();}
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class T> FORCE_WRAPPER* DEFORMABLE_EXAMPLE<T>::
Add_Force(const data_exchange::force* f)
{
    if(const data_exchange::gravity_force* g=dynamic_cast<const data_exchange::gravity_force*>(f)){
        GRAVITY_WRAPPER* wrap = new GRAVITY_WRAPPER(*this);
        force_wrappers.Append(wrap);
        wrap->magnitude=g->magnitude;
        wrap->direction=To_Pb(g->direction);
        return wrap;}
    if(const data_exchange::volumetric_force* g=dynamic_cast<const data_exchange::volumetric_force*>(f)){
        VOLUMETRIC_FORCE_WRAPPER* wrap = new VOLUMETRIC_FORCE_WRAPPER(*this);
        force_wrappers.Append(wrap);
        wrap->stiffness=g->stiffness;
        wrap->poissons_ratio=g->poissons_ratio;
        wrap->damping=g->damping;
        return wrap;}

    return 0;
}
//#####################################################################
// Function Initialize_Simulation
//#####################################################################
template<class T> void DEFORMABLE_EXAMPLE<T>::
Initialize_Simulation()
{
    if(this->restart){driver.current_frame=this->restart_frame;driver.Read_Time(driver.current_frame);}else driver.current_frame=this->first_frame;
    driver.output_number=driver.current_frame;
    driver.time=this->Time_At_Frame(driver.current_frame);
    solid_body_collection.deformable_body_collection.mpi_solids=solid_body_collection.deformable_body_collection.mpi_solids;
    solids_evolution->Set_Solids_Evolution_Callbacks(*this);
}
//#####################################################################
// Function Initialize_After_New_Bodies
//#####################################################################
template<class T> void DEFORMABLE_EXAMPLE<T>::
Initialize_After_New_Bodies()
{
    solids_evolution->time=driver.time;
    solid_body_collection.deformable_body_collection.Initialize(solids_parameters.triangle_collision_parameters);
    if(this->restart){this->Read_Output_Files_Solids(this->restart_frame);solids_evolution->time=driver.time=this->Time_At_Frame(this->restart_frame);}
    solids_evolution->Initialize_Deformable_Objects(this->frame_rate,this->restart);
    solids_evolution->Initialize_Rigid_Bodies(this->frame_rate,this->restart);
    Post_Initialization();
    this->Log_Parameters();
    if(!this->restart) driver.Write_Output_Files(this->first_frame);
}
//#####################################################################
// Function Simulate_Frame
//#####################################################################
template<class T> void DEFORMABLE_EXAMPLE<T>::
Simulate_Frame()
{
    if(added_body){
        Initialize_After_New_Bodies();
        Add_New_Forces();
        added_body=false;}
    else if(new_forces_relations.m)
        Add_New_Forces();

    driver.current_frame++;
    LOG::SCOPE scope("FRAME","Frame %d",driver.current_frame);
    driver.Preprocess_Frame(driver.current_frame);
    this->solids_evolution->kinematic_evolution.Get_Current_Kinematic_Keyframes(this->Time_At_Frame(driver.current_frame)-driver.time,driver.time);
    driver.Advance_To_Target_Time(this->Time_At_Frame(driver.current_frame));
    driver.Postprocess_Frame(driver.current_frame);
    if(this->write_output_files && this->write_substeps_level==-1) Write_Output_Files(driver.current_frame);
    else if(this->write_substeps_level!=-1) driver.Write_Substep(STRING_UTILITIES::string_sprintf("END Frame %d",driver.current_frame),0,this->write_substeps_level);
    LOG::cout<<"TIME = "<<time<<std::endl;
}
//#####################################################################
// Function Add_Simulation_Object
//#####################################################################
template<class T> OBJECT_WRAPPER* DEFORMABLE_EXAMPLE<T>::
Add_Simulation_Object(const data_exchange::simulation_object& body)
{
    added_body=true;

    if(const data_exchange::deformable_body* p=dynamic_cast<const data_exchange::deformable_body*>(&body))
        return Add_Deformable_Body(*p);
    if(const data_exchange::ground_plane* p=dynamic_cast<const data_exchange::ground_plane*>(&body))
        return Add_Rigid_Body(*p);
    if(const data_exchange::scripted_geometry* p=dynamic_cast<const data_exchange::scripted_geometry*>(&body))
        return Add_Rigid_Body(*p);

    return 0;
}
//#####################################################################
// Constructor
//#####################################################################
OBJECT_WRAPPER::OBJECT_WRAPPER(DEFORMABLE_EXAMPLE<float>& de_input,int id_input)
    :id(id_input),de(de_input)
{
    de.object_wrappers.Append(this);
}
//#####################################################################
// Constructor
//#####################################################################
FORCE_WRAPPER::FORCE_WRAPPER(DEFORMABLE_EXAMPLE<float>& de_input,int id_input)
    :id(id_input),de(de_input)
{
    de.force_wrappers.Append(this);
}
template class DEFORMABLE_EXAMPLE<float>;

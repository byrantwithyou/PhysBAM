//#####################################################################
// Copyright 2007-2008, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BINDING_PLASTICITY_EXAMPLE
//##################################################################### 
#ifndef __BINDING_PLASTICITY_EXAMPLE__
#define __BINDING_PLASTICITY_EXAMPLE__

#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/PARTICLE_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BINDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDING.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Heat_Flows/HEAT_UNIFORM.h>
#include <PhysBAM_Dynamics/Meshing/RED_GREEN_TRIANGLES.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class BINDING_PLASTICITY_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::solids_parameters;
    using BASE::fluids_parameters;using BASE::write_last_frame;using BASE::data_directory;using BASE::stream_type;using BASE::solid_body_collection;using BASE::parse_args;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions; // silence -Woverloaded-virtual

    SOLIDS_STANDARD_TESTS<TV> tests;
    RED_GREEN_TRIANGLES<TV>* redgreen;
    int refinement_level;
    IMPLICIT_ZERO_LENGTH_SPRINGS<TV>* spring_force;
    T plastic_yield_threshold;
    T letters_initial_height,letters_scale,base_scale,base_thickness,stamp_depth;
    T final_time;
    ARRAY<FRAME<TV> > letters_frames_save;
    FRAME<TV> base_initial_frame;
    int base_id,ground_id,handle_id;
    FRAME<TV> letters_to_base_frame;
    int version;
    INTERPOLATION_CURVE<T,TV> translation_curve;
    INTERPOLATION_CURVE<T,T> angle_curve;

    BINDING_PLASTICITY_EXAMPLE(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),redgreen(0),refinement_level(3),spring_force(0),plastic_yield_threshold((T).07),
        letters_initial_height(10),letters_scale(3),base_scale((T)2.5),base_thickness((T).5*base_scale),stamp_depth((T).2),final_time(2),letters_frames_save(8),base_id(0),ground_id(0),
        handle_id(0),version(2)
    {
    }

    // unused callbacks
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options()
{
    BASE::Register_Options();
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options()
{
    BASE::Parse_Options();
} 
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Get_Intersecting_Tetrahedron
//#####################################################################
int Get_Intersecting_Tetrahedron(const PARTICLES<TV>& particles,const TV& location,const TETRAHEDRALIZED_VOLUME<T>& dynamic_volume)
{
    ARRAY<int> intersection_list;dynamic_volume.hierarchy->Intersection_List(location,intersection_list);
    for(int i=1;i<=intersection_list.m;i++) if(TETRAHEDRON<T>(particles.X.Subset(dynamic_volume.mesh.elements(intersection_list(i)))).Inside(location)) return intersection_list(i);
    return 0;
}
//#####################################################################
// Function Create_Duplicate_Mesh
//#####################################################################
void Create_Duplicate_Mesh(const TRIANGLE_MESH& mesh,TRIANGLE_MESH& duplicate_mesh,ARRAY<int>& parent_map) const
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;

    ARRAY<int> mesh_nodes;mesh.elements.Flattened().Get_Unique(mesh_nodes);
    ARRAY<int> child_map(particles.array_collection->Size());
    for(int i=1;i<=mesh_nodes.m;i++){int p=mesh_nodes(i);
        child_map(p)=particles.array_collection->Append(*particles.array_collection,p);}
    for(int i=1;i<=mesh.elements.m;i++) duplicate_mesh.elements.Append(VECTOR<int,3>::Map(child_map,mesh.elements(i)));
    duplicate_mesh.Set_Number_Nodes(particles.array_collection->Size());

    parent_map.Resize(particles.array_collection->Size());for(int p=1;p<=child_map.m;p++) if(child_map(p)) parent_map(child_map(p))=p;
}
//#####################################################################
// Function Add_Binding
//#####################################################################
void Add_Binding(const int particle_index,const ARRAY<int>& parents,const ARRAY<T>& weights,const ARRAY<int>& parent_map)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;

    if(parents.m<1||parents.m>3) PHYSBAM_FATAL_ERROR();
    for(int i=1;i<=parents.m;i++) if(!parent_map(parents(i))) PHYSBAM_FATAL_ERROR();

    if(parents.m==1){
        int dynamic_parent=parent_map(parents(1));
        solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new PARTICLE_BINDING<TV>(particles,particle_index,dynamic_parent));}
    else if(parents.m==2){
        VECTOR<int,2> dynamic_parents(parent_map(parents(1)),parent_map(parents(2)));
        solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new LINEAR_BINDING<TV,2>(particles,particle_index,dynamic_parents,VECTOR<T,2>(weights(1),weights(2))));}
    else{
        VECTOR<int,3> dynamic_parents(parent_map(parents(1)),parent_map(parents(2)),parent_map(parents(3)));
        solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new LINEAR_BINDING<TV,3>(particles,particle_index,dynamic_parents,VECTOR<T,3>(weights(1),weights(2),weights(3))));}
}
//#####################################################################
// Function Create_Deformable_Mattress
//#####################################################################
TETRAHEDRALIZED_VOLUME<T>& Create_Deformable_Mattress()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;

    // 1 read in mattress
    TV_INT counts(31,4,7);TV min_corner(-counts.x*(T).25,0,-counts.z*(T).20),max_corner(counts.x*(T).25,counts.y*(T).5,counts.z*(T).20);
    GRID<TV> mattress_grid(counts,BOX<TV>(min_corner,max_corner));
    TETRAHEDRALIZED_VOLUME<T>& mattress=tests.Create_Mattress(mattress_grid);
    mattress.Update_Number_Nodes();

    // 2 create duplicate mattress surface
    mattress.Initialize_Triangulated_Surface();
    EMBEDDING<TV>& embedding=*new EMBEDDING<TV>(particles);
    ARRAY<int> parent_map;Create_Duplicate_Mesh(mattress.triangulated_surface->mesh,embedding.material_surface_mesh,parent_map);
    delete mattress.triangulated_surface;mattress.triangulated_surface=0;
    redgreen=new RED_GREEN_TRIANGLES<TV>(embedding.material_surface);
    deformable_body_collection.deformable_geometry.Add_Structure(&embedding);

    // refine duplicate mattress surface
    for(int i=1;i<=refinement_level;i++){
        ARRAY<int> top_triangles;
        for(int i=1;i<=embedding.material_surface.mesh.elements.m;i++)
        if(particles.X(embedding.material_surface.mesh.elements(i)(1)).y==(T)2 && particles.X(embedding.material_surface.mesh.elements(i)(2)).y==(T)2) top_triangles.Append(i);
        redgreen->Refine_Simplex_List(top_triangles);}
    embedding.Update_Number_Nodes();embedding.material_surface.mesh.Initialize_Incident_Elements();
    
    // add hard bindings between duplicate surface and mattress surface
    redgreen->Initialize_Segment_Index_From_Midpoint_Index();
    ARRAY<int> refined_particles;embedding.material_surface_mesh.elements.Flattened().Get_Unique(refined_particles);
    ARRAY<int> parents;ARRAY<T> weights;
    for(int i=1;i<=refined_particles.m;i++){int p=refined_particles(i);
        redgreen->Unrefined_Parents(p,parents,weights);
        Add_Binding(p,parents,weights,parent_map);}
    
    // 3 duplicate hard bound surface (soft bound surface)
    tests.Substitute_Soft_Bindings_For_Embedded_Nodes(embedding.material_surface,solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.use_impulses_for_collisions.Fill(false);

    return mattress;
}
//#####################################################################
// Function Create_SCA_Stamp
//#####################################################################
void Create_SCA_Stamp()
{
    // SCA2007 rigid body
    int id=solid_body_collection.rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/Triangulated_Surfaces/sca2007",3.5,true,true);
    RIGID_BODY<TV>& sca=solid_body_collection.rigid_body_collection.Rigid_Body(id);sca.Set_Coefficient_Of_Friction(0);
    // base
    RIGID_BODY<TV>& base=tests.Add_Rigid_Body("plank",base_scale,0);
    base_id=base.particle_index;
    base_initial_frame.r=ROTATION<TV>::From_Euler_Angles(0,(T)pi/2,0);
    // handle
    handle_id=tests.Add_Rigid_Body("cylinder",1,0).particle_index;
    // SCA2007 frame
    base.Update_Bounding_Box();sca.Update_Bounding_Box();
    letters_to_base_frame=FRAME<TV>(TV(0,-.625,0))*FRAME<TV>(TV(),ROTATION<TV>::From_Euler_Angles((T)pi/2,-(T)pi/2,0))*
        FRAME<TV>(base.Axis_Aligned_Bounding_Box().Center()-sca.Axis_Aligned_Bounding_Box().Center());
    sca.X()=letters_to_base_frame.t;sca.Rotation()=letters_to_base_frame.r;
}
//#####################################################################
// Function Create_Rigid_Stamp
//#####################################################################
void Create_Rigid_Stamp()
{
    const std::string siggraph("SIGGRAPH");
    T dx=(T)15/(T)26;
    T widths[8]={3,1,4,4,3,3,3,4};
    T y_rotations[8]={(T)pi/12,0,0,(T)pi/4,0,0,0,0};
    for(unsigned int i=0;i<siggraph.length();i++){ // letters
        std::string letter_filename=STRING_UTILITIES::string_sprintf("Letters/%c",siggraph[i]);
        tests.Add_Rigid_Body(letter_filename,letters_scale,0); 
        T offset=-7.5;for(unsigned int l=0;l<i;l++) offset+=dx*widths[l];
        letters_frames_save(i+1).t.x=offset+dx/(T)2*widths[i];
        letters_frames_save(i+1).r=ROTATION<TV>::From_Euler_Angles((T)pi/2,y_rotations[i],(T)pi);}
    base_id=tests.Add_Rigid_Body("plank",base_scale,0).particle_index; // base
    base_initial_frame.r=ROTATION<TV>::From_Euler_Angles(0,(T)pi/2,0);
    handle_id=tests.Add_Rigid_Body("cylinder",1,0).particle_index; // handle
}
//#####################################################################
// Function Postprocess_Solids_Substep
//#####################################################################
void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE
{
    if(time>final_time) return;

    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    TETRAHEDRALIZED_VOLUME<T>& mattress=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();

    bool initialized_tet_hierarchy=false;
    for(int b=1;b<=soft_bindings.bindings.m;b++){
        int particle_index=soft_bindings.bindings(b).x,parent=soft_bindings.bindings(b).y;
        TV dX=particles.X(parent)-particles.X(particle_index);
        if(dX.Magnitude()>plastic_yield_threshold){ // adjust the parent to lie at the threshold
            if(!initialized_tet_hierarchy){mattress.Initialize_Hierarchy();initialized_tet_hierarchy=true;}
            TV new_parent_position=particles.X(particle_index)+plastic_yield_threshold*dX.Normalized();
            int tet=Get_Intersecting_Tetrahedron(particles,new_parent_position,mattress);
            if(tet==0) continue;
            VECTOR<int,4>& parents=mattress.mesh.elements(tet);
            VECTOR<T,3> weights=TETRAHEDRON<T>::First_Three_Barycentric_Coordinates(new_parent_position,particles.X.Subset(parents));
            if(binding_list.binding_index_from_particle_index.m<parent||!binding_list.binding_index_from_particle_index(parent)) PHYSBAM_FATAL_ERROR();
            int binding_index=binding_list.binding_index_from_particle_index(parent);
            delete binding_list.bindings(binding_index); // substitute a new binding
            binding_list.bindings(binding_index)=new LINEAR_BINDING<TV,4>(particles,parent,parents,weights);
            particles.X(parent)=new_parent_position;}}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;

    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=false;
    solids_parameters.cfl=(T)5;
    last_frame=10000;
    frame_rate=24;
    output_directory="Binding_Plasticity/output";
    output_directory+=(version==1)?"_siggraph":"_sca2007";
    std::cout<<"Frame rate: "<<frame_rate<<std::endl;
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-6;
    solids_parameters.implicit_solve_parameters.cg_iterations=500;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.deformable_object_collision_parameters.collide_with_interior=false;
    
    T dz=(T)3;
    T final_height=2+(T).5*base_thickness-2*stamp_depth;
    T dh=letters_initial_height-final_height;
    translation_curve.Add_Control_Point((T)0,TV(0,letters_initial_height,dz));
    translation_curve.Add_Control_Point((T).2*final_time,TV(0,letters_initial_height,dz/2));
    translation_curve.Add_Control_Point((T).4*final_time,TV(0,letters_initial_height,0));
    translation_curve.Add_Control_Point((T).55*final_time,TV(0,letters_initial_height-(T).25*dh,0));
    translation_curve.Add_Control_Point((T).7*final_time,TV(0,letters_initial_height-(T).5*dh,0));
    translation_curve.Add_Control_Point((T).85*final_time,TV(0,letters_initial_height-(T).75*dh,0));
    translation_curve.Add_Control_Point(final_time,TV(0,final_height,0));
    translation_curve.Add_Control_Point((T)1.15*final_time,TV(0,letters_initial_height-(T).75*dh,0));
    translation_curve.Add_Control_Point((T)1.3*final_time,TV(0,letters_initial_height-(T).5*dh,0));
    translation_curve.Add_Control_Point((T)1.45*final_time,TV(0,letters_initial_height-(T).25*dh,0));
    translation_curve.Add_Control_Point((T)1.6*final_time,TV(0,letters_initial_height,0));
    translation_curve.Add_Control_Point((T)1.8*final_time,TV(0,letters_initial_height,dz/2));
    translation_curve.Add_Control_Point(2*final_time,TV(0,letters_initial_height,dz));
    
    angle_curve.Add_Control_Point((T)0,(T)pi/2);
    angle_curve.Add_Control_Point((T).2*final_time,2*(T)pi/6);
    angle_curve.Add_Control_Point((T).4*final_time,(T)pi/6);
    angle_curve.Add_Control_Point((T).55*final_time,0);
    angle_curve.Add_Control_Point((T).7*final_time,0);
    angle_curve.Add_Control_Point((T).85*final_time,0);
    angle_curve.Add_Control_Point(final_time,0);
    angle_curve.Add_Control_Point((T)1.15*final_time,0);
    angle_curve.Add_Control_Point((T)1.3*final_time,0);
    angle_curve.Add_Control_Point((T)1.45*final_time,0);
    angle_curve.Add_Control_Point((T)1.6*final_time,(T)pi/6);
    angle_curve.Add_Control_Point((T)1.8*final_time,2*(T)pi/6);
    angle_curve.Add_Control_Point(2*final_time,(T)pi/2);
    
    
    TETRAHEDRALIZED_VOLUME<T>& mattress=Create_Deformable_Mattress();
    switch(version){
        case 1: Create_Rigid_Stamp();break;
        case 2: Create_SCA_Stamp();break;
        default: PHYSBAM_FATAL_ERROR();}

    ground_id=tests.Add_Ground().particle_index;

    // correct number nodes
    for(int i=1;i<=deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // collisions
    EMBEDDING<TV>& embedding=deformable_body_collection.deformable_geometry.template Find_Structure<EMBEDDING<TV>&>();
    deformable_body_collection.collisions.collision_structures.Append(&embedding);
    deformable_body_collection.collisions.collision_structures.Append(&mattress);

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
//    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();
    for(int b=1;b<=solid_body_collection.deformable_body_collection.soft_bindings.bindings.m;b++){
        VECTOR<int,2>& binding=solid_body_collection.deformable_body_collection.soft_bindings.bindings(b);
        int mass_scale=1<<2*(refinement_level+1);
        particles.effective_mass(binding.x)*=(T)1/(T)mass_scale;
        particles.one_over_effective_mass(binding.x)*=(T)mass_scale;
        particles.mass(binding.x)=particles.effective_mass(binding.x);particles.one_over_mass(binding.x)=particles.one_over_effective_mass(binding.x);}

    // binding springs
    solid_body_collection.deformable_body_collection.soft_bindings.Initialize_Binding_Mesh();
    spring_force=Create_Edge_Binding_Springs(deformable_body_collection.particles,*solid_body_collection.deformable_body_collection.soft_bindings.binding_mesh,(T)1e5,(T)1);
    solid_body_collection.Add_Force(spring_force);
    solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
    solid_body_collection.Add_Force(Create_Finite_Volume(mattress,new NEO_HOOKEAN<T,3>((T)5e6,(T).45,(T).01,(T).25),true,(T).1));

    // create an elastic force on the drifted surface
    solid_body_collection.Add_Force(Create_Edge_Springs(embedding.material_surface,(T)1e4,(T)1));
    TRIANGLE_BENDING_ELEMENTS<T>* bend=Create_Bending_Elements(embedding.material_surface,(T)1e4);
    bend->Set_Area_Cutoff_From_Triangulated_Surface(embedding.material_surface,(T).1);
    bend->use_force_differential=false;
    solid_body_collection.Add_Force(bend);

    solid_body_collection.Update_Simulated_Particles();

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Initialize_Bodies();

    // set the surface of the embedding geometry to collide with the base and ground only
    deformable_body_collection.collisions.Use_Structure_Collide_Collision_Body();
    HASHTABLE<COLLISION_GEOMETRY_ID,void>& structure_collide_collision_body=deformable_body_collection.collisions.structure_collide_collision_body(1);
    structure_collide_collision_body.Insert(solid_body_collection.rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(ground_id));
    structure_collide_collision_body.Insert(solid_body_collection.rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(base_id));
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{
    return false;
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    T handle_height=2;
    T letters_thickness=(T).2*letters_scale;

    TV translation=translation_curve.Value(time);
    T angle=angle_curve.Value(time);

    ROTATION<TV> rotation=ROTATION<TV>::From_Euler_Angles(angle,0,0);
    TV offset_axis=rotation.Rotate(TV::Axis_Vector(2));

    FRAME<TV> base_frame=FRAME<TV>(translation,rotation)*base_initial_frame;

    if(id==base_id) frame=base_frame; // base
    else if(version==1 && id>=int(1) && id<=int(8)){ // letters
        T letters_offset=(T).5*(base_thickness-letters_thickness)+stamp_depth;
        frame.t=translation-offset_axis*letters_offset;
        frame.r=rotation*letters_frames_save(Value(id)).r;
        frame.t.x+=letters_frames_save(Value(id)).t.x;}
    else if(version==2 && id==int(1)) frame=base_frame*letters_to_base_frame; // letters
    else if(id==handle_id){
        FRAME<TV> handle_to_base_frame(TV(0,(T).5*(base_thickness+handle_height)-(T).1,0));
        frame=base_frame*handle_to_base_frame;} // handle
}
//#####################################################################
};
}
#endif

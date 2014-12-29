//#####################################################################
// Copyright 2007-2008, Craig Schroeder, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
//   1. Single cube fracturing
//   1. Plastic mattress
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Math_Tools/constants.h>
#include <Tools/Matrices/ROTATION.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Deformables/Bindings/LINEAR_BINDING.h>
#include <Deformables/Bindings/RIGID_BODY_BINDING.h>
#include <Deformables/Bindings/SOFT_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <Deformables/Constitutive_Models/ROTATED_LINEAR.h>
#include <Deformables/Constitutive_Models/SIMPLE_PLASTICITY.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/FINITE_VOLUME.h>
#include <Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE.h>
#include <Deformables/Particles/FREE_PARTICLES.h>
#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
#include <Solids/Fracture/FRACTURE_OBJECT.h>
#include <Solids/Fracture/FRACTURE_TETRAHEDRALIZED_VOLUME.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS:public SOLIDS_EXAMPLE<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T_input,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    SOLIDS_STANDARD_TESTS<TV> tests;

    typedef SOLIDS_EXAMPLE<TV> BASE;
    using BASE::solids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::restart;
    using BASE::stream_type;using BASE::solid_body_collection;using BASE::test_number;

    T initial_height;
    ROTATION<TV> initial_orientation;
    TV initial_velocity,initial_angular_velocity;
    GRID<TV> mattress_grid;
    bool cube_mesh;
    bool perform_fracture;
    PLASTICITY_MODEL<T,3>* plasticity_model;
    bool fractured_after_rebuild_topology;
    int substeps_before_rebuild;
    bool push_out;
    FRACTURE_OBJECT<TV,3>* fracture_object;

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type),tests(stream_type,data_directory,solid_body_collection),fracture_object(0)
    {
    }

    ~STANDARD_TESTS()
    {}

    // Unused callbacks
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > F,const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Self_Collisions_Begin_Callback(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Advance_One_Time_Step_End_Callback(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE {}
    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return true;}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    tests.data_directory=data_directory;
    output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d",test_number);
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.cfl=1;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.triangle_collision_parameters.allow_intersections=true;
    solids_parameters.triangle_collision_parameters.allow_intersections_tolerance=(T)1e-7;
    solids_parameters.triangle_collision_parameters.collisions_nonrigid_collision_attempts=3;
    solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=true;
    solids_parameters.use_post_cg_constraints=false;
    solids_parameters.triangle_collision_parameters.output_interaction_pairs=true;
    solids_parameters.rigid_body_collision_parameters.use_push_out=true;
    solids_parameters.use_rigid_deformable_contact=true;
    solids_parameters.rigid_body_collision_parameters.collision_bounding_box_thickness=(T)1e-3;
    solids_parameters.triangle_collision_parameters.collisions_output_number_checked=false;
    solids_parameters.verbose_dt=true;
    solids_parameters.fracture=true;
    solids_parameters.write_static_variables_every_frame=true;
    solids_parameters.triangle_collision_parameters.collisions_final_repulsion_youngs_modulus=(T)20;
    solids_parameters.triangle_collision_parameters.repulsions_youngs_modulus=(T)20;

    switch(test_number){
        case 1: last_frame=240;break;
        case 2: last_frame=240;break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    switch(test_number){
        case 1: Plastic_Mattress(2,2,2,ROTATION<TV>(1,TV(1,1,1)),RANGE<TV>(TV(-(T).5,-(T).5,-(T).5),TV((T).5,(T).5,(T).5)),(T)1.2);break;
        case 2: Plastic_Mattress(2,4,4,ROTATION<TV>(0,TV(1,1,1)),RANGE<TV>(TV(-(T).5,-(T)1.5,-(T)1.5),TV((T).5,(T)1.5,(T)1.5)),(T)1.2);break;
      default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    tests.Add_Ground();

    // correct number nodes
    for(int i=0;i<deformable_body_collection.structures.m;i++) deformable_body_collection.structures(i)->Update_Number_Nodes();

    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);soft_bindings.Set_Mass_From_Effective_Mass();

    // disable strain rate CFL for all forces
    for(int i=0;i<solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->limit_time_step_by_strain_rate=false;

    for(int i=0;i<deformable_body_collection.structures.m;i++) if(!dynamic_cast<SEGMENTED_CURVE<TV>*>(deformable_body_collection.structures(i))){
        deformable_body_collection.collisions.collision_structures.Append(deformable_body_collection.structures(i));
        if(solids_parameters.triangle_collision_parameters.perform_self_collision && (!dynamic_cast<FREE_PARTICLES<TV>*>(deformable_body_collection.structures(i))))
            solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(deformable_body_collection.structures(i));}

    // correct mass
    //binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();

    // add forces
    tests.Add_Gravity();

    //((FRACTURE_EVOLUTION_3D<T>*)solids_parameters.fracture_evolution)->fracture_object=fracture_object;
    //solids_parameters.fracture_evolution->Initialize_Bodies();
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    LOG::cout<<"number of fracture initiation points = "<<fracture_object->number_of_fracture_initiations<<std::endl;
    //solids_parameters.fracture_evolution->Rebuild_Topology();
}
//#####################################################################
// Function Plastic_Mattress
//#####################################################################
void Plastic_Mattress(int nx,int ny,int nz,const ROTATION<TV>& rot,const RANGE<TV>& box,T height)
{
    // initialize forces
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>& embedding=*EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::Create(particles);
    EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& embedded_object=embedding.embedded_object;
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=embedded_object.simplicial_object;

    initial_orientation=rot;
    mattress_grid=GRID<TV>(TV_INT(nx,ny,nz),box);
    initial_height=height;

    deformable_body_collection.Add_Structure(&embedding);
    deformable_body_collection.collisions.collision_structures.Append(&embedding);

    tetrahedralized_volume.Initialize_Cube_Mesh_And_Particles(mattress_grid);
    tetrahedralized_volume.Update_Number_Nodes();
    embedded_object.embedded_particles.Update_Subset_Index_From_Element_Index(); // TODO: make this unnecessary

    tetrahedralized_volume.Update_Bounding_Box();
    TV center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->min_corner.y;
    for(int i=0;i<particles.Size();i++){
        particles.V(i)=initial_velocity+TV::Cross_Product(initial_angular_velocity,particles.X(i)-center);
        particles.X(i)=center+initial_orientation.Rotate(particles.X(i)-center);
        particles.X(i).y+=initial_height-bottom;}
    std::cout << "total vertices = " << particles.Size() << std::endl;
    std::cout << "total tets = " << tetrahedralized_volume.mesh.elements.m << std::endl;
    SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(tetrahedralized_volume,1000);

    if(false) plasticity_model=new SIMPLE_PLASTICITY<T,3>(tetrahedralized_volume.mesh.elements.m,2,2);
    solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new ROTATED_LINEAR<T,3>((T)1e6,(T).3,(T).02),true,(T).1,true,true,false,plasticity_model));

    // initialize_fracture_tetrahedralized_volume
    fracture_object=new FRACTURE_TETRAHEDRALIZED_VOLUME<T>(embedded_object);
    fracture_object->fracture_threshold=TV((T)1e1,(T)1e2,(T)1e3);
    fracture_object->compressive_threshold=TV((T)-1e12,(T)-5e12,(T)-1e12);
    fracture_object->number_of_fracture_initiation_points=4;
    embedded_object.Set_Interpolation_Fraction_Threshold((T)1e-1);
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE
{
    std::string f=STRING_UTILITIES::string_sprintf(".%d",frame);

    // viewer and restart output
    SOLIDS_EXAMPLE<TV>::Write_Output_Files(frame);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/fracture_object"+f,*fracture_object);
    if(plasticity_model) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/fp_inverse"+f,plasticity_model->Fp_inverse);

    // debugging output
    ARRAY<TV> spatial_fracture_bias_direction(fracture_object->embedded_object.simplicial_object.mesh.elements.m);
    //for(int t=0;t<spatial_fracture_bias_direction.m;t++) spatial_fracture_bias_direction(t)=solids_parameters.fracture_evolution->Spatial_Fracture_Bias_Direction(t,(T)1e-4);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/fracture_bias_direction"+f,spatial_fracture_bias_direction);

    if(!fracture_object->embedded_object.embedded_subelements_in_parent_element) fracture_object->embedded_object.Initialize_Embedded_Subelements_In_Parent_Element();
    T one_cut=fracture_object->embedded_object.Fraction_Of_Tetrahedra_With_N_Cuts(1),
        two_cuts=fracture_object->embedded_object.Fraction_Of_Tetrahedra_With_N_Cuts(2),
        three_cuts=fracture_object->embedded_object.Fraction_Of_Tetrahedra_With_N_Cuts(3),
        total=fracture_object->embedded_object.Fraction_Of_Elements_With_Embedded_Subelements();
    LOG::cout<<"PERCENT_OF_BROKEN_ELEMENTS = "<<one_cut<<" "<<two_cuts<<" "<<three_cuts<<", with total = "<<total<<std::endl;

    FINITE_VOLUME<TV,3>& finite_volume=solid_body_collection.template Find_Force<FINITE_VOLUME<TV,3>&>();
    ARRAY<VECTOR<T,2> > min_max_stress_eigenvalues(finite_volume.strain_measure.mesh_object.mesh.elements.m);
    ISOTROPIC_CONSTITUTIVE_MODEL<T,3>& isotropic_model=dynamic_cast<ISOTROPIC_CONSTITUTIVE_MODEL<T,3>&>(finite_volume.constitutive_model);
    for(int t=0;t<finite_volume.strain_measure.mesh_object.mesh.elements.m;t++){
        DIAGONAL_MATRIX<T,3> eigenvalues=isotropic_model.P_From_Strain(finite_volume.Fe_hat(t),(T)1,t); // scale for volume too?
        min_max_stress_eigenvalues(t)=VECTOR<T,2>(eigenvalues.Min(),eigenvalues.Max());}
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/min_max_stress_eigenvalue"+f,min_max_stress_eigenvalues);
}
//#####################################################################
// Function Read_Output_Files_Solids
//#####################################################################
void Read_Output_Files_Solids(const int frame) PHYSBAM_OVERRIDE
{
    std::string f=STRING_UTILITIES::string_sprintf(".%d",frame);

    SOLIDS_EXAMPLE<TV>::Read_Output_Files_Solids(frame);
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/fracture_object"+f,*fracture_object);
    if(plasticity_model) FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/fp_inverse"+f,plasticity_model->Fp_inverse);

    //solids_parameters.fracture_evolution->Reinitialize_Bodies();
    //solids_parameters.fracture_evolution->Initialize_Self_Collision();
}
//#####################################################################
};
}
#endif

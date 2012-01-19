//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Neil Molino, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TET_SIM_FRACTURE_EXAMPLE
//#####################################################################
#ifndef __TET_SIM_FRACTURE_EXAMPLE__
#define __TET_SIM_FRACTURE_EXAMPLE__

#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/FRACTURE_EVOLUTION_3D.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class TET_SIM_FRACTURE_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::output_directory;using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::stream_type;

    SOLIDS_STANDARD_TESTS<TV> tests;

    FRACTURE_OBJECT<TV,3>* fracture_object;  // Use pointer here, so that when we combine with Solids_3D_Example, this is optional -- null if not doing fracture !!!
                                                                      // it is set in the individual examples (i.e., see sphere_drop_example)
    bool perform_fracture;
    T perturb_amount_for_collision_freeness;
    PLASTICITY_MODEL<T,3>* plasticity_model;
    bool fractured_after_rebuild_topology;
    int substeps_before_rebuild;
    bool push_out;

    TET_SIM_FRACTURE_EXAMPLE(const STREAM_TYPE stream_type)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >(stream_type,0,fluids_parameters.NONE),tests(*this),fracture_object(0),perform_fracture(true),perturb_amount_for_collision_freeness((T)1e-6),
        plasticity_model(0),fractured_after_rebuild_topology(false),substeps_before_rebuild(1),push_out(false)
    {
        solids_parameters.fracture=true;
        solids_parameters.write_static_variables_every_frame=true;
        solids_parameters.allow_intersections=true;
        solids_parameters.allow_intersections_tolerance=(T)1e-7;
        solids_parameters.collisions_final_repulsion_youngs_modulus=(T)20;
        solids_parameters.repulsions_youngs_modulus=(T)20;
        solids_parameters.collisions_nonrigid_collision_attempts=3;

        FRACTURE_EVOLUTION_3D<T>* fracture_evolution=new FRACTURE_EVOLUTION_3D<T>(solids_parameters);
        solids_parameters.Set_Fracture_Evolution(fracture_evolution);
    }

    virtual ~TET_SIM_FRACTURE_EXAMPLE()
    {}

    // Unused callbacks
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) {}

//#####################################################################
    virtual bool Add_Scripted_Cuts(const T time){return false;}
//#####################################################################

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    ((FRACTURE_EVOLUTION_3D<T>*)solids_parameters.fracture_evolution)->fracture_object=fracture_object;
    solids_parameters.fracture_evolution->Initialize_Bodies();

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Initialize_Bodies();
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    LOG::cout<<"number of fracture initiation points = "<<fracture_object->number_of_fracture_initiations<<std::endl;
    solids_parameters.fracture_evolution->Rebuild_Topology();
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE
{
    const DEFORMABLE_OBJECT<TV>& deformable_object=solid_body_collection.deformable_object;
    std::string f=STRING_UTILITIES::string_sprintf(".%d",frame);

    // viewer and restart output
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Write_Output_Files(frame);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/fracture_object"+f,*fracture_object);
    if(plasticity_model) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/fp_inverse"+f,plasticity_model->Fp_inverse);

    // debugging output
    ARRAY<TV> spatial_fracture_bias_direction(fracture_object->embedded_object.simplicial_object.mesh.elements.m);
    for(int t=0;t<spatial_fracture_bias_direction.m;t++) spatial_fracture_bias_direction(t)=solids_parameters.fracture_evolution->Spatial_Fracture_Bias_Direction(t,(T)1e-4);
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

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Read_Output_Files_Solids(frame);
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/fracture_object"+f,*fracture_object);
    if(plasticity_model) FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/fp_inverse"+f,plasticity_model->Fp_inverse);

    solids_parameters.fracture_evolution->Reinitialize_Bodies();
    solids_parameters.fracture_evolution->Initialize_Self_Collision();
}
//#####################################################################
};
}
#endif

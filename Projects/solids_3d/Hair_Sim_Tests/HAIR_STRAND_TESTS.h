//#####################################################################
// Copyright 2007, Andrew Selle, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Simulation of Hair Stands
//#####################################################################
#ifndef __HAIR_STRAND_TESTS__
#define __HAIR_STRAND_TESTS__
#include <Tools/Parsing/PARAMETER_LIST.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Deformables/Forces/GUIDE_ADHESION.h>
#include <Deformables/Forces/LINEAR_TET_SPRINGS.h>
#include <Deformables/Forces/SEGMENT_ADHESION.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Solids/Examples_And_Drivers/SOLIDS_DRIVER.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>

namespace PhysBAM{

template<class T_input>
class HAIR_STRAND_TESTS:public SOLIDS_EXAMPLE<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_EXAMPLE<TV> BASE;
    using BASE::solids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::viewer_dir;
    using BASE::stream_type;using BASE::write_substeps_level;using BASE::solid_body_collection;using BASE::test_number;
    using BASE::Set_External_Velocities;using BASE::Set_External_Positions;using BASE::Zero_Out_Enslaved_Velocity_Nodes; // silence -Woverloaded-virtual
    using BASE::user_last_frame;
    
    ARRAY<TV> init_positions_start,init_positions_end;
    ARRAY<int> active_particles;
    ARRAY<int> fixed_nodes_start,fixed_nodes_end;
    SOLIDS_STANDARD_TESTS<TV> tests;
    T wind_start_time,wind_stop_time;
    TETRAHEDRALIZED_VOLUME<T> *volume;
    PARAMETER_LIST parameter_list;
    std::string sim_folder,param_file;
    bool use_implicit;
    // Simulation Parameters
    T cfl_strain_rate;
    T overdamping_fraction;
    T restlength_clamp;
    T edge_stiffness;
    T bending_stiffness;
    T torsion_stiffness;
    T altitude_stiffness;
    // Adhesion
    SEGMENT_ADHESION<TV>* segment_adhesion;
    ARRAY<HAIR_ID> particle_to_spring_id;
    bool use_adhesion;
    T adhesion_stiffness;
    T adhesion_start_radius;
    T adhesion_stop_radius;
    T ether_drag_wind;
    int momentum_conserving_projection_iterations;
    bool use_momentum_conserving_before,use_momentum_conserving_after,use_non_momentum_conserving_before,use_non_momentum_conserving_after;
    //special use flags
    bool reset;

    SEGMENT_MESH project_mesh; // (node_to_change, fixed_node)
    ARRAY<T> project_restlengths; // restlengths of the mesh

//#####################################################################
    HAIR_STRAND_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    void Initialize_Bodies() override;
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override;
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override;
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) override;
    void Update_Time_Varying_Material_Properties(const T time) override;
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) override;
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) override;
    void Add_External_Impulses_Helper(ARRAY_VIEW<TV> V,const T time,const T dt,bool use_momentum_conserving,bool use_non_momentum_conserving);
    void Preprocess_Solids_Substep(const T time) override;    
    void Postprocess_Frame(const int frame) override;
    void Write_Output_Files() const override;
//#####################################################################
};
}
#endif //__HAIR_STRAND_TESTS__

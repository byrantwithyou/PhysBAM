//#####################################################################
// Copyright 2007, Andrew Selle, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Simulation of Hair
//#####################################################################
#ifndef __HAIR_SIM_TESTS__
#define __HAIR_SIM_TESTS__
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
class HAIR_SIM_TESTS:public SOLIDS_EXAMPLE<VECTOR<T_input,3> >,TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<VECTOR<T_input,3> >::MASS_MODIFIER
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    typedef SOLIDS_EXAMPLE<TV> BASE;
    using BASE::solids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::viewer_dir;
    using BASE::stream_type;using BASE::restart;using BASE::solid_body_collection;using BASE::test_number;
    using BASE::Set_External_Velocities;using BASE::Set_External_Positions;using BASE::Zero_Out_Enslaved_Velocity_Nodes; // silence -Woverloaded-virtual
    using BASE::user_last_frame;
protected:
    using BASE::write_substeps_level;
public:
    struct COLLISION_PAIR_COMPARATOR{
        COLLISION_PAIR_COMPARATOR(TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>* parameter,RIGID_BODY<TV>* rigid_body){
            geometry=parameter;
            head=rigid_body;
        }

        bool operator()(const VECTOR<int,4>& pair1, const VECTOR<int,4>& pair2){
            const ARRAY<TV>& X=geometry->X_self_collision_free;
            T min1=head->implicit_object->Signed_Distance(X(pair1[0]));
            T min2=head->implicit_object->Signed_Distance(X(pair2[0]));
            for(int i=1;i<4;i++) {
                min1=min(min1,head->implicit_object->Signed_Distance(X(pair1[i])));
                min2=min(min2,head->implicit_object->Signed_Distance(X(pair2[i])));
            }
            return (min1<min2);}
        
        TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>* geometry;
        RIGID_BODY<TV>* head;
    };

    const T start_time;
    ARRAY<int> active_particles;
    FRAME<TV> init_frame;
    FRAME<TV> init_pick;
    ARRAY<TV> bindings1,bindings2;
    ARRAY<TV> binding_velocities;
    T max_binding_speed_squared;
    ARRAY<int> fixed_nodes_indices; // dummyid -> fixed_nodeid 
    ARRAY<int> fixed_nodes_to_apply;
    ARRAY<TV> bindings1_to_apply,bindings2_to_apply;
    ARRAY<TV> interp_points;
    ARRAY<TV> interp_normals;
    SOLIDS_STANDARD_TESTS<TV> tests;
    T wind_start_time,wind_stop_time;
    int sphere_id;
    RIGID_BODY<TV>* head;
    TETRAHEDRALIZED_VOLUME<T> *volume,*guide_volume,*sim_guide_volume;
    SEGMENT_ADHESION<TV>* segment_adhesion;
    GUIDE_ADHESION<TV>* guide_adhesion;
    DEFORMABLE_BODY_COLLECTION<TV> *guide_object1,*guide_object2;
    VIEWER_DIR guide_viewer_dir{"."};
    PARAMETER_LIST parameter_list;
    ARRAY<HAIR_ID> particle_to_spring_id;
    ARRAY<int,HAIR_ID> spring_id_to_particle; // single particle representative of hair
    ARRAY<ARRAY<int>,PARTITION_ID> partition_spring_representative; // for each partition a list of particles that represent the spring
    COLLISION_PAIR_COMPARATOR *comparator;
    RIGID_BODY<TV> *implicit_rigid_body;
    std::string sim_folder,rigid_model,param_file;
    int offset;
    int current_frame,current_levelset;
    T levelset_frequency;
    bool use_implicit;
    // Simulation Parameters
    bool use_wind, use_drag;
    T drag_viscosity;
    T cfl_strain_rate;
    T overdamping_fraction;
    T restlength_clamp;
    T edge_stiffness;
    T bending_stiffness;
    T torsion_stiffness;
    T altitude_stiffness;
    bool use_adhesion;
    T adhesion_stiffness;
    T adhesion_start_radius;
    T adhesion_stop_radius;
    bool use_guide;
    bool use_spring_guide;
    T guide_edge_stiffness;
    T guide_altitude_stiffness;
    T guide_stiffness;
    int max_connections;
    T guide_thickness;
    bool use_deforming_levelsets;
    bool use_collisions_mass_modify;
    bool use_progressive_collision_thickness;
    int momentum_conserving_projection_iterations;
    bool use_momentum_conserving_before,use_momentum_conserving_after,use_non_momentum_conserving_before,use_non_momentum_conserving_after;
    bool project_second_node_to_surface;
    bool use_eulerian_level_set_interpolation;
    T head_friction;
    //mass_alterations
    VECTOR<T,4> saved;
    //epoch data
    HAIR_ID number_of_hairs;
    std::string cameras;

    SEGMENT_MESH project_mesh; // (node_to_change, fixed_node)
    ARRAY<T> project_restlengths; // restlengths of the mesh

    ARRAY<int> distance_to_root;
    ARRAY<T> collision_tolerances;

//#####################################################################
    HAIR_SIM_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    void Initialize_Bodies() override;
    void Update_Keyframed_Parameters_For_Time_Update(const T time);
    template<class T_IMPLICIT_COMBINED> void Update_Keyframed_Parameters_For_Time_Update_Helper(const T time,T_IMPLICIT_COMBINED& combined);
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override;
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override;
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) override;
    void Preprocess_Solids_Substep(const T time) override;    
    void Update_Time_Varying_Material_Properties(const T time) override;
    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) override;
    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) override;
    void Postprocess_Frame(const int frame) override;
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) override;
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) override;
    void Add_External_Impulses_Helper(ARRAY_VIEW<TV> V,const T time,const T dt,bool use_momentum_conserving,bool use_non_momentum_conserving);
    void Write_Output_Files() const override;
    template<class T_IMPLICIT_COMBINED> void Write_Interpolated_Level_Set(const VIEWER_DIR& viewer_dir,T_IMPLICIT_COMBINED& combined) const;
    void Point_Face_Mass(const T attempt_ratio,const VECTOR<int,4>& nodes,const VECTOR<T,3>& weights,T& one_over_mass1,T& one_over_mass2,T& one_over_mass3,T& one_over_mass4);
    void Point_Face_Mass(const T attempt_ratio,const VECTOR<int,4>& nodes,const VECTOR<T,3>& weights,VECTOR<T,4>& one_over_mass) override;
    void Point_Face_Mass(const T attempt_ratio,const VECTOR<int,4>& nodes,const VECTOR<T,3>& weights,ARRAY_VIEW<T>& one_over_mass) override;
    void Point_Face_Mass_Revert(const VECTOR<int,4>& nodes,ARRAY_VIEW<T>& one_over_mass) override;
    void Edge_Edge_Mass(const T attempt_ratio,const VECTOR<int,4>& nodes,const VECTOR<T,2>& weights,T& one_over_mass1,T& one_over_mass2,T& one_over_mass3,T& one_over_mass4);
    void Edge_Edge_Mass(const T attempt_ratio,const VECTOR<int,4>& nodes,const VECTOR<T,2>& weights,VECTOR<T,4>& one_over_mass) override;
    void Edge_Edge_Mass(const T attempt_ratio,const VECTOR<int,4>& nodes,const VECTOR<T,2>& weights,ARRAY_VIEW<T>& one_over_mass) override;
    void Edge_Edge_Mass_Revert(const VECTOR<int,4>& nodes,ARRAY_VIEW<T>& one_over_mass) override;
    void Mass_Revert(const VECTOR<int,4>& nodes,ARRAY_VIEW<T>& one_over_mass);
    void Reorder_Pairs(ARRAY<VECTOR<int,4> >& edge_edge_pairs,ARRAY<VECTOR<int,4> >& point_face_pairs) override;
    void Compute_Binding_Velocities();
//#####################################################################
};
}
#endif //__HAIR_SIM_TESTS__

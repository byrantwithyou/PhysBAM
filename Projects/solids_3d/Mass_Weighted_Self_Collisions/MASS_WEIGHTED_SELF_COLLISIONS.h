//#####################################################################
// Copyright 2007, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MASS_WEIGHTED_SELF_COLLISIONS__
#define __MASS_WEIGHTED_SELF_COLLISIONS__
#include <Tools/Parsing/PARAMETER_LIST.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <Deformables/Forces/LINEAR_TET_SPRINGS.h>
#include <Deformables/Forces/SEGMENT_ADHESION.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Solids/Examples_And_Drivers/SOLIDS_DRIVER.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>

namespace PhysBAM{

template<class T_input>
class MASS_WEIGHTED_SELF_COLLISIONS:public SOLIDS_EXAMPLE<VECTOR<T_input,3> >,TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<VECTOR<T_input,3> >::MASS_MODIFIER
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_EXAMPLE<TV> BASE;
    using BASE::solids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;
    using BASE::stream_type;using BASE::solid_body_collection;using BASE::parse_args;using BASE::test_number;
    using BASE::Set_External_Velocities;using BASE::Set_External_Positions;using BASE::Zero_Out_Enslaved_Velocity_Nodes; // silence -Woverloaded-virtual

    struct COLLISION_PAIR_COMPARATOR{
        COLLISION_PAIR_COMPARATOR(TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>* parameter) {
            geometry=parameter;
        }

        bool operator()(const VECTOR<int,4>& pair1, const VECTOR<int,4>& pair2) {
            const ARRAY<TV>& X=geometry->X_self_collision_free;
            T min1=X(pair1[0]).y;
            T min2=X(pair2[0]).y;
            for(int i=1;i<4;i++) {
                min1=min(min1,X(pair1[i]).y);
                min2=min(min2,X(pair2[i]).y);
            }
            return (min1<min2);}

        TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>* geometry;
    };

    VECTOR<T,4> saved_mass;
    SOLIDS_STANDARD_TESTS<TV> tests;
    COLLISION_PAIR_COMPARATOR *comparator;

//#####################################################################
    MASS_WEIGHTED_SELF_COLLISIONS(const STREAM_TYPE stream_type);
    void Initialize_Bodies() PHYSBAM_OVERRIDE;
    // overrides from MASS_MODIFIER
    void Point_Face_Mass(const T attempt_ratio,const VECTOR<int,4>& nodes,const VECTOR<T,3>& weights,VECTOR<T,4>& one_over_mass);
    void Point_Face_Mass(const T attempt_ratio,const VECTOR<int,4>& nodes,const VECTOR<T,3>& weights,ARRAY_VIEW<T>& one_over_mass);
    void Point_Face_Mass_Revert(const VECTOR<int,4>& nodes,ARRAY_VIEW<T>& one_over_mass);
    void Edge_Edge_Mass(const T attempt_ratio,const VECTOR<int,4>& nodes,const VECTOR<T,2>& weights,VECTOR<T,4>& one_over_mass);
    void Edge_Edge_Mass(const T attempt_ratio,const VECTOR<int,4>& nodes,const VECTOR<T,2>& weights,ARRAY_VIEW<T>& one_over_mass);
    void Edge_Edge_Mass_Revert(const VECTOR<int,4>& nodes,ARRAY_VIEW<T>& one_over_mass);
    void Reorder_Pairs(ARRAY<VECTOR<int,4> >& edge_edge_pairs,ARRAY<VECTOR<int,4> >& point_face_pairs);
    void Register_Options() PHYSBAM_OVERRIDE;
    void Parse_Options() PHYSBAM_OVERRIDE;
    void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
};
}
#endif

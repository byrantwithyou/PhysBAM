//#####################################################################
// Copyright 2018, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PENALTY_FORCE_COLLECTION__
#define __PENALTY_FORCE_COLLECTION__
#include <Core/Data_Structures/CHAINED_ARRAY.h>
#include <Grid_Tools/Grids/GRID.h>
namespace PhysBAM{
template<class TV> class RIGID_DEFORMABLE_PENALTY_WITH_FRICTION;
template<class TV> class IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION;
template<class TV> class RIGID_PENALTY_WITH_FRICTION;
template<class TV> class SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION;
template<class TV> class SOLID_BODY_COLLECTION;
template<class TV> class MOVE_RIGID_BODY_DIFF;
template<class TV> class IMPLICIT_OBJECT;
template<class TV> class TRIANGLE_COLLISION_PARAMETERS;

template<class TV>
class PENALTY_FORCE_COLLECTION
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    GRID<TV> grid;
    int max_resolution=20,max_cells=10000,max_cells_per_object=250;
    SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    const ARRAY<int>& simulated_particles;
    const ARRAY<MOVE_RIGID_BODY_DIFF<TV> >& move_rb_diff;

    CHAINED_ARRAY<int,TV_INT> cell_particles;
    CHAINED_ARRAY<int,TV_INT> cell_objects;
    CHAINED_ARRAY<PAIR<int,int>,TV_INT> cell_vertices; // <rigid_body,vertex>

    // First entry MUST be int, and its value MUST be nonnegative.
    struct RASTERIZED_DATA
    {
        int id;
        T phi;
    };
    CHAINED_ARRAY<RASTERIZED_DATA,TV_INT> rasterized_data;

    // Do not rasterize particles on these rigid bodies (eg, ground)
    HASHTABLE<int> exclude_rigid_body_simplices;

    IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>* di_penalty=0;
    RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>* rd_penalty=0;
    RIGID_PENALTY_WITH_FRICTION<TV>* rr_penalty=0;
    SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>* dd_penalty=0;

    T const_repulsion_thickness=(T).01;
    ARRAY<T> repulsion_thickness; // must be same size as particles.number
    ARRAY<bool> recently_modified; // must be same size as particles.number

    PENALTY_FORCE_COLLECTION(SOLID_BODY_COLLECTION<TV>& solid_body_collection,
        const ARRAY<int>& simulated_particles,
        const ARRAY<MOVE_RIGID_BODY_DIFF<TV> >& move_rb_diff)
        :solid_body_collection(solid_body_collection),
        simulated_particles(simulated_particles),move_rb_diff(move_rb_diff)
    {grid.domain=RANGE<TV>::Empty_Box();}

    void Init(T stiffness,T friction,TRIANGLE_COLLISION_PARAMETERS<TV>* param,
        bool use_di,bool use_dd,bool use_rd,bool use_rr);
    
    PENALTY_FORCE_COLLECTION(const PENALTY_FORCE_COLLECTION&) = delete;
    void operator=(const PENALTY_FORCE_COLLECTION&) = delete;
    ~PENALTY_FORCE_COLLECTION() = default;

    void Update_Collision_Detection_Structures();
    void Update_Cell_Vertices(bool new_grid);
    void Update_Rasterized_Data(bool new_grid);
    void Update_Cell_Particles(bool new_grid);
    void Update_Cell_Objects(bool new_grid);
    void Get_DI_Collision_Candidates();
    void Get_DD_Collision_Candidates();
    void Get_RD_Collision_Candidates();
    void Get_RR_Collision_Candidates();
    void Save_State();
    void Update_Attachments_And_Prune_Pairs();
    bool Update_Grid();
//#####################################################################
};
}
#endif

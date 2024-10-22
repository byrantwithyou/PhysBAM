//#####################################################################
// Copyright 2006-2007, Jon Gretarsson, Avi Robinson-Mosher, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SMOKE_STANDARD_TESTS_2D
//#####################################################################
// Provides 3 tests that are consistent across all smoke codes:
//   1. Plume
//   2. Plume past circle
//   3. Explosion
//   4. Vortex in a box test
// Also supports a variety of resolutions.
//#####################################################################
#ifndef __SMOKE_STANDARD_TESTS_2D__
#define __SMOKE_STANDARD_TESTS_2D__

#include <Core/Math_Tools/RANGE.h>
#include <Core/Matrices/FRAME.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
namespace PhysBAM{

template<class TV> class SOLIDS_FLUIDS_EXAMPLE;
template<class TV> class RIGID_BODY_COLLECTION;
template<class TV> class FLUIDS_PARAMETERS_UNIFORM;
template<class TV> class INCOMPRESSIBLE_FLUID_COLLECTION;
template<class TV> class PROJECTION_DYNAMICS_UNIFORM;
template<class TV>
class SMOKE_STANDARD_TESTS_2D
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    SOLIDS_FLUIDS_EXAMPLE<TV>& example;
    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters;
    INCOMPRESSIBLE_FLUID_COLLECTION<TV>& incompressible_fluid_collection;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;

    int test_number;
    GRID<TV> grid;
    RANGE<TV> source;
    MATRIX<T,3> world_to_source;
    VECTOR<T,2> source_velocity;
    T rho;
    T explosion_divergence,explosion_end_time;
    int oriented_box;
    T rotation_angle;
    FRAME<TV> rotation_frame;
    int left_box,right_box;
    ARRAY<T,FACE_INDEX<TV::m> > beta_face;
    ARRAY<T,FACE_INDEX<TV::m> > divergence_face_weights;
    
    SMOKE_STANDARD_TESTS_2D(SOLIDS_FLUIDS_EXAMPLE<TV>& example,FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters,INCOMPRESSIBLE_FLUID_COLLECTION<TV>& incompressible_fluid_collection,
        RIGID_BODY_COLLECTION<TV>& rigid_body_collection);
    virtual ~SMOKE_STANDARD_TESTS_2D();

//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
void Set_Dirichlet_Boundary_Conditions(const T time)
{
    //TODO This doesn't do anything.  To remove it, we'd need to change something both in Projects and in Public_Library
}
//#####################################################################
    virtual VECTOR<T,2> Initial_Velocity(const VECTOR<T,2>& X) const;
    void Initialize_Bodies();
    void Initialize(const int test_number_input,const int resolution,const T angle_fraction=0);
    void Get_Divergence(ARRAY<T,VECTOR<int,2> >& divergence,const T dt,const T time);
    void Get_Object_Velocities(PROJECTION_DYNAMICS_UNIFORM<TV>& projection,const T dt,const T time);
//#####################################################################
};
}
#endif

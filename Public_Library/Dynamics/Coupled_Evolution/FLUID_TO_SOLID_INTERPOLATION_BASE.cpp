//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_TO_SOLID_INTERPOLATION_BASE
//##################################################################### 
#include <Tools/Arrays/CONSTANT_ARRAY.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY_ID.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Dynamics/Coupled_Evolution/FLUID_TO_SOLID_INTERPOLATION_BASE.h>
#include <Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED.h>
#include <Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUID_TO_SOLID_INTERPOLATION_BASE<TV>::
FLUID_TO_SOLID_INTERPOLATION_BASE(const COLLISION_AWARE_INDEX_MAP<TV>& map)
    :index_map(map)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUID_TO_SOLID_INTERPOLATION_BASE<TV>::
~FLUID_TO_SOLID_INTERPOLATION_BASE()
{
}
//#####################################################################
// Function Times
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION_BASE<TV>::
Times(const ARRAY<T>& fluid_velocity,GENERALIZED_VELOCITY<TV>& solid_velocity) const
{
    solid_velocity*=(T)0;
    Times_Add(fluid_velocity,solid_velocity);
}
//#####################################################################
// Function Transpose_Times
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION_BASE<TV>::
Transpose_Times(const GENERALIZED_VELOCITY<TV>& solid_force,ARRAY<T>& fluid_force) const
{
    // TODO: Careful to zero out enough of the solids state.
    fluid_force.Fill(0);
    Transpose_Times_Add(solid_force,fluid_force);
}
//#####################################################################
// Function Test_Matrix
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION_BASE<TV>::
Test_Matrix(int number_fluid_faces,int number_particles,int number_rigid_particles) const
{
    RANDOM_NUMBERS<T> random;

    ARRAY<TV> V(number_particles),V1(number_particles);
    random.Fill_Uniform(V,-1,1);

    ARRAY<TWIST<TV> > twist(number_rigid_particles),twist2(number_rigid_particles);
    random.Fill_Uniform(twist,-1,1);

    ARRAY<T> U(number_fluid_faces),U1(number_fluid_faces);
    random.Fill_Uniform(U,-1,1);

    ARRAY<int> empty;
    GENERALIZED_VELOCITY<TV> solids(V,empty,twist,empty,empty),solids2(V1,empty,twist2,empty,empty);

    Times(U,solids2);
    Transpose_Times(solids,U1);

    CONSTANT_ARRAY<RIGID_BODY_MASS<TV,true> > rigid_mass(twist.m,RIGID_BODY_MASS<TV,true>(1,DIAGONAL_MATRIX<T,TV::SPIN::m>::Identity_Matrix()));
    T inner_solids=V.Dot(V1)+twist.Inner_Product(rigid_mass,twist2);
    T inner_fluids=U.Dot_Product(U,U1);

    LOG::cout<<"FLUID_TO_SOLID_INTERPOLATION_BASE Test: "<<inner_solids<<"  vs  "<<inner_fluids<<"  relative  "<<
        abs(inner_solids-inner_fluids)/maxabs((T)1e-30,inner_solids,inner_fluids)<<std::endl;
}
//#####################################################################
// Function Store_Maps
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION_BASE<TV>::
Store_Maps(const GENERALIZED_VELOCITY<TV>& G)
{
    V_size=G.V.array.Size();
    V_indices=&G.V.indices;
}
//#####################################################################
namespace PhysBAM{
template class FLUID_TO_SOLID_INTERPOLATION_BASE<VECTOR<float,1> >;
template class FLUID_TO_SOLID_INTERPOLATION_BASE<VECTOR<float,2> >;
template class FLUID_TO_SOLID_INTERPOLATION_BASE<VECTOR<float,3> >;
template class FLUID_TO_SOLID_INTERPOLATION_BASE<VECTOR<double,1> >;
template class FLUID_TO_SOLID_INTERPOLATION_BASE<VECTOR<double,2> >;
template class FLUID_TO_SOLID_INTERPOLATION_BASE<VECTOR<double,3> >;
}

//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED
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
#include <Solids/Solids_Evolution/BACKWARD_EULER_SYSTEM.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Dynamics/Coupled_Evolution/MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED.h>
#include <Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED.h>
#include <Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<TV>::
MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED(const UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& info)
    :MATRIX_SOLID_INTERPOLATION_BASE<TV>(info)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<TV>::
~MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED()
{
}
//#####################################################################
// Function Number_Of_Constraints
//#####################################################################
template<class TV> COUPLING_CONSTRAINT_ID MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<TV>::
Number_Of_Constraints() const
{
    return COUPLING_CONSTRAINT_ID(entries.m*TV::m);
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<TV>::
Compute(const int ghost_cells)
{
}
//#####################################################################
// Function Times_Add
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<TV>::
Times_Add(const GENERALIZED_VELOCITY<TV>& solids,ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints) const
{
    for(int i=0;i<entries.m;i++) for(int axis=0;axis<TV::m;axis++)
        constraints(COUPLING_CONSTRAINT_ID(i*TV::m+axis))+=solids.V.array(entries(i))(axis);
}
//#####################################################################
// Function Transpose_Times_Add
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<TV>::
Transpose_Times_Add(const ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints,GENERALIZED_VELOCITY<TV>& solids) const
{
    for(int i=0;i<entries.m;i++) for(int axis=0;axis<TV::m;axis++)
        solids.V.array(entries(i))(axis)+=constraints(COUPLING_CONSTRAINT_ID(i*TV::m+axis));
}
//#####################################################################
// Function Print_Each_Matrix
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<TV>::
Print_Each_Matrix(int n,GENERALIZED_VELOCITY<TV>& G) const
{
    OCTAVE_OUTPUT<T> oo(STRING_UTILITIES::string_sprintf("J-%i.txt",n).c_str());
    oo.Begin_Sparse_Matrix("J",Value(Number_Of_Constraints()),G.Raw_Size());
    ARRAY<int> reverse_map_deformable(G.V.array.Size());
    reverse_map_deformable.Subset(G.V.indices)=IDENTITY_ARRAY<>(G.V.Size());

    for(int i=0;i<entries.m;i++)
        for(int axis=0;axis<TV::m;axis++)
            oo.Add_Sparse_Entry(i*TV::m+axis,(reverse_map_deformable(entries(i))-1)*TV::m+axis,1);

    oo.End_Sparse_Matrix();
}
//#####################################################################
// Function Add_Raw_Matrix
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<TV>::
Add_Raw_Matrix(ARRAY<TRIPLE<int,int,T> >& data) const
{
    ARRAY<int> reverse_map_deformable(this->V_size);
    reverse_map_deformable.Subset(*this->V_indices)=IDENTITY_ARRAY<>(this->V_size);

    for(int i=0;i<entries.m;i++)
        for(int axis=0;axis<TV::m;axis++)
            data.Append(TRIPLE<int,int,T>(i*TV::m+axis,(reverse_map_deformable(entries(i))-1)*TV::m+axis,1));
}
//#####################################################################
// Function Add_Diagonal
//#####################################################################
template<class TV> void MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<TV>::
Add_Diagonal(ARRAY<T,COUPLING_CONSTRAINT_ID>& diagonal,const GENERALIZED_MASS<TV>& solid_mass) const
{
    for(int i=0;i<entries.m;i++) for(int axis=0;axis<TV::m;axis++)
        diagonal(COUPLING_CONSTRAINT_ID(i*TV::m+axis))+=solid_mass.one_over_mass.array(entries(i));
}
//#####################################################################
namespace PhysBAM{
template class MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<VECTOR<float,1> >;
template class MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<VECTOR<float,2> >;
template class MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<VECTOR<float,3> >;
template class MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<VECTOR<double,1> >;
template class MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<VECTOR<double,2> >;
template class MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED<VECTOR<double,3> >;
}

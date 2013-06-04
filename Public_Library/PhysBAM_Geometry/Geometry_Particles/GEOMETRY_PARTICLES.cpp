//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Particles/PARTICLES.h>
#include <PhysBAM_Tools/Particles/PARTICLES_FORWARD.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GEOMETRY_PARTICLES<TV>::
GEOMETRY_PARTICLES()
    :X(0,0),V(0,0),store_velocity(false)
{
    Add_Array(ATTRIBUTE_ID_X,&X);
}
//#####################################################################
// Destructor 
//#####################################################################
template<class TV> GEOMETRY_PARTICLES<TV>::
~GEOMETRY_PARTICLES()
{}
//#####################################################################
// Function Initialize_Geometry_Particle
//#####################################################################
static int Initialize_Geometry_Particle()
{
    Register_Attribute_Name(ATTRIBUTE_ID_FRAME,"frame");
    Register_Attribute_Name(ATTRIBUTE_ID_TWIST,"twist");
    Register_Attribute_Name(ATTRIBUTE_ID_RIGID_BODY,"rigid_body");
    Register_Attribute_Name(ATTRIBUTE_ID_X,"X");
    Register_Attribute_Name(ATTRIBUTE_ID_V,"V");
    Register_Attribute_Name(ATTRIBUTE_ID_STRUCTURE_IDS,"structure_ids");
    Register_Attribute_Name(ATTRIBUTE_ID_ID,"id");
    Register_Attribute_Name(ATTRIBUTE_ID_COLOR,"color");
    Register_Attribute_Name(ATTRIBUTE_ID_RADIUS,"radius");
    Register_Attribute_Name(ATTRIBUTE_ID_DISPLAY_SIZE,"display_size");

    #define READ_WRITE_VECTOR_HELPER(T,d) \
        Register_Attribute_Sample<VECTOR<T,d> >();         \
        Register_Attribute_Sample<FRAME<VECTOR<T,d> > >(); \
        Register_Attribute_Sample<TWIST<VECTOR<T,d> > >();

    #define READ_WRITE_SCALAR_HELPER(T) \
        Register_Attribute_Sample<T>(); \
        Register_Attribute_Sample<VECTOR<T,0> >(); \
        Register_Attribute_Sample<DIAGONAL_MATRIX<T,1> >(); \
        Register_Attribute_Sample<DIAGONAL_MATRIX<T,2> >(); \
        Register_Attribute_Sample<DIAGONAL_MATRIX<T,3> >(); \
        Register_Attribute_Sample<MATRIX<T,1,1> >(); \
        Register_Attribute_Sample<MATRIX<T,0,0> >(); \
        READ_WRITE_VECTOR_HELPER(T,1);READ_WRITE_VECTOR_HELPER(T,2);READ_WRITE_VECTOR_HELPER(T,3);

    Register_Attribute_Sample<VECTOR<int,1> >();
    Register_Attribute_Sample<VECTOR<int,2> >();
    Register_Attribute_Sample<VECTOR<int,3> >();
    Register_Attribute_Sample<int>();
    Register_Attribute_Sample<bool>();
    Register_Attribute_Sample<unsigned short>();

    READ_WRITE_SCALAR_HELPER(float);
    READ_WRITE_SCALAR_HELPER(double);

    return 1;
}
int initialize_geometry_particle=Initialize_Geometry_Particle();
//#####################################################################
template class GEOMETRY_PARTICLES<VECTOR<float,1> >;
template class GEOMETRY_PARTICLES<VECTOR<float,2> >;
template class GEOMETRY_PARTICLES<VECTOR<float,3> >;
template class GEOMETRY_PARTICLES<VECTOR<double,1> >;
template class GEOMETRY_PARTICLES<VECTOR<double,2> >;
template class GEOMETRY_PARTICLES<VECTOR<double,3> >;
}

//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/FRAME.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Particles/PARTICLES.h>
#include <Tools/Particles/PARTICLES_FORWARD.h>
#include <Tools/Vectors/TWIST.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
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

    #define READ_WRITE_VECTOR_HELPER(d) \
        Register_Attribute_Sample<VECTOR<float,d> ,VECTOR<double,d> >();         \
        Register_Attribute_Sample<FRAME<VECTOR<float,d> > ,FRAME<VECTOR<double,d> > >(); \
        Register_Attribute_Sample<TWIST<VECTOR<float,d> > ,TWIST<VECTOR<double,d> > >();

    Register_Attribute_Sample<VECTOR<int,1> >();
    Register_Attribute_Sample<VECTOR<int,2> >();
    Register_Attribute_Sample<VECTOR<int,3> >();
    Register_Attribute_Sample<int>();
    Register_Attribute_Sample<bool>();
    Register_Attribute_Sample<unsigned short>();
    Register_Attribute_Sample<float,double>();
    Register_Attribute_Sample<VECTOR<float,0> ,VECTOR<double,0> >();
    Register_Attribute_Sample<DIAGONAL_MATRIX<float,1> ,DIAGONAL_MATRIX<double,1> >();
    Register_Attribute_Sample<DIAGONAL_MATRIX<float,2> ,DIAGONAL_MATRIX<double,2> >();
    Register_Attribute_Sample<DIAGONAL_MATRIX<float,3> ,DIAGONAL_MATRIX<double,3> >();
    Register_Attribute_Sample<MATRIX<float,1,1> ,MATRIX<double,1,1> >();
    Register_Attribute_Sample<MATRIX<float,0,0> ,MATRIX<double,0,0> >();
    READ_WRITE_VECTOR_HELPER(1);READ_WRITE_VECTOR_HELPER(2);READ_WRITE_VECTOR_HELPER(3);

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

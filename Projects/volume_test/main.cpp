//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_ORIGIN_AREAS.h>
#include <PhysBAM_Dynamics/Read_Write/EPS_FILE_GEOMETRY.h>
#include <iomanip>
using namespace PhysBAM;

ARRAY<int> trap_cases;

RANDOM_NUMBERS<double> rn;

void Case_Test()
{
    typedef double T;
    typedef VECTOR<double,3> TV;

    TV a,b,c,d,e,f;
    rn.Fill_Uniform(a,-1,1);
    rn.Fill_Uniform(b,-1,1);
    rn.Fill_Uniform(c,-1,1);
    rn.Fill_Uniform(d,-1,1);
    rn.Fill_Uniform(e,-1,1);
    rn.Fill_Uniform(f,-1,1);

    TRIANGLE_ORIGIN_AREAS::DATA<T,1,6> data;
    TRIANGLE_ORIGIN_AREAS::Volume_From_Triangles(data,a,b,c,d,e,f);
}

int main(int argc,char *argv[])
{
    typedef double T;
    typedef double RW;
    typedef VECTOR<T,3> TV;
    LOG::cout<<std::setprecision(16);

    for(int i=0;i<100000000;i++)
    {
        try{Case_Test();}catch(...){}
    }

    return 0;
}

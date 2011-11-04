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
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_ORIGIN_AREAS.h>
#include <PhysBAM_Dynamics/Read_Write/EPS_FILE_GEOMETRY.h>
#include <iomanip>
using namespace PhysBAM;

ARRAY<int> trap_cases;
typedef double T;
typedef VECTOR<double,3> TV;

RANDOM_NUMBERS<double> rn;
int N=200;

T Approximate_Area(TV A,TV B,TV C,TV D,TV E,TV F)
{
    int in=0;
    int sign=1;
    if(TV::Triple_Product(A,B,C)<0) exchange(A,B),sign=-sign;
    if(TV::Triple_Product(D,E,F)<0) exchange(D,E),sign=-sign;

    for(int x=0;x<N;x++)
        for(int y=0;y<N;y++)
            for(int z=0;z<N;z++)
            {
                TV p=TV(x,y,z)/N*2-1;
                if(TETRAHEDRON<T>::Barycentric_Inside(p,TV(),A,B,C) && TETRAHEDRON<T>::Barycentric_Inside(p,TV(),D,E,F)) in++;
            }
    return (T)sign*in/N/N/N*8;
}

void Case_Test()
{
    TV a,b,c,d,e,f;
    rn.Fill_Uniform(a,-1,1);
    rn.Fill_Uniform(b,-1,1);
    rn.Fill_Uniform(c,-1,1);
    rn.Fill_Uniform(d,-1,1);
    rn.Fill_Uniform(e,-1,1);
    rn.Fill_Uniform(f,-1,1);

    TRIANGLE_ORIGIN_AREAS::DATA<T,1,6> data;
    TRIANGLE_ORIGIN_AREAS::Volume_From_Triangles(data,a,b,c,d,e,f);

    
    if(!data.V.x) return;
    T approx=Approximate_Area(a,b,c,d,e,f);

//    T approx=Approximate_Area(TV(1,0,0),TV(0,1,0),TV(0,0,1),TV(1,0,0),TV(0,0,1),TV(0,1,0));

    printf("VOLUMES: %9.6f %9.6f (%.6f)\n", data.V.x, approx, fabs(data.V.x-approx));
}

int main(int argc,char *argv[])
{
    typedef double T;
    typedef double RW;
    typedef VECTOR<T,3> TV;
    LOG::cout<<std::setprecision(16);


    for(int i=0;i<500;i++)
    {
        try{Case_Test();}catch(...){}
    }

    return 0;
}

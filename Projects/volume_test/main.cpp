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

T Approximate_Volume(TV A,TV B,TV C,TV D,TV E,TV F)
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
                rn.Fill_Uniform(p,-1,1);
                if(TETRAHEDRON<T>::Barycentric_Inside(p,TV(),A,B,C) && TETRAHEDRON<T>::Barycentric_Inside(p,TV(),D,E,F)) in++;
            }
    return (T)sign*in/N/N/N*8;
}

void Case_Test()
{
#if 1
    TV a,b,c,d,e,f;
    rn.Fill_Uniform(a,-1,1);
    rn.Fill_Uniform(b,-1,1);
    rn.Fill_Uniform(c,-1,1);
    rn.Fill_Uniform(d,-1,1);
    rn.Fill_Uniform(e,-1,1);
    rn.Fill_Uniform(f,-1,1);
#else
    TV a(-1,-.5,.95);
    TV b(1,-.5,.95);
    TV c(0,1,.95);
    TV d(1,.5,1);
    TV e(-1,.5,1);
    TV f(0,-1,1);
#endif


    static int cnt=0;cnt++;
    for(int i=1;i<=3;i++){
        char file[100];
        sprintf(file, "dump-%c-%i.eps", 'x'+i-1, cnt);
        EPS_FILE_GEOMETRY<T> eps(file);
        eps.Line_Color(VECTOR<T,3>(.5,.5,.5));
        eps.Draw_Point(VECTOR<T,2>(0,0));
        eps.Draw_Point(VECTOR<T,2>(1,1));
        eps.Draw_Point(VECTOR<T,2>(-1,-1));
        eps.Line_Color(VECTOR<T,3>(1,0,0));
        eps.Draw_Line(a.Remove_Index(i),b.Remove_Index(i));
        eps.Draw_Line(b.Remove_Index(i),c.Remove_Index(i));
        eps.Draw_Line(c.Remove_Index(i),a.Remove_Index(i));
        eps.Line_Color(VECTOR<T,3>(0,1,0));
        eps.Draw_Line(d.Remove_Index(i),e.Remove_Index(i));
        eps.Draw_Line(e.Remove_Index(i),f.Remove_Index(i));
        eps.Draw_Line(f.Remove_Index(i),d.Remove_Index(i));
        
        eps.Line_Color(VECTOR<T,3>(1,0,0));
        eps.Draw_Point(a.Remove_Index(i));
        eps.Draw_Point(d.Remove_Index(i));

        eps.Line_Color(VECTOR<T,3>(0,1,0));
        eps.Draw_Point(b.Remove_Index(i));
        eps.Draw_Point(e.Remove_Index(i));

        eps.Line_Color(VECTOR<T,3>(0,0,1));
        eps.Draw_Point(c.Remove_Index(i));
        eps.Draw_Point(f.Remove_Index(i));}

    TRIANGLE_ORIGIN_AREAS::VOL_DATA<T,6> data;
    TRIANGLE_ORIGIN_AREAS::Volume_From_Triangles(data,a,b,c,d,e,f);

//    if(!data.V) return;
    T approx=Approximate_Volume(a,b,c,d,e,f);

    printf("VOLUMES: %9.6f %9.6f (%.6f)\n", data.V, approx, fabs(data.V-approx));
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
        break;
    }

    return 0;
}

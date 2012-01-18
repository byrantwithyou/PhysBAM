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
int N=10000000;

T Approximate_Volume(TV A,TV B,TV C,TV D,TV E,TV F)
{
    int in=0;
    int sign=1;
    if(TV::Triple_Product(A,B,C)<0) exchange(A,B),sign=-sign;
    if(TV::Triple_Product(D,E,F)<0) exchange(D,E),sign=-sign;

    for(int z=0;z<N;z++)
    {
        TV p;
        rn.Fill_Uniform(p,-1,1);
        if(TETRAHEDRON<T>::Barycentric_Inside(p,TV(),A,B,C) && TETRAHEDRON<T>::Barycentric_Inside(p,TV(),D,E,F)) in++;
    }
    return (T)sign*in/N*8;
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

    printf("a=TV(%g,%g,%g);\n",a.x,a.y,a.z);
    printf("b=TV(%g,%g,%g);\n",b.x,b.y,b.z);
    printf("c=TV(%g,%g,%g);\n",c.x,c.y,c.z);
    printf("d=TV(%g,%g,%g);\n",d.x,d.y,d.z);
    printf("e=TV(%g,%g,%g);\n",e.x,e.y,e.z);
    printf("f=TV(%g,%g,%g);\n",f.x,f.y,f.z);

    static int cnt=0;cnt++;
    if(0)
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

    ORIGIN_AREAS::VOL_DATA<T,3,6> data;
    TV pts[6]={a,b,c,d,e,f};
    ORIGIN_AREAS::Volume_From_Simplices(data,TV(),pts);

//    if(!data.V) return;
    T approx=Approximate_Volume(a,b,c,d,e,f);

    printf("VOLUMES: %9.6f %9.6f (%.6f)\n", data.V, approx, fabs(data.V-approx));
}

void Derivative_Test()
{
    T e=1e-8;
    VECTOR<TV,6> a,da,ada;
    for(int i=1;i<=6;i++){
        rn.Fill_Uniform(a(i),-1,1);
        rn.Fill_Uniform(da(i),-e,e);
        ada(i)=a(i)+da(i);}

    for(int i=1;i<=6;i++) printf("a(%i)=TV(%g,%g,%g);\n",i,a(i).x,a(i).y,a(i).z);

    ORIGIN_AREAS::VOL_DATA<T,3,6> data0,data1;
    ORIGIN_AREAS::Volume_From_Simplices(data0,TV(),&a(1));
    ORIGIN_AREAS::Volume_From_Simplices(data1,TV(),&ada(1));

    T dV=(data1.V-data0.V)/e;
    T G=0;
    for(int i=1;i<=6;i++) G+=TV::Dot_Product(da(i),data0.G[i-1]+data1.G[i-1])/2;
    G/=e;
    printf("GRAD %9.6f %9.6f (%.6f)\n", dV, G, fabs(dV-G));

    T H0=0,H1=0,H2=0;
    for(int i=0;i<6;i++){
        TV Hx;
        for(int j=0;j<6;j++) Hx+=(data0.H[i][j]+data1.H[i][j])*da(j+1)/(2*e);
        TV G=(data1.G[i]-data0.G[i])/e;
        H0+=G.Magnitude_Squared();
        H1+=Hx.Magnitude_Squared();
        H2+=(G-Hx).Magnitude_Squared();}
    H0=sqrt(H0);
    H1=sqrt(H1);
    H2=sqrt(H2);

    printf("HESS %9.6f %9.6f (%.6f)\n", H0, H1, H2);
}

int main(int argc,char *argv[])
{
    typedef double T;
    typedef double RW;
    typedef VECTOR<T,3> TV;
    LOG::cout<<std::setprecision(16);

//    if(0)
    for(int i=0;i<20;i++)
    {
        try{Case_Test();}catch(...){}
//        break;
    }

    for(int i=0;i<100;i++) Derivative_Test();

    return 0;
}

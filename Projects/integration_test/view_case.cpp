//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <PhysBAM_Geometry/Images/EPS_FILE.h>
#include <PhysBAM_Geometry/Images/TEX_FILE.h>
#include <map>
using namespace PhysBAM;

typedef float T;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,3> TV_INT;
typedef VECTOR<T,2> V2;

V2 latex_x(1,0),latex_y(0,1),latex_z(.3,.4);
TV points[20];
TV corners[8];

V2 to2d(TV p) {return latex_x*p.x+latex_y*p.y+latex_z*(1-p.z);}

int main(int argc, char* argv[])
{
    int case_number=-1;
    std::string file="case.tex";
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-case",&case_number,"case","case number");
    parse_args.Add("-o",&file,"dir","output filename");
    parse_args.Parse();
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE(0.f),GRID<TV>(TV_INT()+1,RANGE<TV>::Unit_Box(),true),"output");
    const VECTOR<TV_INT,8>& bits=GRID<TV>::Binary_Counts(TV_INT());

    V2 mx_pt=to2d(TV(1,1,0));

    for(int v=0;v<8;v++) corners[v]=TV(v&1,v/2&1,v/4&1);

    for(int v=0;v<8;v++) points[v+12]=corners[v];
    for(int a=0,k=0;a<3;a++){
        int mask=1<<a;
        for(int v=0;v<8;v++){
            if(!(v&mask))
                points[k++]=(corners[v]+corners[v|mask])/2;}}

    for(case_number=0;case_number<256;case_number++){
        ARRAY<VECTOR<TV,TV::m> > surface;
        VECTOR<VECTOR<ARRAY<VECTOR<TV,TV::m> >,2>,2*TV::m> boundary;
        VECTOR<VECTOR<ARRAY<VECTOR<TV,TV::m> >*,2>,2*TV::m> pboundary;
        for(int f=0;f<2*TV::m;f++)
            for(int s=0;s<2;s++)
                pboundary(f)(s)=&boundary(f)(s);

        VECTOR<T,8> phis;
        for(int i=0;i<8;i++)
            phis(i)=case_number&(1<<i)?-1:1;
        MARCHING_CUBES<TV>::Get_Elements_For_Cell(surface,pboundary,phis);

        for(int i=0;i<surface.m;i++)
            Add_Debug_Object(surface(i),TV(0,1,0),TV(1,0,0));

        for(int f=0;f<2*TV::m;f++)
            for(int s=0;s<2;s++)
                for(int i=0;i<boundary(f)(s).m;i++)
                    Add_Debug_Object(boundary(f)(s)(i),TV(1/(s+1),1/(s+1),1),TV(.5,0,.5));

        char buff[100];
        sprintf(buff, "case %i", case_number);
        for(int i=0;i<8;i++)
            Add_Debug_Particle(TV(bits(i)),TV(phis(i)<0,phis(i)>0,0));

        Flush_Frame<TV>(buff);}

    for(case_number=0;case_number<256;case_number++){
        ARRAY<VECTOR<TV,TV::m+1> > interior,exterior;

        VECTOR<T,8> phis;
        for(int i=0;i<8;i++)
            phis(i)=case_number&(1<<i)?-1:1;
        MARCHING_CUBES<TV>::Get_Interior_Elements_For_Cell(&interior,&exterior,phis);

        for(int i=0;i<interior.m;i++){
            TV C(interior(i).Sum()/4);
            interior(i)-=C;
            interior(i)*=(T).9;
            interior(i)+=C;
            TETRAHEDRON<T> tet(interior(i));
            for(int j=0;j<4;j++)
                Add_Debug_Object(tet.triangle(j).X,TV(1,0,0),TV(0,0,1));}

        for(int i=0;i<exterior.m;i++){
            TV C(exterior(i).Sum()/4);
            exterior(i)-=C;
            exterior(i)*=(T).9;
            exterior(i)+=C;
            TETRAHEDRON<T> tet(exterior(i));
            for(int j=0;j<4;j++)
                Add_Debug_Object(tet.triangle(j).X,TV(0,1,0),TV(0,0,1));}

        char buff[100];
        sprintf(buff, "case %i", case_number);
        for(int i=0;i<8;i++)
            Add_Debug_Particle(TV(bits(i)),TV(phis(i)<0,phis(i)>0,0));

        Flush_Frame<TV>(buff);}
    return 0;
}

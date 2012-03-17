//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Images/EPS_FILE.h>
#include <PhysBAM_Geometry/Images/HIDDEN_SURFACE.h>
#include <PhysBAM_Geometry/Images/HIDDEN_SURFACE_PRIMITIVES.h>
#include <PhysBAM_Geometry/Images/TEX_FILE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

#include <iostream>
#include <vector>

#include <boost/foreach.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/io/wkt/wkt.hpp>

using namespace PhysBAM;

typedef double T;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,3> TV_INT;
typedef VECTOR<T,2> TV2;
typedef VECTOR<int,2> TV_INT2;



int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-o","out.tex","output filename");
    parse_args.Add_String_Argument("-i","","input filename");
    parse_args.Parse(argc,argv);
    std::string file=parse_args.Get_String_Value("-o");
    std::string infile=parse_args.Get_String_Value("-i");

    HIDDEN_SURFACE_PRIMITIVES<T> hsp;
    HIDDEN_SURFACE<T> hs(hsp);

    TRIANGULATED_SURFACE<T> ts;
    FILE_UTILITIES::Read_From_File<float>(infile,ts);
    ts.Update_Triangle_List();
    for(int i=0;i<ts.triangle_list->m;i++)
        hsp.Add((*ts.triangle_list)(i));

    LOG::cout<<std::setprecision(16);
    RANDOM_NUMBERS<T> random;
    int N=ts.triangle_list->m;
    ARRAY<TV> colors(N);
    random.Fill_Uniform(colors,0,1);

#if 0
    for(int i=0;i<N;i++){
        TV a,b,c;
        random.Fill_Uniform(a,0,1);
        random.Fill_Uniform(b,0,1);
        random.Fill_Uniform(c,0,1);
        printf("    hsp.Add(TRIANGLE_3D<T>(TV(%.16f,%.16f,%.16f),TV(%.16f,%.16f,%.16f),TV(%.16f,%.16f,%.16f)));\n",a.x,a.y,a.z,b.x,b.y,b.z,c.x,c.y,c.z);
        TRIANGLE_3D<T> tri(a,b,c);
        if(tri.Area()<1e-10) continue;
        hsp.Add(tri);}
#endif
    hs.Compute();

    VECTOR_IMAGE<T>* vi;
    if(file.length()>=4 && file.substr(file.length()-4)==".eps")
        vi=new EPS_FILE<T>(file);
    else vi=new TEX_FILE<T>(file);

//    TV2 mx_pt=TV2(1,1);
//    T margin=.05;
//    vi->Use_Fixed_Bounding_Box(RANGE<TV2>(TV2()-margin,mx_pt+margin));
    vi->cur_format.fill_style=1;
    vi->cur_format.line_style=0;

    for(int i=0;i<hsp.order.m;i++){
        SURFACE_PRIMITIVE<T>& sp=hsp.primitives(hsp.order(i));
        for(int c=0;c<sp.projection.Size();c++){
            SURFACE_PRIMITIVE<T>::POLYGON& poly=sp.projection(c);
            vi->cur_format.fill_color=colors(sp.parent);
            ARRAY<ARRAY_VIEW<TV2> > holes(poly.inners().Size());
            if(holes.m) printf("HOLES: %i\n", holes.m);
            for(int j=0;j<poly.inners().Size();j++){ARRAY_VIEW<TV2> av(poly.inners()(j));
                holes(j).Exchange(av);}
            ARRAY_VIEW<TV2> pts(poly.outer());
            vi->Draw_Object(pts,holes);}}

    delete vi;    

    return 0;
}

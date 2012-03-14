//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Images/EPS_FILE.h>
#include <PhysBAM_Geometry/Images/HIDDEN_SURFACE.h>
#include <PhysBAM_Geometry/Images/HIDDEN_SURFACE_PRIMITIVES.h>
#include <PhysBAM_Geometry/Images/TEX_FILE.h>

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
    parse_args.Parse(argc,argv);
    std::string file=parse_args.Get_String_Value("-o");

    VECTOR_IMAGE<T>* vi;
    if(file.length()>=4 && file.substr(file.length()-4)==".eps")
        vi=new EPS_FILE<T>(file);
    else vi=new TEX_FILE<T>(file);

    TV2 mx_pt=TV2(1,1);
    T margin=.2;
    vi->Use_Fixed_Bounding_Box(RANGE<TV2>(TV2()-margin,mx_pt+margin));
    vi->cur_format.fill_style=1;
    vi->cur_format.line_style=0;

    HIDDEN_SURFACE_PRIMITIVES<T> hsp;
    HIDDEN_SURFACE<T> hs(hsp);
    // hsp.Add_Element(TV(0,0,0),TV(1,0.2,0),TV(0.5,1.1,1),0);
    // hsp.Add_Element(TV(0.1,-0.1,1),TV(1,0.1,1),TV(0,1,0),1);
    hsp.Add(TV(-0.1,0.1,0),TV(1.1,0.3,-1),TV(1,0.6,2));
    hs.Compute();

    for(int i=0;i<hsp.order.m;i++){
        SURFACE_PRIMITIVE<T>& sp=hsp.primitives(hsp.order(i));
        vi->cur_format.fill_color=TV::Axis_Vector(sp.parent);
        PHYSBAM_ASSERT(sp.num_vertices==3);
        vi->Draw_Object(sp.vertices(0).Remove_Index(2),sp.vertices(1).Remove_Index(2),sp.vertices(2).Remove_Index(2));}

    delete vi;    

    return 0;
}

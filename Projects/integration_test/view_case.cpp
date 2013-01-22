//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
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

VECTOR_IMAGE<T>* vi;

V2 to2d(TV p) {return latex_x*p.x+latex_y*p.y+latex_z*(1-p.z);}

void norm(TV a,TV b,TV c,const TV& linecolor, const TV& color, T width, T opacity, T normal_length)
{
    TV u=(a+b+c)/3;
    TV v=u+normal_length*TV::Cross_Product(b-a,c-a).Normalized();
    vi->cur_format.line_width=width;
    vi->cur_format.line_color=color;
    vi->cur_format.arrow_style="->";
    vi->cur_format.line_opacity=opacity;
    vi->Draw_Object(to2d(u),to2d(v));
    vi->cur_format.arrow_style=0;
    vi->cur_format.line_opacity=1;
}

void tri(TV a,TV b,TV c,const TV& linecolor, const TV& color, T width, T opacity, T normal_length)
{
    // This criterion isn't quite right.
    T sg=TV::Cross_Product(c-a,b-a).z;
    if(sg>0) norm(a,b,c,linecolor,color,width,opacity,normal_length);
    vi->cur_format.line_width=width;
    vi->cur_format.line_color=linecolor;
    vi->cur_format.fill_color=color;
    vi->cur_format.fill_opacity=opacity;
    vi->cur_format.fill_style=1;
    vi->Draw_Object(to2d(a),to2d(b),to2d(c));
    vi->cur_format.fill_opacity=1;
    vi->cur_format.fill_style=0;
    if(sg<=0) norm(a,b,c,linecolor,color,width,opacity,normal_length);
}

void cube_edge(int a, int v, T width, int cs, T pt_width)
{
    int mask=1<<a;

    vi->cur_format.line_width=width;
    vi->cur_format.line_color=TV::Axis_Vector(a);
    vi->Draw_Object(to2d(corners[v]),to2d(corners[v|mask]));

    bool A=cs&(1<<v),B=cs&(1<<(v|mask));
    if(A!=B || cs==-1){
        vi->cur_format.line_style=0;
        vi->cur_format.fill_style=1;
        vi->cur_format.fill_color=TV();
        vi->Draw_Object(to2d((corners[v]+corners[v|mask])/2), pt_width);
        vi->cur_format.line_style=1;
        vi->cur_format.fill_style=0;}
}

int main(int argc, char* argv[])
{
    int case_number=-1;
    T corner_radius=(T).05,edge_radius=(T).04,edge_width=(T).02,tri_edge_width=(T).01;
    std::string file="case.tex";
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-case",&case_number,"case","case number");
    parse_args.Add("-corner_radius",&corner_radius,"radius","corner radius");
    parse_args.Add("-edge_radius",&edge_radius,"radius","edge cut marker radius");
    parse_args.Add("-edge_width",&edge_width,"width","cube edge widths");
    parse_args.Add("-tri_edge_width",&tri_edge_width,"width","triangle edge widths");
    parse_args.Add("-o",&file,"dir","output filename");
    parse_args.Parse();
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE(0.f),GRID<TV>(TV_INT()+1,RANGE<TV>::Unit_Box(),true),"output");
    const VECTOR<TV_INT,8>& bits=GRID<TV>::Binary_Counts(TV_INT());

    if(file.length()>=4 && file.substr(file.length()-4)==".eps")
        vi=new EPS_FILE<T>(file);
    else vi=new TEX_FILE<T>(file);

    V2 mx_pt=to2d(TV(1,1,0));
    T margin=.2;
    vi->Use_Fixed_Bounding_Box(RANGE<V2>(V2()-margin,mx_pt+margin));

    for(int v=0;v<8;v++) corners[v]=TV(v&1,v/2&1,v/4&1);

    for(int v=0;v<8;v++) points[v+12]=corners[v];
    for(int a=0,k=0;a<3;a++){
        int mask=1<<a;
        for(int v=0;v<8;v++){
            if(!(v&mask))
                points[k++]=(corners[v]+corners[v|mask])/2;}}

    (void) edge_radius;

    cube_edge(0, 0, edge_width, case_number, edge_radius);
    cube_edge(0, 2, edge_width, case_number, edge_radius);
    cube_edge(1, 0, edge_width, case_number, edge_radius);
    cube_edge(1, 1, edge_width, case_number, edge_radius);

    for(int v=0;v<4;v++){
        vi->cur_format.fill_color=TV((case_number>=0 && case_number&(1<<v)),0,0);
        vi->cur_format.line_style=0;
        vi->cur_format.fill_style=1;
        vi->Draw_Object(to2d(corners[v]),corner_radius);
        vi->cur_format.line_style=1;
        vi->cur_format.fill_style=0;}


    cube_edge(2, 0, edge_width, case_number, edge_radius);

    for(case_number=0;case_number<256;case_number++){
    std::multimap<float, PAIR<VECTOR<TV,3>,TV> > tris;

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

    const VECTOR<TV,4> mc_tri_col(TV(1,0,0), TV(0,1,0), TV(0,0,1), TV(1,0,1));
    const VECTOR<TV,4> ex_tri_col(TV(.3,0,0), TV(0,.3,0), TV(0,0,.3), TV(.3,0,.3));
    typedef std::pair<float, PAIR<VECTOR<TV,3>,TV> > pr;

    for(int i=0;i<surface.m;i++){
        Add_Debug_Object(surface(i),TV(0,1,0),TV(1,0,0));
        tris.insert(pr(surface(i).Sum().z/3,PAIR<VECTOR<TV,3>,TV>(surface(i),TV(1,0,0))));}

    for(int f=0;f<2*TV::m;f++)
        for(int s=0;s<2;s++)
            for(int i=0;i<boundary(f)(s).m;i++){
                Add_Debug_Object(boundary(f)(s)(i),TV(1/(s+1),1/(s+1),1),TV(.5,0,.5));
                tris.insert(pr(boundary(f)(s)(i).Sum().z/3,PAIR<VECTOR<TV,3>,TV>(boundary(f)(s)(i),ex_tri_col[f/2]/(s+1))));}

    for(std::map<float, PAIR<VECTOR<TV,3>,TV> >::iterator it=tris.begin(); it!=tris.end(); it++)
        tri(it->second.x(0),it->second.x(1),it->second.x(2), TV(), it->second.y, tri_edge_width, .5, .1);

    char buff[100];
    sprintf(buff, "case %i", case_number);
    for(int i=0;i<8;i++)
        Add_Debug_Particle(TV(bits(i)),TV(phis(i)<0,phis(i)>0,0));

    Flush_Frame<TV>(buff);}

    cube_edge(2, 1, edge_width, case_number, edge_radius);
    cube_edge(2, 2, edge_width, case_number, edge_radius);
    cube_edge(2, 3, edge_width, case_number, edge_radius);
    cube_edge(0, 4, edge_width, case_number, edge_radius);
    cube_edge(0, 6, edge_width, case_number, edge_radius);
    cube_edge(1, 4, edge_width, case_number, edge_radius);
    cube_edge(1, 5, edge_width, case_number, edge_radius);

    vi->cur_format.fill_style=1;
    vi->cur_format.line_style=0;
    for(int v=4;v<8;v++){
        vi->cur_format.fill_color=TV(case_number>=0 && case_number&(1<<v),0,0);
        vi->Draw_Object(to2d(corners[v]),corner_radius);}
    vi->cur_format.fill_style=0;
    vi->cur_format.line_style=1;

    delete vi;
    return 0;
}

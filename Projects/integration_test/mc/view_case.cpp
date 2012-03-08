//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <map>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
using namespace PhysBAM;

typedef float T;
typedef VECTOR<T,3> TV;
typedef VECTOR<T,2> V2;

V2 latex_x(1,0),latex_y(0,1),latex_z(.3,.4);
const char * rgb[3] = {"red", "green", "blue"};
TV points[20];
TV corners[8];

FILE * F = 0;

V2 to2d(TV p) {return latex_x*p.x+latex_y*p.y+latex_z*(1-p.z);}

void pt(V2 p) {fprintf(F, "(%.3f,%.3f)", p.x, p.y);}

void pt(TV p) {pt(to2d(p));}

void line(V2 a,V2 b,const char* opts)
{
    fprintf(F, "\\psline[%s]", opts);
    pt(a);
    pt(b);
    fprintf(F, "\n");
}

void circ(V2 a,T r,const char* opts)
{
    fprintf(F, "\\pscircle[%s]", opts);
    pt(a);
    fprintf(F, "{%.3f}\n", r);
}

void line(TV a,TV b,const char* opts) {line(to2d(a), to2d(b), opts);}

void circ(TV a,T r,const char* opts) {circ(to2d(a), r, opts);}

void norm(TV a,TV b,TV c,const char* linecolor, const char* color, T width, T opacity, T normal_length)
{
    TV u=(a+b+c)/3;
    TV v=u+normal_length*TV::Cross_Product(b-a,c-a).Normalized();
    fprintf(F, "\\psline[linewidth=%.3f,linecolor=%s]{->}", width, color);
    pt(u);
    pt(v);
    fprintf(F, "\n");
}

void tri(TV a,TV b,TV c,const char* linecolor, const char* color, T width, T opacity, T normal_length)
{
    // This criterion isn't quite right.
    T sg=TV::Cross_Product(c-a,b-a).z;
    if(sg>0) norm(a,b,c,linecolor,color,width,opacity,normal_length);
    fprintf(F, "\\pspolygon[linestyle=none,fillcolor=%s,fillstyle=solid,opacity=%.3f]", color, opacity);
    pt(a);
    pt(b);
    pt(c);
    fprintf(F, "\n");
    fprintf(F, "\\pspolygon[linecolor=%s,linewidth=%.3f]", linecolor, width);
    pt(a);
    pt(b);
    pt(c);
    fprintf(F, "\n");
    if(sg<=0) norm(a,b,c,linecolor,color,width,opacity,normal_length);
}

void cube_edge(int a, int v, T width, int cs, T pt_width)
{
    char buff[100];
    int mask=1<<a;
    sprintf(buff, "linecolor=%s,linewidth=%.3f", rgb[a], width);
    line(corners[v],corners[v|mask],buff);
    bool A=cs&(1<<v),B=cs&(1<<(v|mask));
    sprintf(buff, "fillstyle=solid,fillcolor=black,linestyle=none");
    if(A!=B || cs==-1) circ((corners[v]+corners[v|mask])/2, pt_width, buff);
}

int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args;
    parse_args.Add_Integer_Argument("-case",-1,"case number");
    parse_args.Add_Double_Argument("-corner_radius",.05,"corner radius");
    parse_args.Add_Double_Argument("-edge_radius",.04,"edge cut marker radius");
    parse_args.Add_Double_Argument("-edge_width",.02,"cube edge widths");
    parse_args.Add_Double_Argument("-tri_edge_width",.01,"triangle edge widths");
    parse_args.Add_String_Argument("-o","case.tex","output filename");
    parse_args.Parse(argc,argv);
    int case_number=parse_args.Get_Integer_Value("-case");
    T corner_radius=parse_args.Get_Double_Value("-corner_radius");
    T edge_radius=parse_args.Get_Double_Value("-edge_radius");
    T edge_width=parse_args.Get_Double_Value("-edge_width");
    T tri_edge_width=parse_args.Get_Double_Value("-tri_edge_width");
    std::string file=parse_args.Get_String_Value("-o");

    const ARRAY<MARCHING_CUBES_CASE<3> >& table=MARCHING_CUBES<TV>::Case_Table();
    const MARCHING_CUBES_CASE<3>& cs = table(case_number>=0?case_number:0);

    F = fopen(file.c_str(), "w");
    char buff[1000];

    V2 mx_pt=to2d(TV(1,1,0));
    T margin=.5;
    fprintf(F, "\\documentclass{article}\n");
    fprintf(F, "\\usepackage[margin=0cm,papersize={%.3fcm,%.3fcm}]{geometry}\n", mx_pt.x*5+5*margin, mx_pt.y*5+5*margin);
    fprintf(F, "\\usepackage{pstricks}\n");
    fprintf(F, "\\usepackage{color}\n");
    fprintf(F, "\\definecolor{dkred}{rgb}{0.3,0,0}\n");
    fprintf(F, "\\definecolor{dkgreen}{rgb}{0,0.3,0}\n");
    fprintf(F, "\\definecolor{dkblue}{rgb}{0,0,0.3}\n");
    fprintf(F, "\\definecolor{dkmagenta}{rgb}{0.3,0,0.3}\n");
    fprintf(F, "\\begin{document}\n");
    fprintf(F, "\\psset{unit=5cm}\n");
    fprintf(F, "\\noindent\\begin{pspicture}(%.3f,%.3f)(%.3f,%.3f)\n", -margin/2, -margin/2, mx_pt.x+margin/2, mx_pt.y+margin/2);

    for(int v=0;v<8;v++) corners[v]=TV(v&1,v/2&1,v/4&1);

    for(int v=0;v<8;v++) points[v+12]=corners[v];
    for(int a=0,k=0;a<3;a++){
        int mask=1<<a;
        for(int v=0;v<8;v++){
            if(!(v&mask))
                points[k++]=(corners[v]+corners[v|mask])/2;}}

    (void) cs;
    (void) edge_radius;

    cube_edge(0, 4, edge_width, case_number, edge_radius);
    cube_edge(0, 6, edge_width, case_number, edge_radius);
    cube_edge(1, 4, edge_width, case_number, edge_radius);
    cube_edge(1, 5, edge_width, case_number, edge_radius);

    for(int v=4;v<8;v++){
        sprintf(buff, "fillstyle=solid,fillcolor=%s,linestyle=none", (case_number>=0 && case_number&(1<<v))?"red":"black");
        circ(corners[v],corner_radius,buff);}

    cube_edge(2, 0, edge_width, case_number, edge_radius);
    cube_edge(2, 1, edge_width, case_number, edge_radius);
    cube_edge(2, 2, edge_width, case_number, edge_radius);

    std::multimap<float, PAIR<int,const char*> > tris;

    const char * mc_tri_col[4] = {"red", "green", "blue", "magenta"};
    const char * ex_tri_col[4] = {"dkred", "dkgreen", "dkblue", "dkmagenta"};
    typedef std::pair<float, PAIR<int,const char*> > pr;
    for(int i=0,c=-1;i<MARCHING_CUBES_CASE<3>::max_surface && cs.surface[i];i++){
        if(cs.surface[i]&0x8000) c++;
        tris.insert(pr(((points[cs.surface[i]&31]+points[(cs.surface[i]>>5)&31]+points[(cs.surface[i]>>10)&31])/3).z,PAIR<int,const char*>(cs.surface[i],mc_tri_col[c])));}

    for(int i=0,c=-1;i<MARCHING_CUBES_CASE<3>::max_boundary && cs.boundary[i];i++){
        if(cs.boundary[i]&0x8000) c++;
        tris.insert(pr(((points[cs.boundary[i]&31]+points[(cs.boundary[i]>>5)&31]+points[(cs.boundary[i]>>10)&31])/3).z,PAIR<int,const char*>(cs.boundary[i],ex_tri_col[c])));}

    for(std::map<float, PAIR<int,const char*> >::iterator it=tris.begin(); it!=tris.end(); it++)
        tri(points[it->second.x&31], points[(it->second.x>>5)&31], points[(it->second.x>>10)&31], "black", it->second.y, tri_edge_width, .5, .1);

    cube_edge(2, 3, edge_width, case_number, edge_radius);
    cube_edge(0, 0, edge_width, case_number, edge_radius);
    cube_edge(0, 2, edge_width, case_number, edge_radius);
    cube_edge(1, 0, edge_width, case_number, edge_radius);
    cube_edge(1, 1, edge_width, case_number, edge_radius);

    for(int v=0;v<4;v++){
        sprintf(buff, "fillstyle=solid,fillcolor=%s,linestyle=none", (case_number>=0 && case_number&(1<<v))?"red":"black");
        circ(corners[v],corner_radius,buff);}

    if(case_number>=0){
        TV p(1.2,-.1,1);
        T dx=.2;
        TV a(1,0,0),b(0,0,1);
        if(cs.proj_dir==0) a=TV(0,1,0);
        if(cs.proj_dir==2) b=TV(0,1,0);
        TV q=p+(TV(1,1,1)-a-b)*dx;
        fprintf(F, "\\pspolygon[linestyle=none,fillcolor=%s,fillstyle=solid,opacity=%.3f]", rgb[cs.proj_dir], .5);
        pt(q);
        pt(q+a*dx);
        pt(q+a*dx+b*dx);
        pt(q+b*dx);
        fprintf(F, "\n");
        sprintf(buff, "fillstyle=solid,fillcolor=%s,linestyle=none", cs.enclose_inside?"red":"black");
        circ(TV(1.3,0,.9),corner_radius/2,buff);
        fprintf(F, "\\pspolygon[linestyle=none,fillcolor=%s,fillstyle=solid,opacity=%.3f]", rgb[cs.proj_dir], .5);
        pt(p);
        pt(p+a*dx);
        pt(p+a*dx+b*dx);
        pt(p+b*dx);
        fprintf(F, "\n");}

    fprintf(F, "\\end{pspicture}\n");
    fprintf(F, "\\end{document}\n");
    fclose(F);

    return 0;
}

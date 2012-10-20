//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <map>
using namespace PhysBAM;

typedef float T;
typedef VECTOR<T,2> TV;

const char * rgb[3] = {"red", "green", "blue"};
TV points[8]={TV(0,0),TV(1,0),TV(0,1),TV(1,1),};
TV corners[4]={TV(0,0),TV(1,0),TV(0,1),TV(1,1)};

FILE * F = 0;

void pt(TV p) {fprintf(F, "(%.3f,%.3f)", p.x, p.y);}

void line(TV a,TV b,const char* opts)
{
    fprintf(F, "\\psline[%s]", opts);
    pt(a);
    pt(b);
    fprintf(F, "\n");
}

void circ(TV a,T r,const char* opts)
{
    fprintf(F, "\\pscircle[%s]", opts);
    pt(a);
    fprintf(F, "{%.3f}\n", r);
}

void norm(TV a,TV b,const char* linecolor, const char* color, T width, T normal_length)
{
    TV u=(a+b)/2;
    TV v=u-normal_length*(b-a).Unit_Orthogonal_Vector();
    fprintf(F, "\\psline[linewidth=%.3f,linecolor=%s]", width, color);
    pt(a);
    pt(b);
    fprintf(F, "\n");
    fprintf(F, "\\psline[linewidth=%.3f,linecolor=%s]{->}", width, color);
    pt(u);
    pt(v);
    fprintf(F, "\n");
}

int main(int argc, char* argv[])
{
    int case_number=-1;
    T corner_radius=(T).05,tri_edge_width=(T).025;
    std::string file="case.tex";
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-case",&case_number,"case","case number");
    parse_args.Add("-corner_radius",&corner_radius,"radius","corner radius");
    parse_args.Add("-tri_edge_width",&tri_edge_width,"width","triangle edge widths");
    parse_args.Add("-o",&file,"file","output filename");
    parse_args.Parse();

    for(int v=0;v<4;v++) corners[v]=TV(v&1,v/2&1);
    for(int v=0;v<4;v++) points[v+4]=corners[v];
    for(int a=0,k=0;a<2;a++){
        int mask=1<<a;
        for(int v=0;v<4;v++){
            if(!(v&mask))
                points[k++]=(corners[v]+corners[v|mask])/2;}}

    const ARRAY<MARCHING_CUBES_CASE<2> >& table=MARCHING_CUBES<TV>::Case_Table();
    const MARCHING_CUBES_CASE<2>& cs = table(case_number>=0?case_number:0);

    F = fopen(file.c_str(), "w");
    char buff[1000];

    TV mx_pt(1,1);
    T margin=.5;
    fprintf(F, "\\documentclass{article}\n");
    fprintf(F, "\\usepackage[margin=0cm,papersize={%.3fcm,%.3fcm}]{geometry}\n", mx_pt.x*5+5*margin, mx_pt.y*5+5*margin);
    fprintf(F, "\\usepackage{pstricks}\n");
    fprintf(F, "\\usepackage{color}\n");
    fprintf(F, "\\definecolor{dkred}{rgb}{0.5,0,0}\n");
    fprintf(F, "\\definecolor{dkgreen}{rgb}{0,0.5,0}\n");
    fprintf(F, "\\definecolor{dkblue}{rgb}{0,0,0.5}\n");
    fprintf(F, "\\definecolor{dkmagenta}{rgb}{0.5,0,0.5}\n");
    fprintf(F, "\\begin{document}\n");
    fprintf(F, "\\psset{unit=5cm}\n");
    fprintf(F, "\\noindent\\begin{pspicture}(%.3f,%.3f)(%.3f,%.3f)\n", -margin/2, -margin/2, mx_pt.x+margin/2, mx_pt.y+margin/2);

    fprintf(F, "\\psframe[linewidth=.02](0,0)(1,1)\n");

    const char * mc_tri_col[4] = {"red", "green", "blue", "magenta"};
    const char * ex_tri_col[4] = {"dkred", "dkgreen", "dkblue", "dkmagenta"};
    typedef std::pair<float, PAIR<int,const char*> > pr;
    for(int i=0,c=-1;i<MARCHING_CUBES_CASE<2>::max_surface && cs.surface[i];i++){
        if(cs.surface[i]&0x8000) c++;
        norm(points[cs.surface[i]&31], points[(cs.surface[i]>>5)&31], "black", mc_tri_col[c], tri_edge_width, .2);}

    for(int i=0,c=-1;i<MARCHING_CUBES_CASE<2>::max_boundary && cs.boundary[i];i++){
        if(cs.boundary[i]&0x8000) c++;
        norm(points[cs.boundary[i]&31], points[(cs.boundary[i]>>5)&31], "black", ex_tri_col[c], tri_edge_width, .2);}

    for(int v=0;v<4;v++){
        sprintf(buff, "fillstyle=solid,fillcolor=%s,linestyle=none", (case_number>=0 && case_number&(1<<v))?"red":"black");
        circ(corners[v],corner_radius,buff);}

    if(case_number>=0){
        TV p(.75,-.1);
        T dx=.15;
        p-=dx/2;
        fprintf(F, "\\psline[linecolor=%s,linewidth=.01,opacity=%.3f]", rgb[cs.proj_dir], .5);
        pt(p);
        pt(p+TV::Axis_Vector(1-cs.proj_dir)*dx);
        fprintf(F, "\n");
        sprintf(buff, "fillstyle=solid,fillcolor=%s,linestyle=none", cs.enclose_inside?"red":"black");
        circ(p+dx/2,corner_radius/2,buff);
        fprintf(F, "\\psline[linecolor=%s,linewidth=.01,opacity=%.3f]", rgb[cs.proj_dir], .5);
        pt(p+TV::Axis_Vector(cs.proj_dir)*dx);
        pt(p+TV(1,1)*dx);
        fprintf(F, "\n");}

    fprintf(F, "\\end{pspicture}\n");
    fprintf(F, "\\end{document}\n");
    fclose(F);

    return 0;
}

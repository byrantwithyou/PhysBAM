#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Data_Structures/UNION_FIND.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Images/EPS_FILE.h>
#include <Geometry/Images/TEX_FILE.h>
#include <Geometry/Topology/TRIANGLE_MESH.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include "MARCHING_TETRAHEDRA_CUTTING.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef double T;
    typedef VECTOR<int,2> TV_INT;
    typedef VECTOR<T,2> TV;

    int seed=time(0),size=4;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-s",&seed,"value","Seed");
    parse_args.Add("-n",&size,"value","Domain size");
    parse_args.Parse();

#if 0
    ARRAY<VECTOR<int,3> > m;
    m.Append(VECTOR<int,3>(0,1,2));
    ARRAY<T> p0;
    p0.Append(0.5);
    p0.Append(0);
    p0.Append(-0.1);
    ARRAY<TV> X(3);
    X(0)=TV(0,0);
    X(1)=TV(1,0);
    X(2)=TV(0,1);
    ARRAY<VECTOR<int,3> > c=m,sp;
    ARRAY<PAIR<VECTOR<int,2>,T> > weights;
    LOG::cout<<m<<"    "<<c<<std::endl;
    MARCHING_TETRAHEDRA_CUTTING<TV>::Query_Case(m,c,sp,p0,weights);
    LOG::cout<<m<<"    "<<c<<"    "<<sp<<"    "<<weights<<std::endl;
    for(int i=0;i<weights.m;i++)
        X.Append(X(weights(i).x.x)*(1-weights(i).y)+X(weights(i).x.y)*weights(i).y);

    EPS_FILE<T> eps("out.eps");
    for(int i=0;i<c.m;i++)
        eps.Draw_Object(X(c(i)(0)),X(c(i)(1)),X(c(i)(2)));
    p0.Resize(X.m);
    for(int i=0;i<p0.m;i++){
        eps.cur_format.fill_color=VECTOR<T,3>(p0(i)<0,p0(i)>0,0);
        eps.cur_format.fill_style=1;
        eps.cur_format.line_style=0;
        eps.Draw_Object(X(i),(T).02);}

    return 0;

#else

    TEX_FILE<T> eps("out.tex");
    GRID<TV> grid(TV_INT()+(size+1),RANGE<TV>::Unit_Box());

    TRIANGULATED_AREA<T>& tv=*TRIANGULATED_AREA<T>::Create();
    tv.Initialize_Herring_Bone_Mesh_And_Particles(grid);

    tv.mesh.Initialize_Boundary_Mesh();
    tv.mesh.Initialize_Boundary_Nodes();
    printf("boundary: %i\n",tv.mesh.boundary_mesh->elements.m);

    ARRAY<T> phi0(tv.particles.number),phi1(tv.particles.number);
    RANDOM_NUMBERS<T> random;
    printf("seed %i\n", seed);
    random.Set_Seed(seed);
    random.Fill_Uniform(phi0,-5,5);
    random.Fill_Uniform(phi1,-5,5);
    for(int i=0;i<phi0.m;i++) phi0(i)=::PhysBAM::rint(phi0(i));
    for(int i=0;i<phi1.m;i++) phi1(i)=::PhysBAM::rint(phi1(i));
    phi0.Subset(*tv.mesh.boundary_nodes).Fill(1);
    phi1.Subset(*tv.mesh.boundary_nodes).Fill(1);

    eps.cur_format.line_width=.01;
    eps.cur_format.line_color=VECTOR<T,3>(1,0,0);
    for(int i=0;i<tv.mesh.elements.m;i++){
        VECTOR<int,3> e=tv.mesh.elements(i);
        ARRAY<TV> X;
        for(int j=0;j<3;j++) if(phi0(e(j))==0) X.Append(tv.particles.X(e(j)));
        for(int j=0;j<3;j++){
            T p=phi0(e(j)),q=phi0(e((j+1)%3));
            TV Y=tv.particles.X(e(j)),Z=tv.particles.X(e((j+1)%3));
            if(p*q<0)
                X.Append(Y+(Z-Y)*p/(p-q));}
        if(X.m==2)
            eps.Draw_Object(X(0),X(1));}

    eps.cur_format.line_color=VECTOR<T,3>(0,1,0);
    for(int i=0;i<tv.mesh.elements.m;i++){
        VECTOR<int,3> e=tv.mesh.elements(i);
        ARRAY<TV> X;
        for(int j=0;j<3;j++) if(phi1(e(j))==0) X.Append(tv.particles.X(e(j)));
        for(int j=0;j<3;j++){
            T p=phi1(e(j)),q=phi1(e((j+1)%3));
            TV Y=tv.particles.X(e(j)),Z=tv.particles.X(e((j+1)%3));
            if(p*q<0)
                X.Append(Y+(Z-Y)*p/(p-q));}
        if(X.m==2)
            eps.Draw_Object(X(0),X(1));}

    eps.cur_format.line_color=VECTOR<T,3>();
    eps.cur_format.line_width=.002;
    ARRAY<VECTOR<int,3> > m=tv.mesh.elements,sp,tmp0,tmp1;
    TRIANGLE_MESH tm;
    ARRAY<PAIR<VECTOR<int,2>,T> > weights;
    ARRAY<bool> side;
    MARCHING_TETRAHEDRA_CUTTING<TV>::Query_Case(m,tv.mesh.elements,tmp0,tmp1,sp,side,phi0,weights);
    m=tmp0;
    tv.mesh.elements=tmp1;
    tv.particles.Add_Elements(weights.m);
    for(int i=0;i<weights.m;i++)
        tv.particles.X(i+phi0.m)=tv.particles.X(weights(i).x.x)*(1-weights(i).y)+tv.particles.X(weights(i).x.y)*weights(i).y;
    for(int i=0;i<weights.m;i++)
        phi1.Append(phi1(weights(i).x.x)*(1-weights(i).y)+phi1(weights(i).x.y)*weights(i).y);
    phi0.Resize(tv.particles.X.m);

    weights.Remove_All();
    tmp0.Remove_All();
    tmp1.Remove_All();
    side.Remove_All();
    MARCHING_TETRAHEDRA_CUTTING<TV>::Query_Case(sp,tv.mesh.elements,tmp0,tmp1,tm.elements,side,phi1,weights);
    sp=tmp0;
    tv.mesh.elements=tmp1;
    tm.Set_Number_Nodes(tm.elements.Flattened().Max()+1);
    tv.particles.Add_Elements(weights.m);
    for(int i=0;i<weights.m;i++)
        tv.particles.X(i+phi1.m)=tv.particles.X(weights(i).x.x)*(1-weights(i).y)+tv.particles.X(weights(i).x.y)*weights(i).y;
    for(int i=0;i<weights.m;i++)
        phi0.Append(phi0(weights(i).x.x)*(1-weights(i).y)+phi0(weights(i).x.y)*weights(i).y);

    UNION_FIND<> uf(tm.elements.m);
    ARRAY<VECTOR<T,3> > colors(tm.elements.m);
    random.Fill_Uniform(colors,0,1);
    tm.Initialize_Neighbor_Elements();
    for(int i=0;i<tm.elements.m;i++)
        for(int j=0;j<(*tm.neighbor_elements)(i).m;j++)
            uf.Union(i,(*tm.neighbor_elements)(i)(j));

    eps.cur_format.fill_opacity=.5;
    eps.cur_format.fill_style=1;
    for(int i=0;i<tv.mesh.elements.m;i++){
        eps.cur_format.fill_color=colors(uf.Find(i));
        eps.Draw_Object(tv.particles.X(tv.mesh.elements(i)(0)),tv.particles.X(tv.mesh.elements(i)(1)),tv.particles.X(tv.mesh.elements(i)(2)));}
    phi1.Resize(tv.particles.X.m);
    for(int i=0;i<tv.particles.X.m;i++){
        eps.cur_format.fill_color=VECTOR<T,3>(.5*(1+(phi0(i)<0)-(phi0(i)>0)),.5*(1+(phi1(i)<0)-(phi1(i)>0)),.5);
        eps.cur_format.fill_style=1;
        eps.cur_format.line_style=0;
        eps.Draw_Object(tv.particles.X(i),.01);}

    tv.Update_Number_Nodes();
    tv.mesh.Initialize_Boundary_Nodes();
    printf("boundary: %i\n",tv.mesh.boundary_mesh->elements.m);
#endif
    return 0;
}

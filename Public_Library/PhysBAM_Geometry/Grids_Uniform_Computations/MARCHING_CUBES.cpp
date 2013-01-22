//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#define TRI(a,b,c,s) (((s)<<15) | ((c)<<10) | ((b)<<5) | (a))
#define TRI_ORIENT_MAP(x,f) (x?((x&0x8000) | (f[(x>>10)&31]<<10) | (f[x&31]<<5) | f[(x>>5)&31]):0)
namespace PhysBAM{
template<> int MARCHING_CUBES_CASE<1>::edge_lookup[2][2]={{-1,0},{0,-1}};
template<> int MARCHING_CUBES_CASE<1>::vertex_lookup[1][2]={{0,1}};
template<> int MARCHING_CUBES_CASE<1>::permute_map[2][3]={};
template<> int MARCHING_CUBES_CASE<1>::face_map[2][1]={{1},{2}};
template<> int MARCHING_CUBES_CASE<2>::edge_lookup[4][4]={{-1,0,2,-1},{0,-1,-1,3},{2,-1,-1,1},{-1,3,1,-1}};
template<> int MARCHING_CUBES_CASE<2>::vertex_lookup[4][2]={{0,1},{2,3},{0,2},{1,3}};
template<> int MARCHING_CUBES_CASE<2>::permute_map[2][8]={{2,3,0,1,4,6,5,7},{0,1,3,2,5,4,7,6}};
template<> int MARCHING_CUBES_CASE<2>::face_map[4][3]={{2,4,6},{3,5,7},{0,4,5},{1,6,7}};
template<> int MARCHING_CUBES_CASE<3>::edge_lookup[8][8]={};
template<> int MARCHING_CUBES_CASE<3>::vertex_lookup[12][2]={};
template<> int MARCHING_CUBES_CASE<3>::permute_map[6][20]={};
#define MCC(d,e,f,g,A,b,c) {d,e,f,g,A+12,A+b+12,A+c+12,A+b+c+12}
template<> int MARCHING_CUBES_CASE<3>::face_map[6][8]={MCC(4,6,8,10,0,2,4),MCC(5,7,9,11,1,2,4),MCC(0,2,8,9,0,1,4),MCC(1,3,10,11,2,1,4),MCC(0,1,4,5,0,1,2),MCC(2,3,6,7,4,1,2)};
static void Initialize_Case_Table(ARRAY<MARCHING_CUBES_CASE<1> >& table);
static void Initialize_Case_Table(ARRAY<MARCHING_CUBES_CASE<2> >& table);
static void Initialize_Case_Table(ARRAY<MARCHING_CUBES_CASE<3> >& table);
template<int d> static void Initialize_Interior_Case_Table(ARRAY<MARCHING_CUBES_INTERIOR_CASE<d> >& table);
//#####################################################################
// Function Case_Table
//#####################################################################
template<int d> const ARRAY<MARCHING_CUBES_CASE<d> >& MARCHING_CUBES_CASE<d>::
Case_Table()
{
    static ARRAY<MARCHING_CUBES_CASE<d> > table;
    static bool filled=false;
    if(!filled){
        Initialize_Case_Table(table);
        filled=true;}
    return table;
}
//#####################################################################
// Function Case_Table
//#####################################################################
template<int d> const ARRAY<MARCHING_CUBES_INTERIOR_CASE<d> >& MARCHING_CUBES_INTERIOR_CASE<d>::
Case_Table()
{
    static ARRAY<MARCHING_CUBES_INTERIOR_CASE<d> > table;
    static bool filled=false;
    if(!filled){
        Initialize_Interior_Case_Table(table);
        filled=true;}
    return table;
}
//#####################################################################
// Function has_ambig
//#####################################################################
static int has_ambig(int n)
{
    if((n&0x0f)==0x09 || (n&0x0f)==0x06) return true;
    if((n&0xf0)==0x90 || (n&0xf0)==0x60) return true;
    if((n&0x33)==0x21 || (n&0x33)==0x12) return true;
    if((n&0xaa)==0x82 || (n&0xaa)==0x28) return true;
    if((n&0x55)==0x41 || (n&0x55)==0x14) return true;
    if((n&0xcc)==0x84 || (n&0xcc)==0x48) return true;
    return false;
}
//#####################################################################
// Function Initialize_Case_Table
//#####################################################################
template<int d> static void Initialize_Neighbor_Cases(ARRAY<MARCHING_CUBES_CASE<d> >& table, int c)
{
    for(int p=0;p<4*d-6;p++){
        int b=0;
        for(int i=0;i<MARCHING_CUBES_CASE<d>::num_corners;i++)
            if(c&(1<<i))
                b|=1<<(MARCHING_CUBES_CASE<d>::permute_map[p][i+MARCHING_CUBES_CASE<d>::num_edges]-MARCHING_CUBES_CASE<d>::num_edges);
        if(!table(b).elements[0]){
            for(int i=0;i<MARCHING_CUBES_CASE<d>::max_elements;i++) table(b).elements[i]=TRI_ORIENT_MAP(table(c).elements[i],MARCHING_CUBES_CASE<d>::permute_map[p]);
            Initialize_Neighbor_Cases(table, b);}}

    int b=(1<<MARCHING_CUBES_CASE<d>::num_corners)-c-1;
    if((d==2 || !has_ambig(c)) && !table(b).elements[0]){
        table(b)=table(c);
        MARCHING_CUBES_CASE<d>& cs=table(b);
        for(int i=0;i<MARCHING_CUBES_CASE<d>::max_elements;i++){
            int r=cs.elements[i]&0xFC00,s=cs.elements[i]&0x03E0,t=cs.elements[i]&0x001F;
            cs.elements[i]=r|(s>>5)|(t<<5);}
        Initialize_Neighbor_Cases(table, b);}
}
//#####################################################################
// Function Initialize_Case_Table
//#####################################################################
static void Initialize_Case_Table(ARRAY<MARCHING_CUBES_CASE<1> >& table)
{
    MARCHING_CUBES_CASE<1> d0 = {{}};
    MARCHING_CUBES_CASE<1> d1 = {{TRI(0,0,0,1)}};
    table.Resize(4);
    table(0)=d0;
    table(1)=d1;
    table(2)=d1;
    table(3)=d0;
}
static MARCHING_CUBES_CASE<2> d0 = {{}};
static MARCHING_CUBES_CASE<2> d1 = {{TRI(0,2,0,1)}};
static MARCHING_CUBES_CASE<2> d5 = {{TRI(0,1,0,1)}};
static MARCHING_CUBES_CASE<2> d9 = {{TRI(0,2,0,1), TRI(1,3,0,1)}};
//#####################################################################
// Function Initialize_Case_Table
//#####################################################################
static void Initialize_Case_Table(ARRAY<MARCHING_CUBES_CASE<2> >& table)
{
    table.Resize(16);
    table(0)=d0;
    table(1)=d1;
    table(5)=d5;
    table(9)=d9;
    table(15)=d0;

    Initialize_Neighbor_Cases(table,1);
    Initialize_Neighbor_Cases(table,5);
    Initialize_Neighbor_Cases(table,9);
}
static MARCHING_CUBES_CASE<3> c0 = {{}};
static MARCHING_CUBES_CASE<3> c1 = {{TRI(8,0,4,1)}};
static MARCHING_CUBES_CASE<3> c3 = {{TRI(9,4,8,1), TRI(4,9,5,0)}};
static MARCHING_CUBES_CASE<3> c6 = {{TRI(9,5,0,1), TRI(1,10,4,1)}};
static MARCHING_CUBES_CASE<3> c22 = {{TRI(9,5,0,1), TRI(6,2,8,1), TRI(1,10,4,1)}};
static MARCHING_CUBES_CASE<3> c23 = {{TRI(1,10,6,1), TRI(2,1,6,0), TRI(5,1,2,0), TRI(9,5,2,0)}};
static MARCHING_CUBES_CASE<3> c24 = {{TRI(6,2,8,1), TRI(1,5,11,1)}};
static MARCHING_CUBES_CASE<3> c25 = {{TRI(6,2,4,1), TRI(4,2,0,0), TRI(1,5,11,1)}};
static MARCHING_CUBES_CASE<3> c27 = {{TRI(6,2,4,1), TRI(4,2,1,0), TRI(1,2,9,0), TRI(11,1,9,0)}};
static MARCHING_CUBES_CASE<3> c60 = {{TRI(6,7,8,1), TRI(8,7,9,0), TRI(10,4,5,1), TRI(10,5,11,0)}};
static MARCHING_CUBES_CASE<3> c61 = {{TRI(6,7,9,1), TRI(0,6,9,0), TRI(5,6,0,0), TRI(11,6,5,0), TRI(10,6,11,0)}};
static MARCHING_CUBES_CASE<3> c69 = {{TRI(6,1,3,1), TRI(8,1,6,0), TRI(0,1,8,0)}};
static MARCHING_CUBES_CASE<3> c85 = {{TRI(0,1,2,1), TRI(2,1,3,0)}};
static MARCHING_CUBES_CASE<3> c101 = {{TRI(6,1,3,1), TRI(8,1,6,0), TRI(0,1,8,0), TRI(2,7,9,1)}};
static MARCHING_CUBES_CASE<3> c105 = {{TRI(4,8,0,1), TRI(7,9,2,1), TRI(11,1,5,1), TRI(10,3,6,1)}};
static MARCHING_CUBES_CASE<3> c125 = {{TRI(5,11,3,1), TRI(5,3,0,0), TRI(0,3,9,0), TRI(9,3,7,0)}};
static MARCHING_CUBES_CASE<3> c151 = {{TRI(1,10,6,1), TRI(2,1,6,0), TRI(5,1,2,0), TRI(9,5,2,0), TRI(3,11,7,1)}};
static MARCHING_CUBES_CASE<3> c255 = {{}};
//#####################################################################
// Function Initialize_Case_Table
//#####################################################################
static void Initialize_Case_Table(ARRAY<MARCHING_CUBES_CASE<3> >& table)
{
    for(int a=0,k=0;a<3;a++){
        int mask=1<<a;
        for(int v=0;v<8;v++)
            if(!(v&mask)){
                MARCHING_CUBES_CASE<3>::vertex_lookup[k][0]=v;
                MARCHING_CUBES_CASE<3>::vertex_lookup[k][1]=v|mask;
                MARCHING_CUBES_CASE<3>::edge_lookup[v][v|mask]=MARCHING_CUBES_CASE<3>::edge_lookup[v|mask][v]=k++;}}

    for(int a=0;a<3;a++) for(int v=0;v<8;v++) MARCHING_CUBES_CASE<3>::permute_map[a][v+12]=(v^(1<<a))+12;
    for(int v=0;v<8;v++) MARCHING_CUBES_CASE<3>::permute_map[3][v+12]=((v&1)|(v&2)*2|(v&4)/2)+12;
    for(int v=0;v<8;v++) MARCHING_CUBES_CASE<3>::permute_map[4][v+12]=((v&1)*4|(v&2)|(v&4)/4)+12;
    for(int v=0;v<8;v++) MARCHING_CUBES_CASE<3>::permute_map[5][v+12]=((v&1)*2|(v&2)/2|(v&4))+12;

    for(int i=0;i<6;i++)
        for(int e=0;e<12;e++){
            int v0=MARCHING_CUBES_CASE<3>::permute_map[i][MARCHING_CUBES_CASE<3>::vertex_lookup[e][0]+12]-12;
            int v1=MARCHING_CUBES_CASE<3>::permute_map[i][MARCHING_CUBES_CASE<3>::vertex_lookup[e][1]+12]-12;
            MARCHING_CUBES_CASE<3>::permute_map[i][e]=MARCHING_CUBES_CASE<3>::edge_lookup[v0][v1];}

    table.Resize(256);
    table(0)=c0;
    table(1)=c1;
    table(3)=c3;
    table(6)=c6;
    table(22)=c22;
    table(23)=c23;
    table(24)=c24;
    table(25)=c25;
    table(27)=c27;
    table(60)=c60;
    table(61)=c61;
    table(69)=c69;
    table(85)=c85;
    table(101)=c101;
    table(105)=c105;
    table(125)=c125;
    table(151)=c151;
    table(255)=c255;

    Initialize_Neighbor_Cases(table,1);
    Initialize_Neighbor_Cases(table,3);
    Initialize_Neighbor_Cases(table,6);
    Initialize_Neighbor_Cases(table,22);
    Initialize_Neighbor_Cases(table,23);
    Initialize_Neighbor_Cases(table,24);
    Initialize_Neighbor_Cases(table,25);
    Initialize_Neighbor_Cases(table,27);
    Initialize_Neighbor_Cases(table,60);
    Initialize_Neighbor_Cases(table,61);
    Initialize_Neighbor_Cases(table,69);
    Initialize_Neighbor_Cases(table,85);
    Initialize_Neighbor_Cases(table,101);
    Initialize_Neighbor_Cases(table,105);
    Initialize_Neighbor_Cases(table,125);
    Initialize_Neighbor_Cases(table,151);
}
//#####################################################################
// Function Get_Interface_Elements
//#####################################################################
template<int d> static int Get_Interface_Elements(VECTOR<int,d>* elements,int cs)
{
    const MARCHING_CUBES_CASE<d>& c=MARCHING_CUBES_CASE<d>::Case_Table()(cs);

    int count=0;
    for(int i=0;i<MARCHING_CUBES_CASE<d>::max_elements && c.elements[i];i++){
        VECTOR<int,d>& e=elements[count++];
        for(int j=0;j<d;j++)
            e(j)=(c.elements[i]>>5*j)&31;}

    return count;
}
//#####################################################################
// Function Get_Interior_Elements
//#####################################################################
template<int dp1> static int Get_Interior_Elements(VECTOR<int,dp1>* elements,bool inside,int cs)
{
    const MARCHING_CUBES_INTERIOR_CASE<dp1-1>& c=MARCHING_CUBES_INTERIOR_CASE<dp1-1>::Case_Table()(cs);

    int count=0;
    for(int i=0;i<MARCHING_CUBES_INTERIOR_CASE<dp1-1>::max_elements && c.elements[inside][i];i++){
        VECTOR<int,dp1>& e=elements[count++];
        for(int j=0;j<dp1;j++)
            e(j)=(c.elements[inside][i]>>5*j)&31;}
    return count;
}
//#####################################################################
// Function Get_Interior_Elements
//#####################################################################
static int Get_Interior_Elements(VECTOR<int,1>* elements,bool inside,int cs)
{
    if(inside!=cs) return 0;
    elements[0].x=0;
    return 1;
}
//#####################################################################
// Function Get_Boundary_Elements
//#####################################################################
template<int d> static int Get_Boundary_Elements(VECTOR<int,d>* elements,int axis,int side,bool inside,int cs)
{
    int cs2=0;
    for(int i=0,k=0;i<(1<<d);i++)
        if(((i>>axis)&1)==side)
            cs2|=((cs>>i)&1)<<k++;

    int count=Get_Interior_Elements(elements,inside,cs2);
    const int* face_map=MARCHING_CUBES_CASE<d>::face_map[2*axis+side];

    for(int i=0;i<count;i++){
        for(int j=0;j<d;j++)
            elements[i](j)=face_map[elements[i](j)];
        int flip=d>1 && side==(axis==1);
        exchange(elements[i](0),elements[i](flip));}

    return count;
}
//#####################################################################
// Function Get_Boundary_Elements
//#####################################################################
template<int d> static int Get_Boundary_Elements(VECTOR<int,d>* elements,bool inside,int cs)
{
    int count=0;
    for(int a=0;a<d;a++)
        for(int s=0;s<2;s++)
            count+=Get_Boundary_Elements(elements+count,a,s,inside,cs);
    return count;
}
//#####################################################################
// Function Initialize_Interior_Case_Table
//#####################################################################
template<int d> static void Initialize_Interior_Case_Table(ARRAY<MARCHING_CUBES_INTERIOR_CASE<d> >& table,int cs)
{
    MARCHING_CUBES_INTERIOR_CASE<d> ec={};
    table(cs)=ec;
    MARCHING_CUBES_INTERIOR_CASE<d>& cas=table(cs);
    UNION_FIND<> uf(MARCHING_CUBES_CASE<d>::num_pts);
    for(int in=0;in<2;in++){
        int j=0;
        VECTOR<int,d> elements[2*d*(d-1)+MARCHING_CUBES_CASE<d>::max_elements];
        int count=Get_Interface_Elements(elements,cs);
        if(!in)
            for(int i=0;i<count;i++)
                exchange(elements[i](0),elements[i](1));
        count+=Get_Boundary_Elements(elements+count,in,cs);

        VECTOR<int,MARCHING_CUBES_CASE<d>::num_pts> smallest;
        uf.Clear_Connectivity();
        for(int i=0;i<count;i++) uf.Union(elements[i]);
        for(int i=MARCHING_CUBES_CASE<d>::num_corners-1;i>=0;i--)
            if(((cs>>i)&1)==in)
                smallest(uf.Find(i+MARCHING_CUBES_CASE<d>::num_edges))=i+MARCHING_CUBES_CASE<d>::num_edges;

        for(int i=0;i<count;i++){
            int k=smallest(uf.Find(elements[i][0]));
            if(!elements[i].Contains(k)){
                VECTOR<int,d+1> el=elements[i].Append(k);
                int elem=0;
                for(int t=0;t<d+1;t++)
                    elem|=el(t)<<5*t;
                cas.elements[in][j++]=elem;}}}
}
//#####################################################################
// Function Initialize_Interior_Case_Table
//#####################################################################
template<int d> static void Initialize_Interior_Case_Table(ARRAY<MARCHING_CUBES_INTERIOR_CASE<d> >& table)
{
    int cases=1<<(1<<d);
    table.Resize(cases);
    for(int i=0;i<cases;i++)
        Initialize_Interior_Case_Table(table,i);
}
//#####################################################################
// Function Initialize_Interior_Case_Table
//#####################################################################
static void Initialize_Interior_Case_Table(ARRAY<MARCHING_CUBES_INTERIOR_CASE<1> >& table)
{
    MARCHING_CUBES_INTERIOR_CASE<1> c0={{{TRI(1,2,0,1)},{}}};
    MARCHING_CUBES_INTERIOR_CASE<1> c1={{{TRI(0,2,0,1)},{TRI(1,0,0,1)}}};
    MARCHING_CUBES_INTERIOR_CASE<1> c2={{{TRI(1,0,0,1)},{TRI(0,2,0,1)}}};
    MARCHING_CUBES_INTERIOR_CASE<1> c3={{{},{TRI(1,2,0,1)}}};
    table.Resize(4);
    table(0)=c0;
    table(1)=c1;
    table(2)=c2;
    table(3)=c3;
}
//#####################################################################
// Function Get_Triangles_For_Cell
//#####################################################################
template<class TV> void MARCHING_CUBES<TV>::
Compute_Phis_For_Cell(VECTOR<T,num_corners>& phis,const ARRAY<T,TV_INT>& phi,const TV_INT& cell)
{
    const VECTOR<TV_INT,num_corners>& bits=GRID<TV>::Binary_Counts(TV_INT());
    for(int i=0;i<bits.m;i++)
        phis(i)=phi(cell+bits(i));
}
//#####################################################################
// Function Get_Triangles_For_Cell
//#####################################################################
template<class TV> int MARCHING_CUBES<TV>::
Compute_Points_For_Cell(VECTOR<TV,num_pts>& pts,const VECTOR<T,num_corners>& phis)
{
    MARCHING_CUBES_CASE<TV::m>::Case_Table();

    const VECTOR<TV_INT,num_corners>& bits=GRID<TV>::Binary_Counts(TV_INT());
    int c=0;
    for(int i=0;i<bits.m;i++){
        pts[i+num_edges]=TV(bits(i));
        c|=(phis(i)<0)<<i;}

    for(int i=0;i<num_edges;i++){
        int v0=MARCHING_CUBES_CASE<TV::m>::vertex_lookup[i][0];
        int v1=MARCHING_CUBES_CASE<TV::m>::vertex_lookup[i][1];
        if(((c>>v0)^(c>>v1))&1){
            T t=phis(v0)/(phis(v0)-phis(v1));
            pts[i]=pts[v0+num_edges]+t*(pts[v1+num_edges]-pts[v0+num_edges]);}}
    return c;
}
//#####################################################################
// Function Get_Elements_For_Cell
//#####################################################################
template<class TV> void MARCHING_CUBES<TV>::
Get_Elements_For_Cell(ARRAY<VECTOR<TV,TV::m> >& surface,const VECTOR<VECTOR<ARRAY<VECTOR<TV,TV::m> >*,2>,2*TV::m>& boundary,const ARRAY<T,TV_INT>& phi,const TV_INT& cell)
{
    VECTOR<T,num_corners> phis;
    Compute_Phis_For_Cell(phis,phi,cell);
    Get_Elements_For_Cell(surface,boundary,phis);
}
//#####################################################################
// Function Get_Elements_For_Cell
//#####################################################################
template<class TV> void MARCHING_CUBES<TV>::
Get_Elements_For_Cell(ARRAY<VECTOR<TV,TV::m> >& surface,const VECTOR<VECTOR<ARRAY<VECTOR<TV,TV::m> >*,2>,2*TV::m>& boundary,const VECTOR<T,num_corners>& phis)
{
    VECTOR<TV,num_pts> pts;
    int c=Compute_Points_For_Cell(pts,phis);

    VECTOR<int,TV::m> tmp_elements[MARCHING_CUBES_INTERIOR_CASE<TV::m>::max_elements];
    int len=Get_Interface_Elements(tmp_elements,c);
    for(int i=0;i<len;i++) surface.Append(VECTOR<TV,TV::m>(pts.Subset(tmp_elements[i])));

    for(int a=0;a<TV::m;a++)
        for(int s=0;s<2;s++)
            for(int t=0;t<2;t++)
                if(ARRAY<VECTOR<TV,TV::m> >* ar=boundary(2*a+s)(t)){
                    len=Get_Boundary_Elements(tmp_elements,a,s,t,c);
                    for(int i=0;i<len;i++) ar->Append(VECTOR<TV,TV::m>(pts.Subset(tmp_elements[i])));}
}
//#####################################################################
// Function Create_Surface
//#####################################################################
template<class TV> int MARCHING_CUBES<TV>::
Create_Surface(T_SURFACE& surface,const GRID<TV>& grid,const ARRAY<T,TV_INT>& phi)
{
    const VECTOR<TV_INT,num_corners>& bits=GRID<TV>::Binary_Counts(TV_INT());
    HASHTABLE<FACE_INDEX<TV::m>,int> ht;

    int cut_cells=0;
    for(RANGE_ITERATOR<TV::m> it(phi.domain.To_Closed());it.Valid();it.Next()){
        VECTOR<T,num_corners> phis;
        VECTOR<TV,num_pts> pts;
        Compute_Phis_For_Cell(phis,phi,it.index);
        int c=Compute_Points_For_Cell(pts,phis);

        VECTOR<int,TV::m> tmp_elements[MARCHING_CUBES_CASE<TV::m>::max_elements];
        int len=Get_Interface_Elements(tmp_elements,c);
        if(len) cut_cells++;
        for(int i=0;i<len;i++){
            TV_INT face;
            for(int j=0;j<TV::m;j++){
                int e=tmp_elements[i](j);
                FACE_INDEX<TV::m> fi(e/(1<<(TV::m-1)),it.index+bits(MARCHING_CUBES_CASE<TV::m>::vertex_lookup[e][0]));
                if(!ht.Get(fi,face(j))){
                    int index=surface.particles.Add_Element();
                    face(j)=index;
                    ht.Set(fi,index);
                    surface.particles.X(index)=pts[e];}}
            surface.mesh.elements.Append(face);}}

    surface.Update_Number_Nodes();
    return cut_cells;
}
//#####################################################################
// Function Create_Surface
//#####################################################################
template<class TV,class TV_INT,class T> static int Create_Surface_Helper(POINT_SIMPLICES_1D<T>& surface,const GRID<TV>& grid,const ARRAY<T,TV_INT>& phi)
{
    int cut_cells=0;
    for(RANGE_ITERATOR<TV::m> it(phi.domain.To_Closed());it.Valid();it.Next()){
        T a=phi(it.index),b=phi(it.index+1);
        if((a>0)==(b>0)) continue;
        cut_cells++;
        TV X0=grid.Center(it.index),X1=grid.Center(it.index+1);
        surface.mesh.elements.Append(TV_INT(surface.particles.Add_Element()));
        surface.particles.X.Last()=X0+a/(a-b)*(X1-X0);}
    return cut_cells;
}
template class MARCHING_CUBES<VECTOR<float,1> >;
template class MARCHING_CUBES<VECTOR<float,2> >;
template class MARCHING_CUBES<VECTOR<float,3> >;
template class MARCHING_CUBES<VECTOR<double,1> >;
template class MARCHING_CUBES<VECTOR<double,2> >;
template class MARCHING_CUBES<VECTOR<double,3> >;
}

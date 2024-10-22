//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/PAIR.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Math_Tools/clamp.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Utilities/DEBUG_CAST.h>
#include <Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Basic_Geometry/LINE_2D.h>
#include <Geometry/Basic_Geometry/PLANE.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES_COLOR.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

using namespace PhysBAM;
namespace{
ARRAY<int> interface_case_table;
ARRAY<int> interface_triangle_table;
ARRAY<int> boundary_case_table;
ARRAY<int> boundary_triangle_table;
ARRAY<int> interface_case_table_2d;
ARRAY<int> interface_segment_table;
const int comparison_bit=1<<31;
const int last_tri_bit=1<<30;
#define GET_C(t,i) (((t)>>3*(i))&0x7)
#define GET_V(t,i) (((t)>>(5*(i)+6))&0x1f)
}
void Pr_Tri(int t,const char* str="")
{
    printf("tri %s %i %i  %i %i %i %i %c\n", str, GET_C(t,0), GET_C(t,1), GET_V(t,0), GET_V(t,1), GET_V(t,2), GET_V(t,3), t&last_tri_bit?'L':' ');
}
namespace{
int Encode_Triangle(int c0,int c1,int v0,int v1,int v2,int v3=31)
{
    return (v3<<21)|(v2<<16)|(v1<<11)|(v0<<6)|(c1<<3)|c0;
}
void Add_Triangle(int c0,int c1,int v0,int v1,int v2,int v3=31)
{
    interface_triangle_table.Append(Encode_Triangle(c0,c1,v0,v1,v2,v3));
}
void Add_Mapped_Triangle(int t,const int* mp,const int* cmp,bool flip)
{
    int v3=GET_V(t,3);
    Add_Triangle(cmp[GET_C(t,0)],cmp[GET_C(t,1)],mp[GET_V(t,0)],mp[GET_V(t,1+flip)],mp[GET_V(t,2-flip)],v3==31?31:mp[v3]);
    if(t&last_tri_bit)
        interface_triangle_table.Last()|=last_tri_bit;
}
/* Triangle encoding
 * 0 e 0000 ddddd ccccc bbbbb aaaaa xxx yyy   (triangle)
 * 1 0 aaa bbb ccc ddd ssssss ssssss ssssss (comparison)
 *
 * indexing: [0-11 = edges] [12-17; axis*2+side] [18 center]
 * faces: [0-3 = corners] [edges: 1 0 3 4] [8 center]
 * 
 */
// start, end, inside color, outside color (c0<c1)
struct EDGE {int v0,v1,c0,c1;void Flip(){exchange(c0,c1);exchange(v0,v1);}};
EDGE center_edges[6][4]; // c0,c1 are vertex to take color from
struct DOMINANT_PAIR {EDGE a,b;int dc;} dominant_pair[6][2];
const int permute_rx[19]={1,3,0,2,10,11,8,9,4,5,6,7,12,13,16,17,15,14,18};
const int permute_ry[19]={8,10,9,11,6,4,7,5,2,0,3,1,17,16,14,15,12,13,18};
const int permute_flip[19]={0,1,2,3,5,4,7,6,9,8,11,10,13,12,14,15,16,17,18};
const int permute_rx_corners[8]={2,3,6,7,0,1,4,5};
const int permute_ry_corners[8]={4,0,6,2,5,1,7,3};
const int permute_flip_corners[8]={1,0,3,2,5,4,7,6};
int face_edges[6][4];

inline EDGE Rotate_X(const EDGE& ep)
{
    EDGE e={permute_rx[ep.v0],permute_rx[ep.v1],permute_rx_corners[ep.c0],permute_rx_corners[ep.c1]};
    return e;
}
inline EDGE Rotate_Y(const EDGE& ep)
{
    EDGE e={permute_ry[ep.v0],permute_ry[ep.v1],permute_ry_corners[ep.c0],permute_ry_corners[ep.c1]};
    return e;
}
inline DOMINANT_PAIR Rotate_X(const DOMINANT_PAIR& ep)
{
    DOMINANT_PAIR e={Rotate_X(ep.a),Rotate_X(ep.b),permute_rx_corners[ep.dc]};
    return e;
}
inline DOMINANT_PAIR Rotate_Y(const DOMINANT_PAIR& ep)
{
    DOMINANT_PAIR e={Rotate_Y(ep.a),Rotate_Y(ep.b),permute_ry_corners[ep.dc]};
    return e;
}
void Add_Center_Edges(int* colors,EDGE* edges,int& num_edges,int axis,int s)
{
    for(int i=0;i<4;i++){
        const EDGE& e=center_edges[2*axis+s][i];
        EDGE ec={e.v0,e.v1,colors[e.c0],colors[e.c1]};
        if(ec.c0==ec.c1) continue;
        if(ec.c0>ec.c1) ec.Flip();
        edges[num_edges++]=ec;}
}
void Add_Dominat_Pair(int* colors,EDGE* edges,int& num_edges,int axis,int s,int dc)
{
    for(int i=0;i<2;i++){
        const DOMINANT_PAIR& d=dominant_pair[2*axis+s][i];
        if(colors[d.dc]!=dc) continue;
        EDGE a={d.a.v0,d.a.v1,colors[d.a.c0],colors[d.a.c1]},b={d.b.v0,d.b.v1,colors[d.b.c0],colors[d.b.c1]};
        if(a.c0>a.c1) a.Flip();
        if(b.c0>b.c1) b.Flip();
        edges[num_edges++]=a;
        edges[num_edges++]=b;}
}
void Add_Face_Edges(int* colors,EDGE* edges,int& num_edges,int a,int b,int s,int color_hint)
{
    int mask_a=1<<a,mask_b=1<<b,axis=3-a-b,mask_s=s<<axis;
    int col[4]={colors[mask_s],colors[mask_s|mask_a],colors[mask_s|mask_b],colors[mask_s|mask_a|mask_b]};
    int have=(1<<col[0])|(1<<col[1])|(1<<col[2])|(1<<col[3]);
    if(!(have&=have-1)) return;
    if(!(have&=have-1)){
        if(col[0]==col[3] && col[1]==col[2]){
            if(color_hint&(1<<(col[0]+8*col[1]))) Add_Dominat_Pair(colors,edges,num_edges,axis,s,col[1]);
            else Add_Dominat_Pair(colors,edges,num_edges,axis,s,col[0]);
            return;}
        Add_Center_Edges(colors,edges,num_edges,axis,s);
        num_edges--;
        if(edges[num_edges-1].v0==edges[num_edges].v1) edges[num_edges-1].v0=edges[num_edges].v0;
        else edges[num_edges-1].v1=edges[num_edges].v1;
        return;}
    if(!(have&=have-1)){
        if(col[0]==col[3]) Add_Dominat_Pair(colors,edges,num_edges,axis,s,col[0]);
        else if(col[1]==col[2]) Add_Dominat_Pair(colors,edges,num_edges,axis,s,col[1]);
        else Add_Center_Edges(colors,edges,num_edges,axis,s);
        return;}
    Add_Center_Edges(colors,edges,num_edges,axis,s);
}
bool on_face[19][19]={};
bool Try_Add_Triangle(int (*adj)[2],EDGE* edges,int (*cur_edges)[2],int& num_cur_edges)
{
    bool reduced=true;
    for(int i=0;i<12;i++){
        if(adj[i][0]==-1) continue;
        EDGE& in=edges[adj[i][0]],&out=edges[adj[i][1]];
        PHYSBAM_ASSERT(in.c0==out.c0 && in.c1==out.c1);
        if(in.v0==out.v1) continue;
        if(on_face[in.v0][out.v1] && out.v1<12 && edges[adj[out.v1][1]].v1!=in.v0){
            reduced=false;
            continue;}

        Add_Triangle(in.c0,in.c1,in.v0,in.v1,out.v1);
        in.v1=out.v1;
        int adj_in_v1=-1;
        if(in.v1<12){
            adj_in_v1=adj[in.v1][0];
            adj[in.v1][0]=adj[i][0];}
        int a0=adj[i][0],a1=adj[i][1];
        adj[i][0]=-1;
        adj[i][1]=-1;
        cur_edges[num_cur_edges][0]=in.v0;
        cur_edges[num_cur_edges++][1]=out.v1;

        if(Try_Add_Triangle(adj,edges,cur_edges,num_cur_edges))
            return true;
        reduced=false;

        interface_triangle_table.Pop();
        num_cur_edges--;
        adj[i][0]=a0;
        adj[i][1]=a1;
        if(in.v1<12) adj[in.v1][0]=adj_in_v1;
        in.v1=i;}

    return reduced;
}
//#####################################################################
// Function Emit_Interface_Triangles
//#####################################################################
void Emit_Interface_Triangles(int* colors,int color_hint)
{  
    EDGE edges[25]; // need an extra one
    int num_edges=0;
    for(int a=0;a<3;a++)
        for(int b=a+1;b<3;b++)
            for(int s=0;s<2;s++)
                Add_Face_Edges(colors,edges,num_edges,a,b,s,color_hint);

    int adj[12][2];
    for(int i=0;i<12;i++) for(int j=0;j<2;j++) adj[i][j]=-1;
    for(int i=0;i<num_edges;i++){
        if(edges[i].v0<12) adj[edges[i].v0][1]=i;
        if(edges[i].v1<12) adj[edges[i].v1][0]=i;}

    int cur_edges[100][2],num_cur_edges=0;
    bool loops_only=true;
    for(int i=0;i<12;i++){
        if(adj[i][0]==-1) continue;
        EDGE& in=edges[adj[i][0]];
        if(in.v0<12) continue;
        loops_only=false;
        break;}
    if(!loops_only){
        for(int i=0;i<num_edges;i++)
            Add_Triangle(edges[i].c0,edges[i].c1,edges[i].v0,edges[i].v1,18);
        interface_triangle_table.Last()|=last_tri_bit;
        return;}

    Try_Add_Triangle(adj,edges,cur_edges,num_cur_edges);
    interface_triangle_table.Last()|=last_tri_bit;
}
//#####################################################################
// Function Emit_Interface_Triangles
//#####################################################################
void Emit_Interface_Triangles(int* colors)
{
    int amb[4]={-1,-1,-1,-1};
    for(int a=0;a<3;a++)
        for(int b=a+1;b<3;b++)
            for(int s=0;s<2;s++){
                int mask_a=1<<a,mask_b=1<<b,mask_s=s<<(3-a-b);
                int c0=colors[mask_s],c1=colors[mask_s|mask_a],c2=colors[mask_s|mask_b],c3=colors[mask_s|mask_a|mask_b];
                if(c0==c3 && c1==c2 && c0!=c1){
                    int b=(amb[0]!=-1)*2;
                    amb[b]=c0;
                    amb[b+1]=c1;}}
    if(amb[0]==-1) return Emit_Interface_Triangles(colors,0);
    int test_index=interface_triangle_table.Append(0);
    if(amb[2]==-1){
        int hint=amb[0]+8*amb[1];
        Emit_Interface_Triangles(colors,1<<hint); // amb[0] < amb[1]
        int skip=interface_triangle_table.m-test_index;
        Emit_Interface_Triangles(colors,1<<(amb[1]+8*amb[0])); // amb[1] < amb[0]
        interface_triangle_table(test_index)=comparison_bit|(hint<<24)|(skip<<12);
        return;}
    int hint0=amb[0]+8*amb[1],hint1=amb[2]+8*amb[3],hintn0=amb[1]+8*amb[0],hintn1=amb[3]+8*amb[2];
    Emit_Interface_Triangles(colors,(1<<hint0)|(1<<hint1)); // amb[0] < amb[1]; amb[2] < amb[3]
    int skip0=interface_triangle_table.m-test_index;
    Emit_Interface_Triangles(colors,(1<<hint0)|(1<<hintn1)); // amb[0] < amb[1]; amb[3] < amb[2]
    int skip1=interface_triangle_table.m-test_index;
    Emit_Interface_Triangles(colors,(1<<hintn0)|(1<<hint1)); // amb[1] < amb[0]; amb[2] < amb[3]
    int skip2=interface_triangle_table.m-test_index;
    Emit_Interface_Triangles(colors,(1<<hintn0)|(1<<hintn1)); // amb[1] < amb[0]; amb[3] < amb[2]
    interface_triangle_table(test_index)=comparison_bit|(hint0<<24)|(hint1<<18)|(skip0<<12)|(skip1<<6)|skip2;
}
//#####################################################################
// Function Register_Permutation
//#####################################################################
void Register_Permutation(const int* colors,const int* perm,const int* cperm,int st,int end,bool flip)
{
    int col[8],cs=0,next_color=0,color_map[8]={-1,-1,-1,-1,-1,-1,-1,-1};
    for(int i=0;i<8;i++) col[cperm[i]]=colors[i];
    for(int i=0;i<8;i++){
        int& cm=color_map[col[i]];
        if(cm==-1) cm=next_color++;
        col[i]=cm;
        cs=cs*(i+1)+cm;}
    if(interface_case_table(cs)) return;
    int start=interface_triangle_table.m;;
    interface_case_table(cs)=start;

    int i=st,ct=interface_triangle_table(i);
    if(ct&comparison_bit){
        int t=(ct&0xc003ffff)|(color_map[(ct>>24)&7]<<24)|(color_map[(ct>>27)&7]<<27);
        if(ct&(0x3f<<18))
            t|=(color_map[(ct>>18)&7]<<18)|(color_map[(ct>>21)&7]<<21);
        interface_triangle_table.Append(t);
        i++;}
    for(;i<end;i++)
        Add_Mapped_Triangle(interface_triangle_table(i),perm,color_map,flip);

    Register_Permutation(col,permute_rx,permute_rx_corners,start,end-st+start,false);
    Register_Permutation(col,permute_ry,permute_ry_corners,start,end-st+start,false);
}
//#####################################################################
// Function Enumerate_Interface_Cases_3D
//#####################################################################
void Enumerate_Interface_Cases_3D(int* colors, int i, int mc, int cs)
{
    if(i==8){
        if(!cs || interface_case_table(cs)) return;
        int start=interface_triangle_table.m;
        interface_case_table(cs)=start;
        Emit_Interface_Triangles(colors);
        int end=interface_triangle_table.m;
        Register_Permutation(colors,permute_flip,permute_flip_corners,start,end,true);
        Register_Permutation(colors,permute_rx,permute_rx_corners,start,end,false);
        Register_Permutation(colors,permute_ry,permute_ry_corners,start,end,false);
        return;}

    for(int c=0;c<=mc;c++){
        colors[i]=c;
        Enumerate_Interface_Cases_3D(colors,i+1,mc,cs*(i+1)+c);}

    colors[i]=mc+1;
    Enumerate_Interface_Cases_3D(colors,i+1,mc+1,cs*(i+1)+(mc+1));
}
//#####################################################################
// Function Enumerate_Interface_Cases_3D
//#####################################################################
void Enumerate_Boundary_Cases_3D()
{
#define ST(c) boundary_case_table(c)=boundary_triangle_table.m
#define EMIT(v0,v1,v2,c) boundary_triangle_table.Append(Encode_Triangle(c,0,v0,v1,v2))
#define LA boundary_triangle_table.Last()|=last_tri_bit
    ST(0);EMIT(0,2,3,0);EMIT(0,3,1,0);LA;
    ST(1);EMIT(0,2,1,0);EMIT(2,7,1,0);EMIT(2,5,7,0);EMIT(5,3,7,1);LA;
    ST(4);EMIT(0,3,1,0);EMIT(0,5,3,0);EMIT(0,6,5,0);EMIT(2,5,6,1);LA;
    ST(5);EMIT(1,0,6,0);EMIT(1,6,7,0);EMIT(6,2,7,1);EMIT(2,3,7,1);LA;
    ST(6);EMIT(0,6,8,0);EMIT(8,7,1,0);EMIT(0,8,1,0);EMIT(6,2,5,1);EMIT(6,5,8,1);EMIT(5,7,8,2);EMIT(5,3,7,2);LA;
    ST(12);EMIT(0,2,3,0);EMIT(7,0,3,0);EMIT(4,0,7,0);EMIT(1,4,7,1);LA;
    ST(13);EMIT(0,2,5,0);EMIT(0,5,4,0);EMIT(5,1,4,1);EMIT(5,3,1,1);LA;
    ST(14);EMIT(4,0,8,0);EMIT(0,2,8,0);EMIT(2,5,8,0);EMIT(8,1,4,1);EMIT(8,7,1,1);EMIT(5,7,8,2);EMIT(5,3,7,2);LA;
    boundary_case_table(16)=boundary_triangle_table.m;
    int test_index=boundary_triangle_table.Append(0);
    int hint=0+8*1;
    EMIT(4,0,6,0);EMIT(2,4,6,1);EMIT(1,4,2,1);EMIT(7,1,2,1);EMIT(5,7,2,1);EMIT(5,3,7,0);LA;
    int skip=boundary_triangle_table.m-test_index;
    EMIT(4,7,1,1);EMIT(4,0,7,0);EMIT(0,3,7,0);EMIT(0,5,3,0);EMIT(0,6,5,0);EMIT(6,2,5,1);LA;
    boundary_triangle_table(test_index)=comparison_bit|(hint<<24)|(skip<<12);
    ST(17);EMIT(1,2,3,1);EMIT(1,6,2,1);EMIT(1,4,6,1);EMIT(0,6,4,0);LA;
    ST(18);EMIT(4,0,6,0);EMIT(2,4,6,1);EMIT(1,4,2,1);EMIT(7,1,2,1);EMIT(5,7,2,1);EMIT(5,3,7,2);LA;
    ST(20);EMIT(4,7,1,1);EMIT(4,0,7,0);EMIT(0,3,7,0);EMIT(0,5,3,0);EMIT(0,6,5,0);EMIT(6,2,5,2);LA;
    ST(21);EMIT(0,6,8,0);EMIT(0,8,4,0);EMIT(6,2,5,2);EMIT(6,5,8,2);EMIT(1,4,8,1);EMIT(3,1,8,1);EMIT(5,3,8,1);LA;
    ST(22);EMIT(0,6,8,0);EMIT(0,8,4,0);EMIT(8,1,4,1);EMIT(8,7,1,1);EMIT(3,7,8,2);EMIT(2,8,6,2);EMIT(2,3,8,2);LA;
    ST(23);EMIT(0,6,8,0);EMIT(0,8,4,0);EMIT(6,2,5,2);EMIT(6,5,8,2);EMIT(8,1,4,1);EMIT(8,7,1,1);EMIT(5,7,8,3);EMIT(5,3,7,3);LA;
#undef ST
#undef EMIT
#undef LA
}
//#####################################################################
// Function Initialize_Case_Table_3D
//#####################################################################
void Initialize_Case_Table_3D()
{
    EDGE e={7,13,5,7};
    for(int i=0;i<4;i++){
        center_edges[1][i]=e;
        center_edges[4][i]=Rotate_Y(center_edges[1][i]);
        center_edges[0][i]=Rotate_Y(center_edges[4][i]);
        center_edges[5][i]=Rotate_Y(center_edges[0][i]);
        center_edges[2][i]=Rotate_X(center_edges[5][i]);
        center_edges[3][i]=Rotate_X(center_edges[4][i]);
        e=Rotate_X(e);}

    DOMINANT_PAIR d={{7,9,5,1},{5,11,3,7},1};
    for(int i=0;i<2;i++){
        dominant_pair[1][i]=d;
        dominant_pair[4][i]=Rotate_Y(dominant_pair[1][i]);
        dominant_pair[0][i]=Rotate_Y(dominant_pair[4][i]);
        dominant_pair[5][i]=Rotate_Y(dominant_pair[0][i]);
        dominant_pair[2][i]=Rotate_X(dominant_pair[5][i]);
        dominant_pair[3][i]=Rotate_X(dominant_pair[4][i]);
        d=Rotate_X(d);}

    int face_count[6]={0};
    for(int a=0;a<3;a++)
        for(int b=a+1;b<3;b++)
            for(int c=0;c<2;c++)
                for(int d=0;d<2;d++){
                    int axis=3-a-b,edge=4*axis+2*d+c,f0=2*a+c,f1=2*b+d;
                    face_edges[f0][face_count[f0]++]=edge;
                    face_edges[f1][face_count[f1]++]=edge;}

    interface_case_table.Resize(40320);
    boundary_case_table.Resize(24);
    int colors[8]={0};
    Enumerate_Interface_Cases_3D(colors, 1, 0, 0);
    Enumerate_Boundary_Cases_3D();
}
//#####################################################################
// Function Initialize_Case_Table_2D
//#####################################################################
void Initialize_Case_Table_2D()
{
    interface_case_table_2d.Resize(24);
#define ST(c) interface_case_table_2d(c)=interface_segment_table.m
#define EMIT(v0,v1,c0,c1) interface_segment_table.Append(Encode_Triangle(c0,c1,v0,v1,31));
#define LA interface_segment_table.Last()|=last_tri_bit
    ST(1);EMIT(1,3,0,1);LA;
    ST(4);EMIT(2,1,0,1);LA;
    ST(5);EMIT(2,3,0,1);LA;
    ST(6);EMIT(2,4,0,1);EMIT(4,3,0,2);EMIT(1,4,1,2);LA;
    ST(12);EMIT(3,0,0,1);LA;
    ST(13);EMIT(1,0,0,1);LA;
    ST(14);EMIT(4,0,0,1);EMIT(1,4,0,2);EMIT(4,3,1,2);LA;
    interface_case_table_2d(16)=interface_segment_table.m;
    int test_index=interface_segment_table.Append(0);
    int hint=0+4*1;
    EMIT(0,2,1,0);EMIT(1,3,1,0);LA;
    int skip=interface_segment_table.m-test_index;
    EMIT(0,3,1,0);EMIT(2,1,0,1);LA;
    interface_segment_table(test_index)=comparison_bit|(hint<<24)|(skip<<12);
    ST(17);EMIT(0,2,1,0);LA;
    ST(18);EMIT(0,2,1,0);EMIT(1,3,1,2);LA;
    ST(20);EMIT(0,3,1,0);EMIT(2,1,0,2);LA;
    ST(21);EMIT(2,4,0,2);EMIT(4,0,0,1);EMIT(1,4,2,1);LA;
    ST(22);EMIT(2,4,0,2);EMIT(4,0,0,1);EMIT(4,3,1,2);LA;
    ST(23);EMIT(2,4,0,2);EMIT(4,0,0,1);EMIT(1,4,2,3);EMIT(4,3,1,3);LA;
#undef ST
#undef EMIT
#undef LA
}
const int averaging_order[7][2]={{4,6},{5,7},{0,2},{1,3},{0,1},{2,3},{12,13}};
//#####################################################################
// Function Get_Interface_Elements_For_Cell
//#####################################################################
template<class T,class INTERFACE_ELEMENT> void
Get_Interface_Elements_For_Cell(ARRAY<INTERFACE_ELEMENT>& interface,const VECTOR<int,8>& re_color,
    const VECTOR<int,8>& colors,const VECTOR<T,8>& phi,const int* color_list)
{
    typedef VECTOR<T,3> TV;
    int cs=0;
    for(int i=0;i<8;i++)
        cs=cs*(i+1)+re_color(i);
    if(!cs) return;

    int tri=interface_case_table(cs);
    if(interface_triangle_table(tri)&comparison_bit){
        int pat=interface_triangle_table(tri);
        if(pat&(0x3f<<18)){
            int amb0=(pat>>24)&7,amb1=(pat>>27)&7,amb2=(pat>>18)&7,amb3=(pat>>21)&7;
            if(color_list[amb0]<color_list[amb1]){
                if(color_list[amb2]<color_list[amb3]) tri++;
                else tri+=(pat>>12)&0x3f;}
            else{
                if(color_list[amb2]<color_list[amb3]) tri+=(pat>>6)&0x3f;
                else tri+=pat&0x3f;}}
        else{
            int amb0=(pat>>24)&7,amb1=(pat>>27)&7;
            if(color_list[amb1]<color_list[amb0]) tri+=(pat>>12)&0x3f;
            else tri++;}}

    const VECTOR<VECTOR<int,3>,8>& bits=GRID<TV>::Binary_Counts(VECTOR<int,3>());
    TV pts[19];
    T pts_phi[19];
    for(int a=0,k=0;a<3;a++){
        int mask=1<<a;
        for(int v=0;v<8;v++)
            if(!(v&mask)){
                T theta=phi(v)/(phi(v)+phi(v|mask));
                pts[k]=(1-theta)*TV(bits(v))+theta*TV(bits(v|mask));
                pts_phi[k++]=(1-theta)*phi(v)+theta*phi(v|mask);}}

    for(int i=0;i<7;i++){
        int a=averaging_order[i][0],b=averaging_order[i][1];
        T theta=pts_phi[a]/(pts_phi[a]+pts_phi[b]);
        pts[12+i]=(1-theta)*pts[a]+theta*pts[b];
        pts_phi[12+i]=(1-theta)*pts_phi[a]+theta*pts_phi[b];}

    int pat;
    do{
        pat=interface_triangle_table(tri++);
        TRIANGLE_3D<T> triangle(pts[GET_V(pat,0)],pts[GET_V(pat,2)],pts[GET_V(pat,1)]);
        INTERFACE_ELEMENT ie={triangle,VECTOR<int,2>(color_list[GET_C(pat,0)],color_list[GET_C(pat,1)])};
        if(ie.color_pair.x>ie.color_pair.y){
            exchange(ie.color_pair.x,ie.color_pair.y);
            exchange(ie.face.X(1),ie.face.X(2));}
        interface.Append(ie);
    } while(!(pat&last_tri_bit));
}
//#####################################################################
// Function Get_Boundary_Elements_For_Cell
//#####################################################################
template<class T> void
Get_Boundary_Elements_For_Cell(ARRAY<typename MARCHING_CUBES_COLOR<VECTOR<T,3> >::BOUNDARY_ELEMENT>& boundary,const int* re_color,
    const int* colors,const T* phi,int s,const int* color_list)
{
    typedef VECTOR<T,3> TV;
    typedef typename MARCHING_CUBES_COLOR<VECTOR<T,3> >::BOUNDARY_ELEMENT BOUNDARY_ELEMENT;
    int cs=0;
    for(int i=0;i<4;i++)
        cs=cs*(i+1)+re_color[i];

    int tri=boundary_case_table(cs);
    if(boundary_triangle_table(tri)&comparison_bit){
        int pat=boundary_triangle_table(tri);
        int amb0=(pat>>24)&7,amb1=(pat>>27)&7;
        if(color_list[amb1]<color_list[amb0])
            tri+=(pat>>12)&0x3f;
        else tri++;}

    const int num_pts=9;
    TV pts[num_pts]={TV(0,0,s),TV(1,0,s),TV(0,1,s),TV(1,1,s)};
    if(colors[0]!=colors[1]){T theta=phi[0]/(phi[0]+phi[1]);pts[4]=(1-theta)*pts[0]+theta*pts[1];}
    if(colors[2]!=colors[3]){T theta=phi[2]/(phi[2]+phi[3]);pts[5]=(1-theta)*pts[2]+theta*pts[3];}
    if(colors[0]!=colors[2]){T theta=phi[0]/(phi[0]+phi[2]);pts[6]=(1-theta)*pts[0]+theta*pts[2];}
    if(colors[1]!=colors[3]){T theta=phi[1]/(phi[1]+phi[3]);pts[7]=(1-theta)*pts[1]+theta*pts[3];}
    T total=0;
    for(int v=0;v<4;v++){
        total+=1/phi[v];
        pts[8]+=pts[v]/phi[v];}
    pts[8]/=total;

    int pat;
    do{
        pat=boundary_triangle_table(tri++);
        TRIANGLE_3D<T> triangle(pts[GET_V(pat,0)],pts[GET_V(pat,1)],pts[GET_V(pat,2)]);
        BOUNDARY_ELEMENT be={triangle,color_list[GET_C(pat,0)]};
        if(s) exchange(be.face.X.y,be.face.X.z);
        boundary.Append(be);
    } while(!(pat&last_tri_bit));
}
//#####################################################################
// Function Get_Interface_Elements_For_Cell
//#####################################################################
template<class T,class INTERFACE_ELEMENT> void
Get_Interface_Elements_For_Cell(ARRAY<INTERFACE_ELEMENT>& interface,const VECTOR<int,4>& re_color,
    const VECTOR<int,4>& colors,const VECTOR<T,4>& phi,const int* color_list)
{
    typedef VECTOR<T,2> TV;
    int cs=0;
    for(int i=0;i<4;i++)
        cs=cs*(i+1)+re_color(i);
    if(!cs) return;

    int seg=interface_case_table_2d(cs);
    if(interface_segment_table(seg)&comparison_bit){
        int pat=interface_segment_table(seg);
        int amb0=(pat>>24)&7,amb1=(pat>>27)&7;
        if(color_list[amb1]<color_list[amb0]) seg+=(pat>>12)&0x3f;
        else seg++;}

    TV corners[4]={TV(0,0),TV(1,0),TV(0,1),TV(1,1)};
    TV pts[5];
    int num=0;
    if(colors[0]!=colors[1]){T theta=phi[0]/(phi[0]+phi[1]);pts[0]=(1-theta)*corners[0]+theta*corners[1];num++;}
    if(colors[2]!=colors[3]){T theta=phi[2]/(phi[2]+phi[3]);pts[1]=(1-theta)*corners[2]+theta*corners[3];num++;}
    if(colors[0]!=colors[2]){T theta=phi[0]/(phi[0]+phi[2]);pts[2]=(1-theta)*corners[0]+theta*corners[2];num++;}
    if(colors[1]!=colors[3]){T theta=phi[1]/(phi[1]+phi[3]);pts[3]=(1-theta)*corners[1]+theta*corners[3];num++;}
    if(num>=3){
        T total=0;
        for(int v=0;v<4;v++){
            total+=1/phi[v];
            pts[4]+=corners[v]/phi[v];}
        pts[4]/=total;}

    int pat;
    do{
        pat=interface_segment_table(seg++);
        SEGMENT_2D<T> segment(pts[GET_V(pat,0)],pts[GET_V(pat,1)]);
        INTERFACE_ELEMENT ie={segment,VECTOR<int,2>(color_list[GET_C(pat,0)],color_list[GET_C(pat,1)])};
        if(ie.color_pair.x>ie.color_pair.y){
            exchange(ie.color_pair.x,ie.color_pair.y);
            exchange(ie.face.X(0),ie.face.X(1));}
        interface.Append(ie);
    } while(!(pat&last_tri_bit));
}
//#####################################################################
// Function Get_Boundary_Elements_For_Cell
//#####################################################################
template<class T> void
Get_Boundary_Elements_For_Cell(ARRAY<typename MARCHING_CUBES_COLOR<VECTOR<T,2> >::BOUNDARY_ELEMENT>& boundary,const int* re_color,
    const int* colors,const T* phi,int s,const int* color_list)
{
    typedef VECTOR<T,2> TV;
    typedef typename MARCHING_CUBES_COLOR<VECTOR<T,2> >::BOUNDARY_ELEMENT BOUNDARY_ELEMENT;
    if(colors[0]==colors[1]){
        BOUNDARY_ELEMENT be={SEGMENT_2D<T>(TV(s,s),TV(1-s,s)),colors[0]};
        boundary.Append(be);
        return;}

    T theta=phi[0]/(phi[0]+phi[1]);
    BOUNDARY_ELEMENT be0={SEGMENT_2D<T>(TV(theta,s),TV(1-s,s)),colors[1-s]};
    BOUNDARY_ELEMENT be1={SEGMENT_2D<T>(TV(s,s),TV(theta,s)),colors[s]};
    boundary.Append(be0);
    boundary.Append(be1);
}
}
//#####################################################################
// Function Initialize_Case_Table
//#####################################################################
template<class TV> void MARCHING_CUBES_COLOR<TV>::
Initialize_Case_Table()
{
    static bool first=true;
    if(!first) return;
    first=false;
    if(TV::m==3) Initialize_Case_Table_3D();
    if(TV::m==2) Initialize_Case_Table_2D();
}
//#####################################################################
// Function Get_Elements_For_Cell
//#####################################################################
template<class TV> void MARCHING_CUBES_COLOR<TV>::
Get_Elements_For_Cell(ARRAY<INTERFACE_ELEMENT>& interface,ARRAY<BOUNDARY_ELEMENT>& boundary,
    const VECTOR<int,num_corners>& colors,const VECTOR<T,num_corners>& phi)
{
    int next_color=0;
    HASHTABLE<int,int> color_map;
    VECTOR<int,num_corners> re_color;
    int color_list[8];
    for(int i=0;i<num_corners;i++)
        if(!color_map.Get(colors(i),re_color(i))){
            re_color(i)=next_color;
            color_list[next_color]=colors(i);
            color_map.Set(colors(i),next_color++);}

    Get_Interface_Elements_For_Cell(interface,re_color,colors,phi,color_list);
    Get_Boundary_Elements_For_Cell(boundary,re_color.array,colors.array,phi.array,0,color_list);

    next_color=0;
    color_map.Remove_All();
    int re_color2[num_corners/2];
    for(int i=0;i<num_corners/2;i++)
        if(!color_map.Get(colors(i+num_corners/2),re_color2[i])){
            re_color2[i]=next_color;
            color_list[next_color]=colors(i+num_corners/2);
            color_map.Set(colors(i+num_corners/2),next_color++);}
    Get_Boundary_Elements_For_Cell(boundary,re_color2,colors.array+num_corners/2,phi.array+num_corners/2,1,color_list);
}
int vertex_lookup_2d[4][2]={{0,1},{2,3},{0,2},{1,3}};
template<int n>
int Rename_Colors(const VECTOR<int,n> &phi_color,VECTOR<int,n>& re_color,VECTOR<int,n>& color_list,bool& junction)
{
    int next_color=0;
    HASHTABLE<int,int> color_map; // maps original colors to renamed ones
    for(int i=0;i<n;i++)
        if(!color_map.Get(phi_color(i),re_color(i))){
            re_color(i)=next_color;
            color_list[next_color]=phi_color(i);
            color_map.Set(phi_color(i),next_color++);}
    junction=(next_color>=3);

    int cs=0;
    for(int i=0;i<n;i++)
        cs=cs*(i+1)+re_color(i);
    return cs;
}
//#####################################################################
// Function Get_Hashed_Interface_Elements_For_Cell
//#####################################################################
template<class T,class TV_INT,class TV,class T_SURFACE,class HASH_INTERFACE,class HASH_CELL_INTERFACE> void
Get_Hashed_Interface_Elements_For_Cell(const TV_INT& cell_index,const VECTOR<TV,4>& corners,const VECTOR<int,4>& phi_color,const VECTOR<T,4>& phi_value,
    HASHTABLE<FACE_INDEX<2>,int>& edge_vertices,HASHTABLE<FACE_INDEX<2>,int>& face_vertices,HASHTABLE<TV_INT,int>& cell_vertices,
    HASH_INTERFACE& interface,GEOMETRY_PARTICLES<TV>& particles,HASH_CELL_INTERFACE& interface_cell_elements,bool& junction)
{
    VECTOR<int,4> re_color;    // renamed colors of cell corners
    VECTOR<int,4> color_list;  // maps renamed colors to original ones

    const int cs=Rename_Colors(phi_color,re_color,color_list,junction);
    PHYSBAM_ASSERT(cs);
    
    int seg=interface_case_table_2d(cs);
    if(interface_segment_table(seg)&comparison_bit){
        int pat=interface_segment_table(seg);
        int amb0=(pat>>24)&7,amb1=(pat>>27)&7;
        if(color_list[amb1]<color_list[amb0]) seg+=(pat>>12)&0x3f;
        else seg++;}

    const VECTOR<VECTOR<int,2>,4>& bits=GRID<TV>::Binary_Counts(TV_INT());
    VECTOR<int,5> p_index;
    p_index.Fill(-1);
    int num=0;
    for(int e=0;e<4;e++){
        int v0=vertex_lookup_2d[e][0],v1=vertex_lookup_2d[e][1];
        if(phi_color(v0)!=phi_color(v1)){
            FACE_INDEX<TV::m> fi(e/2,cell_index+bits(v0));
            num++;
            if(!edge_vertices.Get(fi,p_index[e])){
                int index=particles.Add_Element();
                p_index[e]=index;
                edge_vertices.Set(fi,index);
                T theta=phi_value(v0)/(phi_value(v0)+phi_value(v1));
                particles.X(index)=(1-theta)*corners[v0]+theta*corners[v1];}}}
    if(num>=3){
        if(!cell_vertices.Get(cell_index,p_index[4])){
            int index=particles.Add_Element();
            p_index[4]=index;
            cell_vertices.Set(cell_index,index);
            T total=0;
            TV X;
            for(int v=0;v<4;v++){
                total+=1/phi_value[v];
                X+=corners[v]/phi_value[v];}
            particles.X(index)=X/total;}}

    int pat;
    do{
        pat=interface_segment_table(seg++);
        VECTOR<int,2> c(color_list[GET_C(pat,0)],color_list[GET_C(pat,1)]);
        VECTOR<int,2> v(p_index[GET_V(pat,1)],p_index[GET_V(pat,0)]);
        if(c.x>c.y){
            exchange(c.x,c.y);
            exchange(v.x,v.y);}

        T_SURFACE* s=0;
        if(!interface.Get(c,s)){
            s=T_SURFACE::Create(particles);
            interface.Set(c,s);}
        s->mesh.elements.Append(v);
        interface_cell_elements.Get_Or_Insert(c,INTERVAL<int>(s->mesh.elements.m-1,s->mesh.elements.m)).max_corner=s->mesh.elements.m;
    } while(!(pat&last_tri_bit));
}
//#####################################################################
// Function Get_Hashed_Interface_Elements_For_Cell
//#####################################################################
template<class T,class TV_INT,class TV,class T_SURFACE,class HASH_INTERFACE,class HASH_CELL_INTERFACE> void
Get_Hashed_Interface_Elements_For_Cell(const TV_INT& cell_index,const VECTOR<TV,8>& corners,const VECTOR<int,8>& phi_color,const VECTOR<T,8>& phi_value,
    HASHTABLE<FACE_INDEX<3>,int>& edge_vertices,HASHTABLE<FACE_INDEX<3>,int>& face_vertices,HASHTABLE<TV_INT,int>& cell_vertices,
    HASH_INTERFACE& interface,GEOMETRY_PARTICLES<TV>& particles,HASH_CELL_INTERFACE& interface_cell_elements,bool& junction)
{
    VECTOR<int,8> re_color;    // renamed colors of cell corners
    VECTOR<int,8> color_list;  // maps renamed colors to original ones

    const int cs=Rename_Colors(phi_color,re_color,color_list,junction);
    PHYSBAM_ASSERT(cs);
    
    int tri=interface_case_table(cs);
    if(interface_triangle_table(tri)&comparison_bit){
        int pat=interface_triangle_table(tri);
        if(pat&(0x3f<<18)){
            int amb0=(pat>>24)&7,amb1=(pat>>27)&7,amb2=(pat>>18)&7,amb3=(pat>>21)&7;
            if(color_list[amb0]<color_list[amb1]){
                if(color_list[amb2]<color_list[amb3]) tri++;
                else tri+=(pat>>12)&0x3f;}
            else{
                if(color_list[amb2]<color_list[amb3]) tri+=(pat>>6)&0x3f;
                else tri+=pat&0x3f;}}
        else{
            int amb0=(pat>>24)&7,amb1=(pat>>27)&7;
            if(color_list[amb1]<color_list[amb0]) tri+=(pat>>12)&0x3f;
            else tri++;}}

    VECTOR<TV,19> pts;
    VECTOR<T,19> pts_phi;
    for(int a=0,k=0;a<3;a++){
        int mask=1<<a;
        for(int v=0;v<8;v++)
            if(!(v&mask)){
                T theta=phi_value(v)/(phi_value(v)+phi_value(v|mask));
                pts[k]=(1-theta)*corners[v]+theta*corners[v|mask];
                pts_phi[k++]=(1-theta)*phi_value(v)+theta*phi_value(v|mask);}}

    for(int i=0;i<7;i++){
        int a=averaging_order[i][0],b=averaging_order[i][1];
        T theta=pts_phi[a]/(pts_phi[a]+pts_phi[b]);
        pts[12+i]=(1-theta)*pts[a]+theta*pts[b];
        pts_phi[12+i]=(1-theta)*pts_phi[a]+theta*pts_phi[b];}

    int pat;
    VECTOR<int,19> p_index;
    p_index.Fill(-1);
    do{
        pat=interface_triangle_table(tri++);
        const TV_INT triangle(GET_V(pat,0),GET_V(pat,1),GET_V(pat,2));
        for(int p=0;p<TV::m;p++){
            const int e=triangle(p);
            if(e<12){
                const int axis=e>>2;
                const int offset1=e&1;
                const int offset2=(e>>1)&1;
                FACE_INDEX<TV::m> fi(axis,cell_index+VECTOR<int,2>(offset1,offset2).Insert(0,axis));
                if(!edge_vertices.Get(fi,p_index[e])){
                    int index=particles.Add_Element();
                    p_index[e]=index;
                    edge_vertices.Set(fi,index);
                    particles.X(index)=pts[e];}}
            else if(e<18){
                const int axis=(e-12)>>1;
                const int sign=(e-12)&1;
                FACE_INDEX<TV::m> fi(axis,cell_index+TV_INT::Axis_Vector(axis)*sign);
                if(!face_vertices.Get(fi,p_index[e])){
                    int index=particles.Add_Element();
                    p_index[e]=index;
                    face_vertices.Set(fi,index);
                    particles.X(index)=pts[e];}}
            else{
                if(!cell_vertices.Get(cell_index,p_index[e])){
                    int index=particles.Add_Element();
                    p_index[e]=index;
                    cell_vertices.Set(cell_index,index);
                    particles.X(index)=pts[e];}}}

        TV_INT v(p_index[triangle(0)],p_index[triangle(1)],p_index[triangle(2)]);
        VECTOR<int,2> c(color_list[GET_C(pat,0)],color_list[GET_C(pat,1)]);
        if(c.x>c.y){exchange(c.x,c.y);exchange(v(1),v(2));}
        
        T_SURFACE* s=0;
        if(!interface.Get(c,s)){
            s=T_SURFACE::Create(particles);
            interface.Set(c,s);}
        s->mesh.elements.Append(v);
        interface_cell_elements.Get_Or_Insert(c,INTERVAL<int>(s->mesh.elements.m-1,s->mesh.elements.m)).max_corner=s->mesh.elements.m;
    }while(!(pat&last_tri_bit));
}
//#####################################################################
// Function Get_Hashed_Boundary_Elements_For_Cell
//#####################################################################
template<class T,class TV_INT,class TV,class T_SURFACE,class HASH_BOUNDARY,class HASH_CELL_BOUNDARY> void
Get_Hashed_Boundary_Elements_For_Cell(const TV_INT& cell_index,const VECTOR<TV,4>& corners,const VECTOR<int,4>& phi_color,const VECTOR<T,4>& phi_value,
    HASHTABLE<FACE_INDEX<2>,int>& edge_vertices,HASHTABLE<FACE_INDEX<2>,int>& face_vertices,HASHTABLE<TV_INT,int>& node_vertices,
    HASH_BOUNDARY& boundary,GEOMETRY_PARTICLES<TV>& particles,VECTOR<HASH_CELL_BOUNDARY,2>& boundary_elements_vector)
{
    for(int a=0;a<2;a++){
        const int b=1-a;
        const int oa=1<<a;
        const int ob=1<<b;
        HASH_CELL_BOUNDARY& boundary_elements=boundary_elements_vector(b);
        for(int sb=0;sb<2;sb++){
            
            // INITIALIZE CORNER POINTS
            
            VECTOR<TV,2> pts;
            VECTOR<int,2> p_index;
            p_index.Fill(-1);

            for(int sa=0;sa<2;sa++)
                pts[sa]=corners[sa*oa+sb*ob];

            for(int sa=0;sa<2;sa++){
                TV_INT node_index=cell_index+VECTOR<int,1>(sa).Insert(sb,b);
                if(!node_vertices.Get(node_index,p_index[sa])){
                    int index=particles.Add_Element();
                    p_index[sa]=index;
                    node_vertices.Set(node_index,index);
                    particles.X(index)=pts[sa];}
                else{
                    assert(p_index[sa]==p_index[sa] || p_index[sa]==-1);
                    p_index[sa]=p_index[sa];}}

            // PROCESS CUT AND UNCUT SEGMENT CASES

            if(phi_color[sb*ob]==phi_color[oa+sb*ob]){
                VECTOR<int,2> v(p_index[sb],p_index[1-sb]);
                if(a) exchange(v.x,v.y);
                int color=phi_color[sb*ob];
                T_SURFACE* surf=0;
                if(!boundary.Get(color,surf)){
                    surf=T_SURFACE::Create(particles);
                    boundary.Set(color,surf);}
                surf->mesh.elements.Append(v);
                boundary_elements.Get_Or_Insert(color,INTERVAL<int>(surf->mesh.elements.m-1,surf->mesh.elements.m)).max_corner=surf->mesh.elements.m;}
            else{
                int p_index_m=-1;
                T theta=phi_value[sb*ob]/(phi_value[sb*ob]+phi_value[oa+sb*ob]);
                TV pts_m=pts[1]*theta+pts[0]*(1-theta);
                FACE_INDEX<TV::m> fi(a,cell_index+TV_INT::Axis_Vector(b)*sb);
                if(!edge_vertices.Get(fi,p_index_m)){
                    int index=particles.Add_Element();
                    p_index_m=index;
                    edge_vertices.Set(fi,index);
                    particles.X(index)=pts_m;}
                
                for(int t=0;t<2;t++){
                    VECTOR<int,2> v(p_index[t^sb],p_index_m);
                    if(a^t) exchange(v.x,v.y);
                    int color=phi_color[sb*ob+(t^sb)*oa];
                    
                    T_SURFACE* surf=0;
                    if(!boundary.Get(color,surf)){
                        surf=T_SURFACE::Create(particles);
                        boundary.Set(color,surf);}
                    surf->mesh.elements.Append(v);
                    boundary_elements.Get_Or_Insert(color,INTERVAL<int>(surf->mesh.elements.m-1,surf->mesh.elements.m)).max_corner=surf->mesh.elements.m;}}}}
}
//#####################################################################
// Function Get_Hashed_Boundary_Elements_For_Cell
//#####################################################################
template<class T,class TV_INT,class TV,class T_SURFACE,class HASH_BOUNDARY,class HASH_CELL_BOUNDARY> void
Get_Hashed_Boundary_Elements_For_Cell(const TV_INT& cell_index,const VECTOR<TV,8>& corners,const VECTOR<int,8>& phi_color,const VECTOR<T,8>& phi_value,
    HASHTABLE<FACE_INDEX<3>,int>& edge_vertices,HASHTABLE<FACE_INDEX<3>,int>& face_vertices,HASHTABLE<TV_INT,int>& node_vertices,
    HASH_BOUNDARY& boundary,GEOMETRY_PARTICLES<TV>& particles,VECTOR<HASH_CELL_BOUNDARY,3>& boundary_elements_vector)
{
    for(int a=0;a<3;a++)
    for(int b=a+1;b<3;b++){
        const int c=3-a-b;
        const int oa=1<<a;
        const int ob=1<<b;
        const int oc=1<<c;
        HASH_CELL_BOUNDARY& boundary_elements=boundary_elements_vector(c);
        for(int sc=0;sc<2;sc++){
            
            // RENAME COLORS

            int count=-1;
            int next_color=0;
            VECTOR<int,4> re_color;
            VECTOR<int,8> color_list;
            HASHTABLE<int,int> color_map;

            for(int sb=0;sb<2;sb++)
            for(int sa=0;sa<2;sa++){
                count++;
                const int index=sa*oa+sb*ob+sc*oc;
                if(!color_map.Get(phi_color(index),re_color[count])){
                    re_color[count]=next_color;
                    color_list[next_color]=phi_color(index);
                    color_map.Set(phi_color(index),next_color++);}}

            // RETRIEVE CASE FROM THE MC TABLE
            
            int cs=0;
            for(int i=0;i<4;i++)
                cs=cs*(i+1)+re_color[i];
            int tri=boundary_case_table(cs);
            if(boundary_triangle_table(tri)&comparison_bit){
                int pat=boundary_triangle_table(tri);
                int amb0=(pat>>24)&7,amb1=(pat>>27)&7;
                if(color_list[amb1]<color_list[amb0])
                    tri+=(pat>>12)&0x3f;
                else tri++;}

            // CALCULATE POINTS POSITIONS

            VECTOR<TV,9> pts;
            VECTOR<int,9> p_index;
            p_index.Fill(-1);
            for(int sb=0;sb<2;sb++)
            for(int sa=0;sa<2;sa++)
                pts[sa+(sb<<1)]=corners[sa*oa+sb*ob+sc*oc];
            
            static const int side_pair[4][2][2]={{{0,0},{0,1}},{{1,0},{1,1}},{{0,0},{1,0}},{{0,1},{1,1}}};

            for(int i=0;i<4;i++){
                int sa[2]={side_pair[i][0][1],side_pair[i][1][1]};
                int sb[2]={side_pair[i][0][0],side_pair[i][1][0]};
                if(phi_color[sa[0]*oa+sb[0]*ob+sc*oc]!=phi_color[sa[1]*oa+sb[1]*ob+sc*oc]){
                    T theta=phi_value[sa[0]*oa+sb[0]*ob+sc*oc]/(phi_value[sa[0]*oa+sb[0]*ob+sc*oc]+phi_value[sa[1]*oa+sb[1]*ob+sc*oc]);
                    pts[4+i]=(1-theta)*pts[sa[0]+(sb[0]<<1)]+theta*pts[sa[1]+(sb[1]<<1)];}}

            T total=0;
            for(int sb=0;sb<2;sb++)
            for(int sa=0;sa<2;sa++){
                const int index=sa*oa+sb*ob+sc*oc;
                total+=1/phi_value[index];
                pts[8]+=pts[sa+(sb<<1)]/phi_value[index];}
            pts[8]/=total;

            // FILL IN PARTICLE HASH TABLES
            
            int pat;
            do{
                pat=boundary_triangle_table(tri++);
                TV_INT triangle(GET_V(pat,0),GET_V(pat,1),GET_V(pat,2));
                for(int p=0;p<TV::m;p++){
                    const int e=triangle(p);
                    if(e<4){
                        const int sa=e&1;
                        const int sb=(e>>1)&1;
                        TV_INT node_index=cell_index+VECTOR<int,2>(sa,sb).Insert(sc,c);
                        if(!node_vertices.Get(node_index,p_index[e])){
                            int index=particles.Add_Element();
                            p_index[e]=index;
                            node_vertices.Set(node_index,index);
                            particles.X(index)=pts[e];}}
                    else if(e<8){
                        const int axis=((e>>1)&1)?b:a;
                        FACE_INDEX<TV::m> fi(axis,cell_index+VECTOR<int,2>(e==7,e==5).Insert(sc,c));
                        if(!edge_vertices.Get(fi,p_index[e])){
                            int index=particles.Add_Element();
                            p_index[e]=index;
                            edge_vertices.Set(fi,index);
                            particles.X(index)=pts[e];}}
                    else{
                        FACE_INDEX<TV::m> fi(c,cell_index+TV_INT::Axis_Vector(c)*sc);
                        if(!face_vertices.Get(fi,p_index[e])){
                            int index=particles.Add_Element();
                            p_index[e]=index;
                            face_vertices.Set(fi,index);
                            particles.X(index)=pts[e];}}}
                
                TV_INT v(p_index[triangle(0)],p_index[triangle(1)],p_index[triangle(2)]);
                int color=color_list[GET_C(pat,0)];
                if(sc^(b-a==2)) exchange(v.y,v.z);
                
                T_SURFACE* surf=0;
                if(!boundary.Get(color,surf)){
                    surf=T_SURFACE::Create(particles);
                    boundary.Set(color,surf);}
                surf->mesh.elements.Append(v);
                boundary_elements.Get_Or_Insert(color,INTERVAL<int>(surf->mesh.elements.m-1,surf->mesh.elements.m)).max_corner=surf->mesh.elements.m;
            }while(!(pat&last_tri_bit));}}
}
//#####################################################################
// Function Get_Elements
//#####################################################################
template<class TV> void MARCHING_CUBES_COLOR<TV>::
Get_Elements(HASHTABLE<TV_INT,CELL_ELEMENTS>& index_to_cell_elements,const GRID<TV>& grid,
    const ARRAY<int,TV_INT>& phi_color,const ARRAY<T,TV_INT>& phi_value,
    const int iterations,const bool verbose)
{
    // CELL CONSTANTS

    const VECTOR<VECTOR<int,TV::m>,(1<<TV::m)>& bits=GRID<TV>::Binary_Counts(TV_INT());

    // HASHTABLES FOR DIFFERENTLY CONSTRAINED PARTICLES

    HASHTABLE<FACE_INDEX<TV::m>,int> edge_vertices;
    HASHTABLE<FACE_INDEX<TV::m>,int> face_vertices;
    HASHTABLE<TV_INT,int> cell_vertices;
    HASHTABLE<TV_INT,int> node_vertices;

    // TEMPORARY STRUCTURES FOR THE MESH

    HASH_BOUNDARY boundary;
    HASH_INTERFACE interface;
    HASHTABLE<TV_INT,HASH_CELL_DATA> index_to_cell_data;
    GEOMETRY_PARTICLES<TV>& particles=*new GEOMETRY_PARTICLES<TV>;
    HASHTABLE<TV_INT> junction_cells;

    // INITIALIZE MESHES

    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){

        VECTOR<T,num_corners> cell_phi_value;
        VECTOR<int,num_corners> cell_phi_color;
        for(int i=0;i<num_corners;i++) cell_phi_value(i)=phi_value(it.index+bits(i));
        for(int i=0;i<num_corners;i++) cell_phi_color(i)=phi_color(it.index+bits(i));

        bool cut_cell=false,junction=false;
        const int base_color=cell_phi_color(0);
        for(int b=0;b<num_corners;b++) cut_cell|=(cell_phi_color(b)!=base_color);
        if(!cut_cell) continue; // uncut cell

        assert(!index_to_cell_data.Contains(it.index));
        HASH_CELL_DATA& cell_data=index_to_cell_data.Get_Or_Insert(it.index);
        HASH_CELL_INTERFACE& interface_cell_elements=cell_data.interface;
        VECTOR<HASH_CELL_BOUNDARY,TV::m>& boundary_cell_elements=cell_data.boundary;

        VECTOR<TV,num_corners> corners;
        grid.Cell_Domain(it.index).Corners(corners);
        Get_Hashed_Interface_Elements_For_Cell<T,TV_INT,TV,T_SURFACE,HASH_INTERFACE,HASH_CELL_INTERFACE>
            (it.index,corners,cell_phi_color,cell_phi_value,edge_vertices,face_vertices,
                cell_vertices,interface,particles,interface_cell_elements,junction);
        Get_Hashed_Boundary_Elements_For_Cell<T,TV_INT,TV,T_SURFACE,HASH_BOUNDARY,HASH_CELL_BOUNDARY>
            (it.index,corners,cell_phi_color,cell_phi_value,edge_vertices,face_vertices,
                node_vertices,boundary,particles,boundary_cell_elements);

        if(!junction) continue; // not a junction
        junction_cells.Insert(it.index);}

    // FIX AND SAVE MESHES
    Save_Mesh(index_to_cell_elements,grid,interface,boundary,index_to_cell_data,particles);
}
//#####################################################################
// Function Cut_Elements
//#####################################################################
template<class TV_INT,class TV,class T_FACE,class T_ELEMENT> void
Cut_Elements(ARRAY<ARRAY<T_ELEMENT>,TV_INT>& cut_elements,const ARRAY<T_ELEMENT>& elements,
    const RANGE<TV_INT>& range,const RANGE<TV>& domain)
{
    TV_INT size=range.Edge_Lengths();
    for(int a=0;a<TV::m;a++)
        if(size(a)>1){
            TV_INT new_size=size;
            new_size(a)/=2;
            ARRAY<T_FACE> t[2];
            ARRAY<T_ELEMENT> array[2];
            TV pt=domain.min_corner+domain.Edge_Lengths()*TV(new_size)/TV(size);
            typename BASIC_GEOMETRY_POLICY<TV>::HYPERPLANE plane(TV::Axis_Vector(a),pt);
            for(int i=0;i<elements.m;i++){
                for(int s=0;s<2;s++) t[s].Remove_All();
                T_FACE::Cut_With_Hyperplane(elements(i).face,plane,t[0],t[1],1e-14);
                for(int s=0;s<2;s++)
                    for(int j=0;j<t[s].m;j++){
                        array[s].Append(elements(i));
                        array[s].Last().face=t[s](j);}}
            RANGE<TV> domain0(domain),domain1(domain);
            RANGE<TV_INT> range0(range),range1(range);
            domain0.max_corner(a)=domain1.min_corner(a)=pt(a);
            range0.max_corner(a)=range1.min_corner(a)=range.min_corner(a)+new_size(a);
            Cut_Elements<TV_INT,TV,T_FACE,T_ELEMENT>(cut_elements,array[0],range0,domain0);
            Cut_Elements<TV_INT,TV,T_FACE,T_ELEMENT>(cut_elements,array[1],range1,domain1);
            return;}
    const TV_INT cell=range.min_corner;
    for(int i=0;i<elements.m;i++) cut_elements(cell).Append(elements(i));
}
//#####################################################################
// Function Save_Mesh
//#####################################################################
template<class TV> void MARCHING_CUBES_COLOR<TV>::
Save_Mesh(HASHTABLE<TV_INT,CELL_ELEMENTS>& index_to_cell_elements,const GRID<TV>& grid,const HASH_INTERFACE& interface,
    const HASH_BOUNDARY& boundary,const HASHTABLE<TV_INT,HASH_CELL_DATA>& index_to_cell_data,const GEOMETRY_PARTICLES<TV>& particles,
    const bool recut,const ARRAY<int>* const particle_dofs,const HASHTABLE<TV_INT>* const variable_cells)
{
    for(typename HASHTABLE<TV_INT,HASH_CELL_DATA>::CONST_ITERATOR it(index_to_cell_data);it.Valid();it.Next()){

        // HASH CELL DATA

        const TV_INT& cell_index=it.Key();
        const HASH_CELL_DATA& hash_cell_data=it.Data();
        const HASH_CELL_INTERFACE hash_cell_interface=hash_cell_data.interface;
        const VECTOR<HASH_CELL_BOUNDARY,TV::m> hash_cell_boundary=hash_cell_data.boundary;

        // IF CELL HAS VARIABLE PARTICLES THEN DO RECUTTING

        if(recut && variable_cells->Contains(cell_index)){

            // GET BOUNDARY AND INTERFACE ELEMENTS FORM HASHTABLES

            ARRAY<INTERFACE_ELEMENT> interface_elements;
            ARRAY<BOUNDARY_ELEMENT> boundary_elements;
            HASHTABLE<int> variable_indices;

            for(typename HASH_CELL_INTERFACE::CONST_ITERATOR it(hash_cell_interface);it.Valid();it.Next()){
                const VECTOR<int,2> color_pair=it.Key();
                const INTERVAL<int> interval=it.Data();
                const T_SURFACE& color_pair_interface=*interface.Get(color_pair);
                for(int e=interval.min_corner;e<interval.max_corner;e++){
                    const TV_INT& e_index=color_pair_interface.mesh.elements(e);
                    interface_elements.Add_End();
                    INTERFACE_ELEMENT& element=interface_elements.Last();
                    for(int i=0;i<TV::m;i++){
                        const int p_index=e_index(i);
                        element.face.X(i)=particles.X(p_index);
                        if((*particle_dofs)(p_index) && !variable_indices.Contains(p_index))
                            variable_indices.Insert(p_index);}
                    element.color_pair=color_pair;}}

            for(int v=0;v<TV::m;v++)
            for(typename HASH_CELL_BOUNDARY::CONST_ITERATOR it(hash_cell_boundary(v));it.Valid();it.Next()){
                const int color=it.Key();
                const INTERVAL<int> interval=it.Data();
                const T_SURFACE& color_boundary=*boundary.Get(color);
                for(int e=interval.min_corner;e<interval.max_corner;e++){
                    const TV_INT& e_index=color_boundary.mesh.elements(e);
                    boundary_elements.Add_End();
                    BOUNDARY_ELEMENT& element=boundary_elements.Last();
                    for(int i=0;i<TV::m;i++) element.face.X(i)=particles.X(e_index(i));
                    element.color=color;}}

            // DETERMINE THE RECUTTING RANGE

            ARRAY<TV> variable_particles;
            for(typename HASHTABLE<int>::ITERATOR it(variable_indices);it.Valid();it.Next())
                variable_particles.Append(particles.X(it.Key()));
            variable_particles.Append(grid.Center(cell_index));
            const RANGE<TV_INT> cell_range(grid.Clamp_To_Cell(RANGE<TV>::Bounding_Box(variable_particles)));
            const RANGE<TV> cell_domain(grid.Cell_Domain(cell_range));

            // CUT ELEMENTS

            ARRAY<ARRAY<INTERFACE_ELEMENT>,TV_INT> cut_interface(cell_range);
            ARRAY<ARRAY<BOUNDARY_ELEMENT>,TV_INT> cut_boundary(cell_range);
            
            Cut_Elements<TV_INT,TV,T_FACE,INTERFACE_ELEMENT>(cut_interface,interface_elements,cell_range,cell_domain);
            Cut_Elements<TV_INT,TV,T_FACE,BOUNDARY_ELEMENT>(cut_boundary,boundary_elements,cell_range,cell_domain);

            // COPY WITH CLAMPING

            const int h_min=cell_range.min_corner(TV::m-1);
            const int h_max=cell_range.max_corner(TV::m-1);
            for(RANGE_ITERATOR<TV::m> it(cell_range);it.Valid();it.Next()){

                // COPY ELEMENTS FOR CURRENT CELL

                const TV_INT& current_cell=it.index;
                const int current_h=it.index(TV::m-1); 

                CELL_ELEMENTS& cell_elements=index_to_cell_elements.Get_Or_Insert(current_cell);
                ARRAY<INTERFACE_ELEMENT>& cell_interface=cell_elements.interface;
                ARRAY<BOUNDARY_ELEMENT>& cell_boundary=cell_elements.boundary;
                
                const ARRAY<INTERFACE_ELEMENT>& current_cell_cut_interface=cut_interface(current_cell);
                const ARRAY<BOUNDARY_ELEMENT>& current_cell_cut_boundary=cut_boundary(current_cell);
                for(int i=0;i<current_cell_cut_boundary.m;i++) cell_boundary.Append(current_cell_cut_boundary(i));
                for(int i=0;i<current_cell_cut_interface.m;i++) cell_interface.Append(current_cell_cut_interface(i));

                // CLAMP ELEMENTS FROM OTHER CELLS

                const RANGE<TV> current_domain=grid.Cell_Domain(current_cell);
                for(int h=h_min;h<h_max;h++){

                    if(h==current_h) continue;
                    TV_INT other_cell(current_cell);
                    other_cell(TV::m-1)=h;

                    const ARRAY<INTERFACE_ELEMENT>& other_cell_cut_interface=cut_interface(other_cell);
                    const ARRAY<BOUNDARY_ELEMENT>& other_cell_cut_boundary=cut_boundary(other_cell);

                    for(int i=0;i<other_cell_cut_boundary.m;i++){
                        BOUNDARY_ELEMENT be=other_cell_cut_boundary(i);
                        for(int v=0;v<TV::m;v++) be.face.X(v)=current_domain.Clamp(be.face.X(v));
                        cell_boundary.Append(be);}
                    
                    for(int i=0;i<other_cell_cut_interface.m;i++){
                        INTERFACE_ELEMENT ie=other_cell_cut_interface(i);
                        for(int v=0;v<TV::m;v++) ie.face.X(v)=current_domain.Clamp(ie.face.X(v));
                        BOUNDARY_ELEMENT be0={ie.face,ie.color_pair(0)};
                        BOUNDARY_ELEMENT be1={ie.face,ie.color_pair(1)};
                        exchange(be0.face.X.x,be0.face.X.y);
                        cell_boundary.Append(be0);
                        cell_boundary.Append(be1);}}}}

        // IF NO MOVING PARTICLES THEN JUST COPY OVER TO NEW STRUCTURES

        else{

            CELL_ELEMENTS& cell_elements=index_to_cell_elements.Get_Or_Insert(cell_index);
            ARRAY<INTERFACE_ELEMENT>& cell_interface=cell_elements.interface;
            ARRAY<BOUNDARY_ELEMENT>& cell_boundary=cell_elements.boundary;

            for(typename HASH_CELL_INTERFACE::CONST_ITERATOR it(hash_cell_interface);it.Valid();it.Next()){
                const VECTOR<int,2> color_pair=it.Key();
                const INTERVAL<int> interval=it.Data();
                const T_SURFACE& color_pair_interface=*interface.Get(color_pair);
                for(int e=interval.min_corner;e<interval.max_corner;e++){
                    const TV_INT& e_index=color_pair_interface.mesh.elements(e);
                    cell_interface.Add_End();
                    INTERFACE_ELEMENT& element=cell_interface.Last();
                    for(int i=0;i<TV::m;i++) element.face.X(i)=particles.X(e_index(i));
                    element.color_pair=color_pair;}}

            for(typename HASH_CELL_BOUNDARY::CONST_ITERATOR it(hash_cell_boundary(TV::m-1));it.Valid();it.Next()){
                const int color=it.Key();
                const INTERVAL<int> interval=it.Data();
                const T_SURFACE& color_boundary=*boundary.Get(color);
                for(int e=interval.min_corner;e<interval.max_corner;e++){
                    const TV_INT& e_index=color_boundary.mesh.elements(e);
                    cell_boundary.Add_End();
                    BOUNDARY_ELEMENT& element=cell_boundary.Last();
                    for(int i=0;i<TV::m;i++) element.face.X(i)=particles.X(e_index(i));
                    element.color=color;}}}}
}
//#####################################################################
// Function Get_Elements_For_Cell
//#####################################################################
template<class TV> void MARCHING_CUBES_COLOR<TV>::
Get_Elements_For_Cell(ARRAY<INTERFACE_ELEMENT>& interface,ARRAY<BOUNDARY_ELEMENT>& boundary,
    const VECTOR<int,num_corners>& colors,const VECTOR<T,num_corners>& phi,const RANGE<TV>& cell_range)
{
    int si=interface.m,sb=boundary.m;
    Get_Elements_For_Cell(interface,boundary,colors,phi);
    for(int i=si;i<interface.m;i++)
        interface(i).face.X=interface(i).face.X*cell_range.Edge_Lengths()+cell_range.min_corner;
    for(int i=sb;i<boundary.m;i++)
        boundary(i).face.X=boundary(i).face.X*cell_range.Edge_Lengths()+cell_range.min_corner;
}
namespace PhysBAM{
template class MARCHING_CUBES_COLOR<VECTOR<float,2> >;
template class MARCHING_CUBES_COLOR<VECTOR<float,3> >;
template class MARCHING_CUBES_COLOR<VECTOR<double,2> >;
template class MARCHING_CUBES_COLOR<VECTOR<double,3> >;
}

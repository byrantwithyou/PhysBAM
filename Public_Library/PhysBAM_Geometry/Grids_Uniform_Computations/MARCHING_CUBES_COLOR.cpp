//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES_COLOR.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#define ENABLE_TIMING
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
#ifdef ENABLE_TIMING
#define rdtscll(val) do { \
     unsigned int __a,__d; \
     asm volatile("rdtsc" : "=a" (__a), "=d" (__d)); \
     (val) = ((unsigned long long)__a) | (((unsigned long long)__d)<<32); \
} while(0)
inline unsigned long long rdtsc(){unsigned long long x;rdtscll(x);return x;}
#endif

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

    int face_graph[64][2];
    for(int i=0;i<64;i++) for(int k=0;k<2;k++) face_graph[i][k]=-1;

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
    int test_index=interface_triangle_table.m;
    interface_triangle_table.Append(0);
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
template<class T> void
Get_Interface_Elements_For_Cell(ARRAY<TRIPLE<TRIANGLE_3D<T>,int,int> >& surface,const VECTOR<int,8>& re_color,
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
        TRIPLE<TRIANGLE_3D<T>,int,int> triple(triangle,color_list[GET_C(pat,0)],color_list[GET_C(pat,1)]);
        if(triple.y>triple.z){
            exchange(triple.y,triple.z);
            exchange(triple.x.X(1),triple.x.X(2));}
        surface.Append(triple);
    } while(!(pat&last_tri_bit));
}
//#####################################################################
// Function Get_Boundary_Elements_For_Cell
//#####################################################################
template<class T> void
Get_Boundary_Elements_For_Cell(ARRAY<PAIR<TRIANGLE_3D<T>,int> >& boundary,const int* re_color,
    const int* colors,const T* phi,int s,const int* color_list)
{
    typedef VECTOR<T,3> TV;
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
        PAIR<TRIANGLE_3D<T>,int> pair(triangle,color_list[GET_C(pat,0)]);
        if(s) exchange(pair.x.X.y,pair.x.X.z);
        boundary.Append(pair);
    } while(!(pat&last_tri_bit));
}
//#####################################################################
// Function Get_Interface_Elements_For_Cell
//#####################################################################
template<class T> void
Get_Interface_Elements_For_Cell(ARRAY<TRIPLE<SEGMENT_2D<T>,int,int> >& surface,const VECTOR<int,4>& re_color,
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
        TRIPLE<SEGMENT_2D<T>,int,int> triple(segment,color_list[GET_C(pat,0)],color_list[GET_C(pat,1)]);
        if(triple.y>triple.z){
            exchange(triple.y,triple.z);
            exchange(triple.x.X(0),triple.x.X(1));}
        surface.Append(triple);
    } while(!(pat&last_tri_bit));
}
//#####################################################################
// Function Get_Boundary_Elements_For_Cell
//#####################################################################
template<class T> void
Get_Boundary_Elements_For_Cell(ARRAY<PAIR<SEGMENT_2D<T>,int> >& boundary,const int* re_color,
    const int* colors,const T* phi,int s,const int* color_list)
{
    typedef VECTOR<T,2> TV;
    if(colors[0]==colors[1]){
        boundary.Append(PAIR<SEGMENT_2D<T>,int>(SEGMENT_2D<T>(TV(s,s),TV(1-s,s)),colors[0]));
        return;}

    T theta=phi[0]/(phi[0]+phi[1]);
    boundary.Append(PAIR<SEGMENT_2D<T>,int>(SEGMENT_2D<T>(TV(theta,s),TV(1-s,s)),colors[1-s]));
    boundary.Append(PAIR<SEGMENT_2D<T>,int>(SEGMENT_2D<T>(TV(s,s),TV(theta,s)),colors[s]));
}
}
//#####################################################################
// Function Initialize_Case_Table
//#####################################################################
template<class TV> void MARCHING_CUBES_COLOR<TV>::
Initialize_Case_Table()
{
#ifdef ENABLE_TIMING
    unsigned long long t0=rdtsc();
#endif
    if(TV::m==3) Initialize_Case_Table_3D();
    if(TV::m==2) Initialize_Case_Table_2D();
#ifdef ENABLE_TIMING
    unsigned long long t1=rdtsc();
    printf("setup: %.2f\n", (t1-t0)/3059.107);
#endif
}
//#####################################################################
// Function Get_Elements_For_Cell
//#####################################################################
template<class TV> void MARCHING_CUBES_COLOR<TV>::
Get_Elements_For_Cell(ARRAY<TRIPLE<T_FACE,int,int> >& surface,ARRAY<PAIR<T_FACE,int> >& boundary,
    const VECTOR<int,num_corners>& colors,const VECTOR<T,num_corners>& phi)
{
#ifdef ENABLE_TIMING
    unsigned long long t0=rdtsc();
#endif
    int next_color=0;
    HASHTABLE<int,int> color_map;
    VECTOR<int,num_corners> re_color;
    int color_list[8];
    for(int i=0;i<num_corners;i++)
        if(!color_map.Get(colors(i),re_color(i))){
            re_color(i)=next_color;
            color_list[next_color]=colors(i);
            color_map.Set(colors(i),next_color++);}

    Get_Interface_Elements_For_Cell(surface,re_color,colors,phi,color_list);
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
#ifdef ENABLE_TIMING
    unsigned long long t5=rdtsc();
    printf("query: %.2f\n", (t5-t0)/3059.107);
#endif
}
int vertex_lookup_2d[4][2]={{0,1},{2,3},{0,2},{1,3}};
//#####################################################################
// Function Get_Interface_Elements_For_Cell
//#####################################################################
template<class T,class TV_INT,class TV,class T_SURFACE> void
Get_Interface_Elements_For_Cell(const RANGE<TV>& range,HASHTABLE<VECTOR<int,2>,T_SURFACE*>& surface,const VECTOR<int,8>& re_color,
    const VECTOR<int,8>& colors,const VECTOR<T,8>& phi,const int* color_list,HASHTABLE<FACE_INDEX<3>,int>& edge_vertices,
    HASHTABLE<FACE_INDEX<3>,int>& face_vertices,HASHTABLE<TV_INT,int>& cell_vertices,const TV_INT& cell_index,
    GEOMETRY_PARTICLES<TV>& particles)
{
}
//#####################################################################
// Function Get_Interface_Elements_For_Cell
//#####################################################################
template<class T,class TV_INT,class TV,class T_SURFACE> void
Get_Interface_Elements_For_Cell(const RANGE<TV>& range,HASHTABLE<VECTOR<int,2>,T_SURFACE*>& surface,const VECTOR<int,4>& re_color,
    const VECTOR<int,4>& colors,const VECTOR<T,4>& phi,const int* color_list,HASHTABLE<FACE_INDEX<2>,int>& edge_vertices,
    HASHTABLE<FACE_INDEX<2>,int>& face_vertices,HASHTABLE<TV_INT,int>& cell_vertices,const TV_INT& cell_index,
    GEOMETRY_PARTICLES<TV>& particles)
{
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

    const VECTOR<VECTOR<int,2>,4>& bits=GRID<TV>::Binary_Counts(TV_INT());
    TV corners[4]={TV(range.min_corner.x,range.min_corner.y),TV(range.max_corner.x,range.min_corner.y),
        TV(range.min_corner.x,range.max_corner.y),TV(range.max_corner.x,range.max_corner.y)};
    int pts[5]={-1,-1,-1,-1,-1};
    int num=0;
    for(int e=0;e<4;e++){
        int v0=vertex_lookup_2d[e][0],v1=vertex_lookup_2d[e][1];
        if(colors[v0]!=colors[v1]){
            FACE_INDEX<TV::m> fi(e/2,cell_index+bits(v0));
            num++;
            if(!edge_vertices.Get(fi,pts[e])){
                int index=particles.Add_Element();
                pts[e]=index;
                edge_vertices.Set(fi,index);
                T theta=phi[v0]/(phi[v0]+phi[v1]);
                particles.X(index)=(1-theta)*corners[v0]+theta*corners[v1];}}}
    if(num>=3){
        if(!cell_vertices.Get(cell_index,pts[4])){
            int index=particles.Add_Element();
            pts[4]=index;
            cell_vertices.Set(cell_index,index);
            T total=0;
            TV X;
            for(int v=0;v<4;v++){
                total+=1/phi[v];
                X+=corners[v]/phi[v];}
            particles.X(index)=X/total;}}

    int pat;
    do{
        pat=interface_segment_table(seg++);
        VECTOR<int,2> c(color_list[GET_C(pat,0)],color_list[GET_C(pat,1)]);
        VECTOR<int,2> v(pts[GET_V(pat,0)],pts[GET_V(pat,1)]);
        printf("%i %i %i %i\n",GET_V(pat,0),GET_V(pat,1),v.x,v.y);
        if(c.x>c.y){
            exchange(c.x,c.y);
            exchange(v.x,v.y);}

        T_SURFACE* s=0;
        if(!surface.Get(c,s)){
            s=T_SURFACE::Create(particles);
            surface.Set(c,s);}
        s->mesh.elements.Append(v);
    } while(!(pat&last_tri_bit));
}
//#####################################################################
// Function Get_Elements_For_Cell
//#####################################################################
template<class TV> void MARCHING_CUBES_COLOR<TV>::
Get_Elements(const GRID<TV>& grid,HASHTABLE<VECTOR<int,2>,T_SURFACE*>& surface,HASHTABLE<int,T_SURFACE*>& boundary,
    const ARRAY<int,TV_INT>& color,const ARRAY<T,TV_INT>& phi)
{
    const int num_corners=1<<TV::m;
    HASHTABLE<FACE_INDEX<TV::m>,int> edge_vertices;
    HASHTABLE<FACE_INDEX<TV::m>,int> face_vertices;
    HASHTABLE<TV_INT,int> cell_vertices;
    GEOMETRY_PARTICLES<TV>& particles=*new GEOMETRY_PARTICLES<TV>;

    const VECTOR<VECTOR<int,TV::m>,(1<<TV::m)>& bits=GRID<TV>::Binary_Counts(TV_INT());
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        VECTOR<int,num_corners> cell_color;
        VECTOR<T,num_corners> cell_phi;
        for(int i=0;i<num_corners;i++) cell_color(i)=color(it.index+bits(i));
        for(int i=0;i<num_corners;i++) cell_phi(i)=phi(it.index+bits(i));
        int next_color=0;
        HASHTABLE<int,int> color_map;
        VECTOR<int,num_corners> re_color;
        int color_list[num_corners];
        for(int i=0;i<num_corners;i++)
            if(!color_map.Get(cell_color(i),re_color(i))){
                re_color(i)=next_color;
                color_list[next_color]=cell_color(i);
                color_map.Set(cell_color(i),next_color++);}
        Get_Interface_Elements_For_Cell(grid.Cell_Domain(it.index),surface,re_color,cell_color,cell_phi,color_list,edge_vertices,
            face_vertices,cell_vertices,it.index,particles);}
}
template class MARCHING_CUBES_COLOR<VECTOR<float,2> >;
template class MARCHING_CUBES_COLOR<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MARCHING_CUBES_COLOR<VECTOR<double,2> >;
template class MARCHING_CUBES_COLOR<VECTOR<double,3> >;
#endif

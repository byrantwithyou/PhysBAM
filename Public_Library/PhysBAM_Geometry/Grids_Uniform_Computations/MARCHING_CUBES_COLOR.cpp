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
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Math_Tools/clamp.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
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
// Function Get_Hashed_Elements_For_Cell
//#####################################################################
/*template<class T,class TV_INT,class TV,class T_SURFACE> void
Get_Hashed_Elements_For_Cell(const RANGE<TV>& range,HASHTABLE<VECTOR<int,2>,T_SURFACE*>& interface,HASHTABLE<int,T_SURFACE*>& boundary,
    const VECTOR<int,8>& re_color,const VECTOR<int,8>& colors,const VECTOR<T,8>& phi,const int* color_list,HASHTABLE<FACE_INDEX<3>,int>& edge_vertices,
    HASHTABLE<FACE_INDEX<3>,int>& face_vertices,HASHTABLE<TV_INT,int>& cell_vertices,HASHTABLE<TV_INT,int>& node_vertices,
    const TV_INT& cell_index,GEOMETRY_PARTICLES<TV>& particles,HASHTABLE<VECTOR<int,2>,INTERVAL<int> >& interface_cell_elements,
    HASHTABLE<int,INTERVAL<int> >& boundary_cell_elements)
{
    // RENAME COLORS

    int cs=0;
    for(int i=0;i<8;i++)
        cs=cs*(i+1)+re_color(i);
    if(!cs) return;

    TV corners[8]={TV(range.min_corner.x,range.min_corner.y,range.min_corner.z),TV(range.max_corner.x,range.min_corner.y,range.min_corner.z),
                   TV(range.min_corner.x,range.max_corner.y,range.min_corner.z),TV(range.max_corner.x,range.max_corner.y,range.min_corner.z),
                   TV(range.min_corner.x,range.min_corner.y,range.max_corner.z),TV(range.max_corner.x,range.min_corner.y,range.max_corner.z),
                   TV(range.min_corner.x,range.max_corner.y,range.max_corner.z),TV(range.max_corner.x,range.max_corner.y,range.max_corner.z)};

    // GET INTERFACE ELEMENTS

    {int tri=interface_case_table(cs);
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

    TV pts[19];
    T pts_phi[19];
    for(int a=0,k=0;a<3;a++){
        int mask=1<<a;
        for(int v=0;v<8;v++)
            if(!(v&mask)){
                T theta=phi(v)/(phi(v)+phi(v|mask));
                pts[k]=(1-theta)*corners[v]+theta*corners[v|mask];
                pts_phi[k++]=(1-theta)*phi(v)+theta*phi(v|mask);}}

    for(int i=0;i<7;i++){
        int a=averaging_order[i][0],b=averaging_order[i][1];
        T theta=pts_phi[a]/(pts_phi[a]+pts_phi[b]);
        pts[12+i]=(1-theta)*pts[a]+theta*pts[b];
        pts_phi[12+i]=(1-theta)*pts_phi[a]+theta*pts_phi[b];}

    int pat;
    int points[19]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    do{
        pat=interface_triangle_table(tri++);
        const TV_INT triangle(GET_V(pat,0),GET_V(pat,1),GET_V(pat,2));
        int tmp;
        for(int p=0;p<TV::m;p++){
            const int e=triangle(p);
            if(e<12){
                const int axis=e>>2;
                const int offset1=e&1;
                const int offset2=(e>>1)&1;
                FACE_INDEX<TV::m> fi(axis,cell_index+VECTOR<int,2>(offset1,offset2).Insert(0,axis));
                if(!edge_vertices.Get(fi,tmp)){
                    int index=particles.Add_Element();
                    points[e]=index;
                    edge_vertices.Set(fi,index);
                    particles.X(index)=pts[e];}
                else{
                    assert(points[e]==tmp || points[e]==-1);
                    points[e]=tmp;}}
            else if(e<18){
                const int axis=(e-12)>>1;
                const int sign=(e-12)&1;
                FACE_INDEX<TV::m> fi(axis,cell_index+TV_INT::Axis_Vector(axis)*sign);
                if(!face_vertices.Get(fi,tmp)){
                    int index=particles.Add_Element();
                    points[e]=index;
                    face_vertices.Set(fi,index);
                    particles.X(index)=pts[e];}
                else{
                    assert(points[e]==tmp || points[e]==-1);
                    points[e]=tmp;}}
            else{
                if(!cell_vertices.Get(cell_index,tmp)){
                    int index=particles.Add_Element();
                    points[e]=index;
                    cell_vertices.Set(cell_index,index);
                    particles.X(index)=pts[e];}
                else{
                    assert(points[e]==tmp || points[e]==-1);
                    points[e]=tmp;}}}

        TV_INT v(points[triangle(0)],points[triangle(2)],points[triangle(1)]);
        VECTOR<int,2> c(color_list[GET_C(pat,0)],color_list[GET_C(pat,1)]);
        if(c.x>c.y){
            exchange(c.x,c.y);
            exchange(v(1),v(2));}
        
        T_SURFACE* s=0;
        if(!interface.Get(c,s)){
            s=T_SURFACE::Create(particles);
            interface.Set(c,s);}
        s->mesh.elements.Append(v);
        interface_cell_elements.Get_Or_Insert(c,INTERVAL<int>(s->mesh.elements.m-1,s->mesh.elements.m)).max_corner=s->mesh.elements.m;
    }while(!(pat&last_tri_bit));}
    
    // GET BOUNDARY ELEMENTS

    for(int s=0;s<2;s++){
        const int offset=s<<2;

        int next_color=0;
        HASHTABLE<int,int> color_map;
        int color_list[8];
        int re_color2[4];

        for(int i=0;i<4;i++)
            if(!color_map.Get(colors(i+offset),re_color2[i])){
                re_color2[i]=next_color;
                color_list[next_color]=colors(i+offset);
                color_map.Set(colors(i+offset),next_color++);}

        cs=0;
        for(int i=0;i<4;i++)
            cs=cs*(i+1)+re_color2[i];

        int tri=boundary_case_table(cs);
        if(boundary_triangle_table(tri)&comparison_bit){
            int pat=boundary_triangle_table(tri);
            int amb0=(pat>>24)&7,amb1=(pat>>27)&7;
            if(color_list[amb1]<color_list[amb0])
                tri+=(pat>>12)&0x3f;
            else tri++;}
        
        const int num_pts=9;
        int points[num_pts]={-1,-1,-1,-1,-1,-1,-1,-1,-1};
        TV pts[num_pts]={corners[offset],corners[offset+1],corners[offset+2],corners[offset+3]};
        if(colors[offset+0]!=colors[offset+1]){T theta=phi[offset+0]/(phi[offset+0]+phi[offset+1]);pts[4]=(1-theta)*pts[0]+theta*pts[1];}
        if(colors[offset+2]!=colors[offset+3]){T theta=phi[offset+2]/(phi[offset+2]+phi[offset+3]);pts[5]=(1-theta)*pts[2]+theta*pts[3];}
        if(colors[offset+0]!=colors[offset+2]){T theta=phi[offset+0]/(phi[offset+0]+phi[offset+2]);pts[6]=(1-theta)*pts[0]+theta*pts[2];}
        if(colors[offset+1]!=colors[offset+3]){T theta=phi[offset+1]/(phi[offset+1]+phi[offset+3]);pts[7]=(1-theta)*pts[1]+theta*pts[3];}
        T total=0;
        for(int v=0;v<4;v++){
            total+=1/phi[offset+v];
            pts[8]+=pts[v]/phi[offset+v];}
        pts[8]/=total;
        
        int pat;
        do{
            pat=boundary_triangle_table(tri++);
            TV_INT triangle(GET_V(pat,0),GET_V(pat,1),GET_V(pat,2));
            int tmp;
            for(int p=0;p<TV::m;p++){
                const int e=triangle(p);
                if(e<4){
                    const int x_shift=e&1;
                    const int y_shift=(e>>1)&1;
                    TV_INT node_index=cell_index+TV_INT(x_shift,y_shift,s);
                    if(!node_vertices.Get(node_index,tmp)){
                        int index=particles.Add_Element();
                        points[e]=index;
                        node_vertices.Set(node_index,index);
                        particles.X(index)=pts[e];}
                    else{
                        assert(points[e]==tmp || points[e]==-1);
                        points[e]=tmp;}}
                else if(e<8){
                    const int axis=(e>>1)&1;
                    FACE_INDEX<TV::m> fi(axis,cell_index+TV_INT(e==7,e==5,s));
                    if(!edge_vertices.Get(fi,tmp)){
                        int index=particles.Add_Element();
                        points[e]=index;
                        edge_vertices.Set(fi,index);
                        particles.X(index)=pts[e];}
                    else{
                        assert(points[e]==tmp || points[e]==-1);
                        points[e]=tmp;}}
                else{
                    FACE_INDEX<TV::m> fi(TV::m-1,cell_index+TV_INT(0,0,s));
                    if(!face_vertices.Get(fi,tmp)){
                        int index=particles.Add_Element();
                        points[e]=index;
                        face_vertices.Set(fi,index);
                        particles.X(index)=pts[e];}
                    else{
                        assert(points[e]==tmp || points[e]==-1);
                        points[e]=tmp;}}}

            TV_INT v(points[triangle(0)],points[triangle(2)],points[triangle(1)]);
            int c=color_list[GET_C(pat,0)];
            if(s) exchange(v.y,v.z);
            
            T_SURFACE* surf=0;
            if(!boundary.Get(c,surf)){
                surf=T_SURFACE::Create(particles);
                boundary.Set(c,surf);}
            surf->mesh.elements.Append(v);
            boundary_cell_elements.Get_Or_Insert(c,INTERVAL<int>(surf->mesh.elements.m-1,surf->mesh.elements.m)).max_corner=surf->mesh.elements.m;
        } while(!(pat&last_tri_bit));}
}
//#####################################################################
// Function Get_Hashed_Elements_For_Cell
//#####################################################################
template<class T,class TV_INT,class TV,class T_SURFACE> void
Get_Hashed_Elements_For_Cell(const RANGE<TV>& range,HASHTABLE<VECTOR<int,2>,T_SURFACE*>& interface,HASHTABLE<int,T_SURFACE*>& boundary,
    const VECTOR<int,4>& re_color,const VECTOR<int,4>& colors,const VECTOR<T,4>& phi,const int* color_list,HASHTABLE<FACE_INDEX<2>,int>& edge_vertices,
    HASHTABLE<FACE_INDEX<2>,int>& face_vertices,HASHTABLE<TV_INT,int>& cell_vertices,HASHTABLE<TV_INT,int>& node_vertices,
    const TV_INT& cell_index,GEOMETRY_PARTICLES<TV>& particles,HASHTABLE<VECTOR<int,2>,INTERVAL<int> >& interface_cell_elements,
    HASHTABLE<int,INTERVAL<int> >& boundary_cell_elements)
{
    // RENAME COLORS

    int cs=0;
    for(int i=0;i<4;i++)
        cs=cs*(i+1)+re_color(i);
    if(!cs) return;

    TV corners[4]={TV(range.min_corner.x,range.min_corner.y),TV(range.max_corner.x,range.min_corner.y),
                   TV(range.min_corner.x,range.max_corner.y),TV(range.max_corner.x,range.max_corner.y)};

    // GET INTERFACE ELEMENTS

    {int seg=interface_case_table_2d(cs);
    if(interface_segment_table(seg)&comparison_bit){
        int pat=interface_segment_table(seg);
        int amb0=(pat>>24)&7,amb1=(pat>>27)&7;
        if(color_list[amb1]<color_list[amb0]) seg+=(pat>>12)&0x3f;
        else seg++;}

    const VECTOR<VECTOR<int,2>,4>& bits=GRID<TV>::Binary_Counts(TV_INT());
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
        if(c.x>c.y){
            exchange(c.x,c.y);
            exchange(v.x,v.y);}

        T_SURFACE* s=0;
        if(!interface.Get(c,s)){
            s=T_SURFACE::Create(particles);
            interface.Set(c,s);}
        s->mesh.elements.Append(v);
        interface_cell_elements.Get_Or_Insert(c,INTERVAL<int>(s->mesh.elements.m-1,s->mesh.elements.m)).max_corner=s->mesh.elements.m;
    } while(!(pat&last_tri_bit));}

    // GET BOUNDARY ELEMENTS

    for(int s=0;s<2;s++){
        const int offset=s<<1;
        int points[2]={-1,-1};
        TV pts[2]={corners[offset],corners[offset+1]};
        int tmp;
        for(int e=0;e<2;e++){
            TV_INT node_index=cell_index+TV_INT(e,s);
            if(!node_vertices.Get(node_index,tmp)){
                int index=particles.Add_Element();
                points[e]=index;
                node_vertices.Set(node_index,index);
                particles.X(index)=pts[e];}
            else{
                assert(points[e]==tmp || points[e]==-1);
                points[e]=tmp;}}
        if(colors[offset]==colors[offset+1]){
            VECTOR<int,2> v(points[s],points[1-s]);
            int c=colors[offset];

            T_SURFACE* surf=0;
            if(!boundary.Get(c,surf)){
                surf=T_SURFACE::Create(particles);
                boundary.Set(c,surf);}
            surf->mesh.elements.Append(v);
            boundary_cell_elements.Get_Or_Insert(c,INTERVAL<int>(surf->mesh.elements.m-1,surf->mesh.elements.m)).max_corner=surf->mesh.elements.m;}
        else{
            int points_m=-1;
            T theta=phi[offset]/(phi[offset]+phi[offset+1]);
            TV pts_m=pts[1]*theta+pts[0]*(1-theta);
            int tmp;
            FACE_INDEX<TV::m> fi(0,cell_index+TV_INT(0,s));
            if(!edge_vertices.Get(fi,tmp)){
                int index=particles.Add_Element();
                points_m=index;
                edge_vertices.Set(fi,index);
                particles.X(index)=pts_m;}
            else{
                assert(points_m==tmp || points_m==-1);
                points_m=tmp;}

            {VECTOR<int,2> v(points_m,points[1-s]);
            int c=colors[offset+1-s];
            
            T_SURFACE* surf=0;
            if(!boundary.Get(c,surf)){
                surf=T_SURFACE::Create(particles);
                boundary.Set(c,surf);}
            surf->mesh.elements.Append(v);
            boundary_cell_elements.Get_Or_Insert(c,INTERVAL<int>(surf->mesh.elements.m-1,surf->mesh.elements.m)).max_corner=surf->mesh.elements.m;}

            {VECTOR<int,2> v(points[s],points_m);
            int c=colors[offset+s];
            
            T_SURFACE* surf=0;
            if(!boundary.Get(c,surf)){
                surf=T_SURFACE::Create(particles);
                boundary.Set(c,surf);}
            surf->mesh.elements.Append(v);
            boundary_cell_elements.Get_Or_Insert(c,INTERVAL<int>(surf->mesh.elements.m-1,surf->mesh.elements.m)).max_corner=surf->mesh.elements.m;}}}
}*/

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
        VECTOR<int,2> v(p_index[GET_V(pat,0)],p_index[GET_V(pat,1)]);
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
    HASH_INTERFACE& interface,GEOMETRY_PARTICLES<TV>& particles,HASH_CELL_INTERFACE& interface_cell_elements,bool junction)
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

        TV_INT v(p_index[triangle(0)],p_index[triangle(2)],p_index[triangle(1)]);
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
    HASH_BOUNDARY& boundary,GEOMETRY_PARTICLES<TV>& particles,HASH_CELL_BOUNDARY& boundary_cell_elements)
{
    for(int s=0;s<2;s++){
        const int offset=s<<1;
        int points[2]={-1,-1};
        TV pts[2]={corners[offset],corners[offset+1]};
        int tmp;
        for(int e=0;e<2;e++){
            TV_INT node_index=cell_index+TV_INT(e,s);
            if(!node_vertices.Get(node_index,tmp)){
                int index=particles.Add_Element();
                points[e]=index;
                node_vertices.Set(node_index,index);
                particles.X(index)=pts[e];}
            else{
                assert(points[e]==tmp || points[e]==-1);
                points[e]=tmp;}}
        if(phi_color[offset]==phi_color[offset+1]){
            VECTOR<int,2> v(points[s],points[1-s]);
            int c=phi_color[offset];
            
            T_SURFACE* surf=0;
            if(!boundary.Get(c,surf)){
                surf=T_SURFACE::Create(particles);
                boundary.Set(c,surf);}
            surf->mesh.elements.Append(v);
            boundary_cell_elements.Get_Or_Insert(c,INTERVAL<int>(surf->mesh.elements.m-1,surf->mesh.elements.m)).max_corner=surf->mesh.elements.m;}
        else{
            int points_m=-1;
            T theta=phi_value[offset]/(phi_value[offset]+phi_value[offset+1]);
            TV pts_m=pts[1]*theta+pts[0]*(1-theta);
            int tmp;
            FACE_INDEX<TV::m> fi(0,cell_index+TV_INT(0,s));
            if(!edge_vertices.Get(fi,tmp)){
                int index=particles.Add_Element();
                points_m=index;
                edge_vertices.Set(fi,index);
                particles.X(index)=pts_m;}
            else{
                assert(points_m==tmp || points_m==-1);
                points_m=tmp;}
            
            {VECTOR<int,2> v(points_m,points[1-s]);
            int c=phi_color[offset+1-s];
            
            T_SURFACE* surf=0;
            if(!boundary.Get(c,surf)){
                surf=T_SURFACE::Create(particles);
                boundary.Set(c,surf);}
            surf->mesh.elements.Append(v);
            boundary_cell_elements.Get_Or_Insert(c,INTERVAL<int>(surf->mesh.elements.m-1,surf->mesh.elements.m)).max_corner=surf->mesh.elements.m;}

            {VECTOR<int,2> v(points[s],points_m);
            int c=phi_color[offset+s];
            
            T_SURFACE* surf=0;
            if(!boundary.Get(c,surf)){
                surf=T_SURFACE::Create(particles);
                boundary.Set(c,surf);}
            surf->mesh.elements.Append(v);
            boundary_cell_elements.Get_Or_Insert(c,INTERVAL<int>(surf->mesh.elements.m-1,surf->mesh.elements.m)).max_corner=surf->mesh.elements.m;}}}
}
//#####################################################################
// Function Get_Hashed_Boundary_Elements_For_Cell
//#####################################################################
template<class T,class TV_INT,class TV,class T_SURFACE,class HASH_BOUNDARY,class HASH_CELL_BOUNDARY> void
Get_Hashed_Boundary_Elements_For_Cell(const TV_INT& cell_index,const VECTOR<TV,8>& corners,const VECTOR<int,8>& phi_color,const VECTOR<T,8>& phi_value,
    HASHTABLE<FACE_INDEX<3>,int>& edge_vertices,HASHTABLE<FACE_INDEX<3>,int>& face_vertices,HASHTABLE<TV_INT,int>& node_vertices,
    HASH_BOUNDARY& boundary,GEOMETRY_PARTICLES<TV>& particles,HASH_CELL_BOUNDARY& boundary_cell_elements)
{
    for(int s=0;s<2;s++){
        const int offset=s<<2;
        
        int next_color=0;
        HASHTABLE<int,int> color_map;
        int color_list[8];
        int re_color2[4];
        
        for(int i=0;i<4;i++)
            if(!color_map.Get(phi_color(i+offset),re_color2[i])){
                re_color2[i]=next_color;
                color_list[next_color]=phi_color(i+offset);
                color_map.Set(phi_color(i+offset),next_color++);}
        
        int cs=0;
        for(int i=0;i<4;i++)
            cs=cs*(i+1)+re_color2[i];
        
        int tri=boundary_case_table(cs);
        if(boundary_triangle_table(tri)&comparison_bit){
            int pat=boundary_triangle_table(tri);
            int amb0=(pat>>24)&7,amb1=(pat>>27)&7;
            if(color_list[amb1]<color_list[amb0])
                tri+=(pat>>12)&0x3f;
            else tri++;}
        
        const int num_pts=9;
        int points[num_pts]={-1,-1,-1,-1,-1,-1,-1,-1,-1};
        TV pts[num_pts]={corners[offset],corners[offset+1],corners[offset+2],corners[offset+3]};
        if(phi_color[offset+0]!=phi_color[offset+1]){T theta=phi_value[offset+0]/(phi_value[offset+0]+phi_value[offset+1]);pts[4]=(1-theta)*pts[0]+theta*pts[1];}
        if(phi_color[offset+2]!=phi_color[offset+3]){T theta=phi_value[offset+2]/(phi_value[offset+2]+phi_value[offset+3]);pts[5]=(1-theta)*pts[2]+theta*pts[3];}
        if(phi_color[offset+0]!=phi_color[offset+2]){T theta=phi_value[offset+0]/(phi_value[offset+0]+phi_value[offset+2]);pts[6]=(1-theta)*pts[0]+theta*pts[2];}
        if(phi_color[offset+1]!=phi_color[offset+3]){T theta=phi_value[offset+1]/(phi_value[offset+1]+phi_value[offset+3]);pts[7]=(1-theta)*pts[1]+theta*pts[3];}
        T total=0;
        for(int v=0;v<4;v++){
            total+=1/phi_value[offset+v];
            pts[8]+=pts[v]/phi_value[offset+v];}
        pts[8]/=total;
        
        int pat;
        do{
            pat=boundary_triangle_table(tri++);
            TV_INT triangle(GET_V(pat,0),GET_V(pat,1),GET_V(pat,2));
            int tmp;
            for(int p=0;p<TV::m;p++){
                const int e=triangle(p);
                if(e<4){
                    const int x_shift=e&1;
                    const int y_shift=(e>>1)&1;
                    TV_INT node_index=cell_index+TV_INT(x_shift,y_shift,s);
                    if(!node_vertices.Get(node_index,tmp)){
                        int index=particles.Add_Element();
                        points[e]=index;
                        node_vertices.Set(node_index,index);
                        particles.X(index)=pts[e];}
                    else{
                        assert(points[e]==tmp || points[e]==-1);
                        points[e]=tmp;}}
                else if(e<8){
                    const int axis=(e>>1)&1;
                    FACE_INDEX<TV::m> fi(axis,cell_index+TV_INT(e==7,e==5,s));
                    if(!edge_vertices.Get(fi,tmp)){
                        int index=particles.Add_Element();
                        points[e]=index;
                        edge_vertices.Set(fi,index);
                        particles.X(index)=pts[e];}
                    else{
                        assert(points[e]==tmp || points[e]==-1);
                        points[e]=tmp;}}
                else{
                    FACE_INDEX<TV::m> fi(TV::m-1,cell_index+TV_INT(0,0,s));
                    if(!face_vertices.Get(fi,tmp)){
                        int index=particles.Add_Element();
                        points[e]=index;
                        face_vertices.Set(fi,index);
                        particles.X(index)=pts[e];}
                    else{
                        assert(points[e]==tmp || points[e]==-1);
                        points[e]=tmp;}}}

            TV_INT v(points[triangle(0)],points[triangle(2)],points[triangle(1)]);
            int c=color_list[GET_C(pat,0)];
            if(s) exchange(v.y,v.z);
            
            T_SURFACE* surf=0;
            if(!boundary.Get(c,surf)){
                surf=T_SURFACE::Create(particles);
                boundary.Set(c,surf);}
            surf->mesh.elements.Append(v);
            boundary_cell_elements.Get_Or_Insert(c,INTERVAL<int>(surf->mesh.elements.m-1,surf->mesh.elements.m)).max_corner=surf->mesh.elements.m;
        } while(!(pat&last_tri_bit));}
}
//#####################################################################
// Function Get_Elements
//#####################################################################
template<class TV> void MARCHING_CUBES_COLOR<TV>::
Get_Elements(const GRID<TV>& grid,HASH_INTERFACE& interface,HASH_BOUNDARY& boundary,HASH_CELL_TO_ELEMENT& cell_to_element,
    const ARRAY<int,TV_INT>& phi_color,const ARRAY<T,TV_INT>& phi_value,const int iterations,const bool verbose)
{
    // CELL CONSTANTS

    const VECTOR<VECTOR<int,TV::m>,(1<<TV::m)>& bits=GRID<TV>::Binary_Counts(TV_INT());

    // HASHTABLES FOR DIFFERENTLY CONSTRAINED PARTICLES

    HASHTABLE<FACE_INDEX<TV::m>,int> edge_vertices;
    HASHTABLE<FACE_INDEX<TV::m>,int> face_vertices;
    HASHTABLE<TV_INT,int> cell_vertices;
    HASHTABLE<TV_INT,int> node_vertices;

    // GLOBAL PARTICLES ARRAY AND JUNCTION CELL FLAGS

    GEOMETRY_PARTICLES<TV>& particles=*new GEOMETRY_PARTICLES<TV>;
    ARRAY<bool,TV_INT> junction_cell(grid.Domain_Indices());
    ARRAY<TV_INT> junction_cells_list;
    
    // INITIALIZE MESHES

    int fit_count=0;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){

        VECTOR<T,num_corners> cell_phi_value;
        VECTOR<int,num_corners> cell_phi_color;
        for(int i=0;i<num_corners;i++) cell_phi_value(i)=phi_value(it.index+bits(i));
        for(int i=0;i<num_corners;i++) cell_phi_color(i)=phi_color(it.index+bits(i));

        bool cut_cell=false,junction=false;
        const int base_color=cell_phi_color(0);
        for(int b=0;b<num_corners;b++) cut_cell|=(cell_phi_color(b)!=base_color);
        if(!cut_cell) continue; // uncut cell

        assert(!cell_to_element.Contains(it.index));
        PAIR<HASH_CELL_INTERFACE,HASH_CELL_BOUNDARY>& cell_elements=cell_to_element.Get_Or_Insert(it.index);
        HASH_CELL_INTERFACE& interface_cell_elements=cell_elements.x;
        HASH_CELL_BOUNDARY& boundary_cell_elements=cell_elements.y;

        VECTOR<TV,num_corners> corners;
        grid.Cell_Domain(it.index).Corners(corners);
        Get_Hashed_Interface_Elements_For_Cell<T,TV_INT,TV,T_SURFACE,HASH_INTERFACE,HASH_CELL_INTERFACE>
            (it.index,corners,cell_phi_color,cell_phi_value,edge_vertices,face_vertices,
                cell_vertices,interface,particles,interface_cell_elements,junction);
        Get_Hashed_Boundary_Elements_For_Cell<T,TV_INT,TV,T_SURFACE,HASH_BOUNDARY,HASH_CELL_BOUNDARY>
            (it.index,corners,cell_phi_color,cell_phi_value,edge_vertices,face_vertices,
                node_vertices,boundary,particles,boundary_cell_elements);

        if(!junction) continue; // not a junction
        junction_cell(it.index)=true;
        junction_cells_list.Append(it.index);
        for(typename HASH_CELL_INTERFACE::CONST_ITERATOR it(interface_cell_elements);it.Valid();it.Next()) fit_count++;}

    // DOFS AND NORMALS

    ARRAY<int> particle_dofs(particles.number);
    ARRAY<TV> normals(fit_count);
    ARRAY<TV> midpoints(fit_count);
    
    // INITIALIZE DOFS

    for(typename HASHTABLE<TV_INT,int>::ITERATOR it(cell_vertices);it.Valid();it.Next()){
        particle_dofs(it.Data())=(1<<TV::m)-1;}

    for(typename HASHTABLE<FACE_INDEX<TV::m>,int>::ITERATOR it(face_vertices);it.Valid();it.Next()){
        particle_dofs(it.Data())=(1<<TV::m)-1-(1<<it.Key().axis);}

    for(typename HASHTABLE<FACE_INDEX<TV::m>,int>::ITERATOR it(edge_vertices);it.Valid();it.Next()){
        RANGE<TV_INT> range(RANGE<TV_INT>::Centered_Box());
        range.max_corner(it.Key().axis)=2;
        bool trusted=true;
        for(RANGE_ITERATOR<TV::m> it2(range);it2.Valid();it2.Next()){
            TV_INT index=it2.index+it.Key().index;
            if(junction_cell.Valid_Index(index)&&junction_cell(index)){
                trusted=false;break;}}
        if(!trusted) particle_dofs(it.Data())=1<<it.Key().axis;}

    // PARTICLES TO DOFS CORRESPONDENCE

    ARRAY<int> index_map;                          // maps dofs to particles
    ARRAY<int> reverse_index_map(particle_dofs.m); // maps particles to dofs
    for(int i=0;i<particle_dofs.m;i++){
        reverse_index_map(i)=-1;
        if(particle_dofs(i))
            reverse_index_map(i)=index_map.Append(i);}

    // NORMALS TO PARTICLES CORRESPONDENCE

    ARRAY<ARRAY<int> > particle_to_fit(index_map.m);
    ARRAY<ARRAY<int> > fit_to_particle(fit_count);

    // INITIALIZE NORMALS TO PARTICLES CORRESPONDENCE

    int fit_index=-1;
    for(int c=0;c<junction_cells_list.m;c++){
        const TV_INT& junction_cell_index=junction_cells_list(c);
        const HASH_CELL_INTERFACE& junction_interface_cell_elements=cell_to_element.Get(junction_cell_index).x;
        for(typename HASH_CELL_INTERFACE::CONST_ITERATOR it(junction_interface_cell_elements);it.Valid();it.Next()){
            fit_index++;
            HASHTABLE<int> particle_indices_ht;
            const VECTOR<int,2>& color_pair=it.Key();
            const T_SURFACE& color_pair_interface=*interface.Get(color_pair);
            RANGE<TV_INT> range(RANGE<TV_INT>::Centered_Box()*2);range.max_corner+=1;
            for(RANGE_ITERATOR<TV::m> it2(range);it2.Valid();it2.Next()){
                const TV_INT& cell_index=junction_cell_index+it2.index;
                if(cell_to_element.Contains(cell_index)){
                    const HASH_CELL_INTERFACE& interface_cell_elements=cell_to_element.Get(cell_index).x;
                    INTERVAL<int> interval;
                    if(interface_cell_elements.Get(color_pair,interval))
                        for(int i=interval.min_corner;i<interval.max_corner;i++){
                            TV_INT element=color_pair_interface.mesh.elements(i);
                            for(int j=0;j<TV::m;j++)
                                if(!particle_indices_ht.Contains(element(j)))
                                    particle_indices_ht.Insert(element(j));}}}
            ARRAY<int> particle_indices;
            for(typename HASHTABLE<int>::ITERATOR it(particle_indices_ht);it.Valid();it.Next()){
                fit_to_particle(fit_index).Append(it.Key());
                const int reduced_index=reverse_index_map(it.Key());
                if(reduced_index>=0) particle_to_fit(reduced_index).Append(fit_index);}}}
    PHYSBAM_ASSERT(fit_index==fit_count-1); // make sure all normals got processed 

    // INITIALIZE PARTICLE ADJACENCY MAP

    ARRAY<ARRAY<int> > adjacency_map(index_map.m);    
    for(typename HASH_INTERFACE::ITERATOR it(interface);it.Valid();it.Next()){
        T_SURFACE& surf=*it.Data();
        surf.Update_Number_Nodes();
        surf.mesh.Initialize_Adjacent_Elements();
        const ARRAY<ARRAY<int> >& adjacent_elements=*surf.mesh.adjacent_elements;
        for(int i=0;i<adjacent_elements.m;i++){
            for(int j=0;j<adjacent_elements(i).m;j++){
                int k=adjacent_elements(i)(j);
                if(i<k){
                    VECTOR<int,TV::m+1> nodes;
                    TV_INT ei=surf.mesh.elements(i);
                    TV_INT ek=surf.mesh.elements(k);
                    int u=-1;
                    for(int m=0;m<TV::m;m++)if(!ek.Contains(ei(m))){u=m;break;}
                    for(int m=0;m<TV::m;m++){nodes(m)=ei(u++);if(u==TV::m) u=0;}
                    nodes(TV::m)=ek.Sum()-ei.Sum()+nodes(0);
                    if(TV::m==2){
                        int p=reverse_index_map(nodes(1));
                        if(p!=-1){
                            adjacency_map(p).Append(nodes(0));
                            adjacency_map(p).Append(nodes(2));}}
                    else if(TV::m==3){
                        int p1=reverse_index_map(nodes(1));
                        int p2=reverse_index_map(nodes(2));
                        if(p1!=-1) adjacency_map(p1).Append(nodes(2));
                        if(p2!=-1) adjacency_map(p2).Append(nodes(1));}
                    else PHYSBAM_FATAL_ERROR();}}}}

    // ADJUST THE MESH

    if(verbose) LOG::cout<<"Adjusting mesh...";
    for(int i=0;i<iterations*2;i++){
        if(i&1){ // update positions
            for(int p=0;p<index_map.m;p++){
                const ARRAY<int>& fit_array=particle_to_fit(p);
                const int index=index_map(p);
                TV dx; TV& x=particles.X(index);
                for(int f=0;f<fit_array.m;f++){
                    const int fit=fit_array(f);
                    const TV& n=normals(fit);
                    dx+=n.Dot(x-midpoints(fit))*n;}
                    dx/=fit_array.m; x-=dx;}
            for(int p=0;p<index_map.m;p++){
                const int index=index_map(p);
                ARRAY<int>& adj=adjacency_map(p);
                const int dof=particle_dofs(index);
                if(adj.m&&(dof==1||dof==2||dof==4)){
                    TV dx; TV& x=particles.X(index);
                    for(int i=0;i<adj.m;i++)
                        dx+=particles.X(adj(i));
                    dx/=adj.m; dx-=x; dx*=0.1; x+=dx;}}}
        else{ // update best fit manifolds
            for(int f=0;f<fit_count;f++){
                const ARRAY<int>& particle_array=fit_to_particle(f);
                midpoints(f)=particles.X.Subset(particle_array).Average();
                SYMMETRIC_MATRIX<T,TV::m> fit_matrix;
                for(int j=0;j<particle_array.m;j++)
                    fit_matrix+=SYMMETRIC_MATRIX<T,TV::m>::Outer_Product(particles.X(particle_array(j))-midpoints(f));
                DIAGONAL_MATRIX<T,TV::m> eigenvalues;
                MATRIX<T,TV::m> eigenvectors;
                fit_matrix.Solve_Eigenproblem(eigenvalues,eigenvectors);
                normals(f)=eigenvectors.Column(eigenvalues.To_Vector().Arg_Min());}}
        if(verbose) LOG::cout<<".";}
    if(verbose) LOG::cout<<std::endl;

    // DRAW DEBUG PARTICLES

    for(int p=0;p<particle_dofs.m;p++){
        switch(particle_dofs(p)){
            case 0: Add_Debug_Particle(particles.X(p),VECTOR<T,3>(.5,.5,.5)); break;
            case 1: case 2: case 4: Add_Debug_Particle(particles.X(p),VECTOR<T,3>(1,1,1)); break;
            case 3: case 5: case 6: Add_Debug_Particle(particles.X(p),VECTOR<T,3>(1,1,0)); break;
            case 7: Add_Debug_Particle(particles.X(p),VECTOR<T,3>(1,0,0)); break;
            default: PHYSBAM_FATAL_ERROR();}}
    for(int f=0;f<fit_count;f++){
        const TV& normal_vector=normals(f);
        Add_Debug_Particle(midpoints(f),VECTOR<T,3>(1,0,1));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,normal_vector);
        Add_Debug_Particle(midpoints(f),VECTOR<T,3>(1,0,1));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,-normal_vector);}
}
template class MARCHING_CUBES_COLOR<VECTOR<float,2> >;
template class MARCHING_CUBES_COLOR<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MARCHING_CUBES_COLOR<VECTOR<double,2> >;
template class MARCHING_CUBES_COLOR<VECTOR<double,3> >;
#endif

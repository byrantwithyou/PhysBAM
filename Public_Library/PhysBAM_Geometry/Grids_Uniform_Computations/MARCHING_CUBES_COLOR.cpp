//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES_COLOR.h>
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
const int pts_bit=1<<29;
/* Triangle encoding
 * 0 e 000000000 xxx yyy aaaaa bbbbb ccccc   (triangle)
 * 1 0 aaa bbb ccc ddd ssssss ssssss ssssss (comparison)
 *
 * indexing: [0-11 = edges] [12-17; axis*2+side] [18 center]
 * faces: [0-3 = corners] [edges: 1 0 3 4] [8 center]
 * 
 * 0 0 1 0000 x aaaa bbbb cccc dddd eeee ffff
 * 
 */
// start, end, inside color, outside color (c0<c1)
struct EDGE {int v0,v1,c0,c1;void Flip(){exchange(c0,c1);exchange(v0,v1);}};
EDGE center_edges[6][4]; // c0,c1 are vertex to take color from
struct DOMINANT_PAIR {EDGE a,b;int dc;} dominant_pair[6][2];
const int permute_rx[19]={1,3,0,2,10,11,8,9,4,5,6,7,12,13,16,17,15,14,18};
const int permute_ry[19]={8,10,9,11,6,4,7,5,2,0,3,1,17,16,14,15,12,13,18};
const int permute_rx_corners[8]={2,3,6,7,0,1,4,5};
const int permute_ry_corners[8]={4,0,6,2,5,1,7,3};
int face_edges[6][4];
const int greedy=2;

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
bool Merge_Edges(EDGE& e0,EDGE& e1)
{
    if(e0.c0==e1.c1){e0.c0=e1.c0;return true;}
    if(e0.c1==e1.c0){e0.c1=e1.c1;return true;}
    return false;
}
void Insert_Face_Graph_Edge(int (*face_graph)[2],EDGE* edges,int e)
{
    if(edges[e].v0>edges[e].v1) edges[e].Flip();
    int* fg=face_graph[(1<<edges[e].v0)|(1<<edges[e].v1)];
    if(fg[0]==-1){fg[0]=e;return;}
    if(fg[1]!=-1){
        EDGE& e0=edges[fg[0]],&e1=edges[fg[1]];
        int s=Merge_Edges(e1,edges[e]);
        PHYSBAM_ASSERT(s);
        if(e1.c0==e1.c1){fg[1]=-1;return;}
        int t=Merge_Edges(e0,e1);
        PHYSBAM_ASSERT(t);
        fg[1]=-1;
        if(e0.c0==e0.c1){fg[0]=-1;return;}
        return;}
    EDGE& e0=edges[fg[0]];
    if(!Merge_Edges(e0,edges[e])){fg[1]=e;return;}
    if(e0.c0==e0.c1) fg[0]=-1;
}
void Emit_Loop_Triangles(int* vertices,int n,int c0,int c1)
{
    for(int i=0;i<n-2;i++){
        interface_triangle_table.Append((c0<<18)|(c1<<15)|(vertices[2*i]<<10)|(vertices[2*i+1]<<5)|vertices[2*i+2]);
        vertices[n+i]=vertices[2*i];}
}
#if 0
int edge_to_face_mask[12]={0x14,0x18,0x24,0x28,0x11,0x12,0x21,0x22,0x05,0x06,0x09,0x0a};

// 0 = no intersection, 1 = can intersect, 2 = coplanar
int Add_Face_Edges(int a,int b,int c,int d)
{
    // Share an endpoint
    if(a>b) exchange(a,b);
    if(c>d) exchange(c,d);
    if(a>c){exchange(a,c);exchange(b,d);}
    if(b==c && b==d) return 0;

    // Edge on cube face
    if((edge_to_face_mask[a]&edge_to_face_mask[b]) || (edge_to_face_mask[c]&edge_to_face_mask[d])) return 0;

    
    return 88;
}
#endif
//#####################################################################
// Function Emit_Interface_Triangles
//#####################################################################
void Emit_Interface_Triangles(int* colors,int color_hint)
{
    EDGE edges[25]; // need an extra one
    int num_edges=0,table_size=interface_triangle_table.m;
    for(int a=0;a<3;a++)
        for(int b=a+1;b<3;b++)
            for(int s=0;s<2;s++)
                Add_Face_Edges(colors,edges,num_edges,a,b,s,color_hint);

    int adj[12][2];
    for(int i=0;i<12;i++) for(int j=0;j<2;j++) adj[i][j]=-1;
    for(int i=0;i<num_edges;i++){
        if(edges[i].v0<12) adj[edges[i].v0][1]=i;
        if(edges[i].v1<12) adj[edges[i].v1][0]=i;}

    long long pt_mask_helper=0;
    int pt_mask=pts_bit;
    for(int i=0;i<num_edges;i++){
        if(edges[i].v0>=12)
            pt_mask_helper|=1LL<<(12*((edges[i].v0-12)/2)+edges[i].v1);
        if(edges[i].v1>=12)
            pt_mask_helper|=1LL<<(12*((edges[i].v1-12)/2)+edges[i].v0);}

    for(int f=0;f<6;f++)
        for(int i=0;i<4;i++)
            if(pt_mask_helper&(1LL<<((f/2)*12+face_edges[f][i])))
                pt_mask|=1<<(4*f+i);
    int mask_tri_index=-1;
    if(pt_mask!=pts_bit)
        mask_tri_index=interface_triangle_table.Append(pt_mask);

    int face_graph[64][2];
    for(int i=0;i<64;i++) for(int k=0;k<2;k++) face_graph[i][k]=-1;

    if(greedy==0){
        for(int i=0;i<12;i++){
            if(adj[i][0]==-1) continue;
            EDGE& in=edges[adj[i][0]],&out=edges[adj[i][1]];
            PHYSBAM_ASSERT(in.c0==out.c0 && in.c1==out.c1);
            if(in.v0==out.v1){
                adj[in.v0][0]=-1;
                adj[in.v0][1]=-1;
                adj[i][0]=-1;
                adj[i][1]=-1;
                continue;}
            interface_triangle_table.Append((in.c0<<18)|(in.c1<<15)|(in.v0<<10)|(in.v1<<5)|out.v1);
            in.v1=out.v1;
            if(in.v1<12) adj[in.v1][0]=adj[i][0];
            else if(in.v0>=12){
                in.v0-=12;
                in.v1-=12;
                Insert_Face_Graph_Edge(face_graph,edges,adj[i][0]);}
            adj[i][0]=-1;
            adj[i][1]=-1;}}
    else if(greedy==1){
        int vertices[30];
        for(int i=0;i<12;i++){
            if(adj[i][0]==-1) continue;
            EDGE& in=edges[adj[i][0]],&out=edges[adj[i][1]];
            PHYSBAM_ASSERT(in.c0==out.c0 && in.c1==out.c1);
            vertices[10]=i;
            int L=9,R=11;
            for(int w=i;;L--){
                int e=adj[w][0];
                adj[w][0]=-1;
                assert(edges[e].v1==w);
                w=edges[e].v0;
                vertices[L]=w;
                assert(w>=0);
                if(w==i || w>=12) break;}
            if(vertices[L]==i){
                Emit_Loop_Triangles(vertices+L,10-L+1,in.c0,in.c1);
                continue;}
            for(int v=i;;R++){
                int e=adj[v][1];
                adj[v][0]=-1;
                assert(edges[e].v0==v);
                v=edges[e].v1;
                vertices[R]=v;
                assert(v>=0);
                if(v>=12) break;}
            int M=(L+R)/2;
            interface_triangle_table.Append((in.c0<<18)|(in.c1<<15)|(vertices[L]<<10)|(vertices[M]<<5)|vertices[R]);
            out.v0=vertices[L]-12;
            out.v1=vertices[R]-12;
            Insert_Face_Graph_Edge(face_graph,edges,adj[i][1]);
            Emit_Loop_Triangles(vertices+M,R-M+1,in.c0,in.c1);
            Emit_Loop_Triangles(vertices+L,M-L+1,in.c0,in.c1);}}
    else{
        struct CURVE
        {
            int vertices[30];
            int c0,c1;
            int L,R;
        };
        int degree[6];
        CURVE curves[20];
        int curve_cnt=0;
        for(int i=0;i<12;i++){
            if(adj[i][0]==-1) continue;
            EDGE& in=edges[adj[i][0]],&out=edges[adj[i][1]];
            PHYSBAM_ASSERT(in.c0==out.c0 && in.c1==out.c1);
            curves[curve_cnt].vertices[10]=i;
            curves[curve_cnt].L=9;
            curves[curve_cnt].R=11;
            for(int w=i;;curves[curve_cnt].L--){
                int e=adj[w][0];
                adj[w][0]=-1;
                assert(edges[e].v1==w);
                w=edges[e].v0;
                curves[curve_cnt].vertices[curves[curve_cnt].L]=w;
                assert(w>=0);
                if(w==i || w>=12) break;}
            if(curves[curve_cnt].vertices[curves[curve_cnt].L]==i){
                Emit_Loop_Triangles(curves[curve_cnt].vertices+curves[curve_cnt].L,10-curves[curve_cnt].L+1,in.c0,in.c1);
                continue;}
            for(int v=i;;curves[curve_cnt].R++){
                int e=adj[v][1];
                adj[v][0]=-1;
                assert(edges[e].v0==v);
                v=edges[e].v1;
                curves[curve_cnt].vertices[curves[curve_cnt].R]=v;
                assert(v>=0);
                if(v>=12) break;}
            curves[curve_cnt].c0=in.c0;
            curves[curve_cnt].c1=in.c1;
            out.v0=curves[curve_cnt].vertices[curves[curve_cnt].L]-12;
            out.v1=curves[curve_cnt].vertices[curves[curve_cnt].R]-12;
            Insert_Face_Graph_Edge(face_graph,edges,adj[i][1]);
            degree[curves[curve_cnt].vertices[curves[curve_cnt].L]-12]=max(degree[curves[curve_cnt].vertices[curves[curve_cnt].L]-12],curves[curve_cnt].R-curves[curve_cnt].L);
            degree[curves[curve_cnt].vertices[curves[curve_cnt].R]-12]=max(degree[curves[curve_cnt].vertices[curves[curve_cnt].R]-12],curves[curve_cnt].R-curves[curve_cnt].L);
            curve_cnt++;}
        
        for(int c=0;c<curve_cnt;c++){
            CURVE& curve=curves[c];
            int M=(curve.L+curve.R)/2;
            int i;
            for(i=curve.L+2;i<=M;i++){
                interface_triangle_table.Append((curve.c0<<18)|(curve.c1<<15)|(curve.vertices[curve.L]<<10)|(curve.vertices[i-1]<<5)|curve.vertices[i]);
                interface_triangle_table.Append((curve.c0<<18)|(curve.c1<<15)|(curve.vertices[curve.R]<<10)|(curve.vertices[curve.R-(i-curve.L)]<<5)|curve.vertices[curve.R-(i-curve.L)+1]);}
            i--;
            if((i-curve.L)==(curve.R-i))
                interface_triangle_table.Append((curve.c0<<18)|(curve.c1<<15)|(curve.vertices[curve.L]<<10)|(curve.vertices[i]<<5)|curve.vertices[curve.R]);
            else{
                if((degree[curve.vertices[curve.L]-12]>degree[curve.vertices[curve.R]-12])||((degree[curve.vertices[curve.L]-12]==degree[curve.vertices[curve.R]-12])&&(curve.vertices[curve.L]>curve.vertices[curve.R]))){
                    interface_triangle_table.Append((curve.c0<<18)|(curve.c1<<15)|(curve.vertices[curve.L]<<10)|(curve.vertices[i]<<5)|curve.vertices[i+1]);
                    interface_triangle_table.Append((curve.c0<<18)|(curve.c1<<15)|(curve.vertices[curve.L]<<10)|(curve.vertices[i+1]<<5)|curve.vertices[curve.R]);}
                else{
                    interface_triangle_table.Append((curve.c0<<18)|(curve.c1<<15)|(curve.vertices[curve.L]<<10)|(curve.vertices[i]<<5)|curve.vertices[curve.R]);
                    interface_triangle_table.Append((curve.c0<<18)|(curve.c1<<15)|(curve.vertices[i]<<10)|(curve.vertices[i+1]<<5)|curve.vertices[curve.R]);
                }}}}

    int face_graph_copy[64][2];
    for(int i=0;i<64;i++) for(int s=0;s<2;s++) face_graph_copy[i][s]=face_graph[i][s];

    bool progress=true,make_center=false;
    while(progress){
        progress=false;
        make_center=false;
        for(int i=0;i<6;i++)
            for(int j=0;j<6;j++)
                for(int m=0;m<2;m++){
                    int* fgij=face_graph[(1<<i)|(1<<j)],fg0=fgij[m];
                    if(fg0>=0){
                        make_center=true;
                        for(int k=i+1;k<6;k++)
                            for(int n=0;n<2;n++){
                                int* fgjk=face_graph[(1<<j)|(1<<k)],fg1=fgjk[n];
                                if(fg1>=0){
                                    if(((1<<edges[fg0].c0)|(1<<edges[fg0].c1))==((1<<edges[fg1].c0)|(1<<edges[fg1].c1))){
                                        if(edges[fg1].v0!=j) edges[fg1].Flip();
                                        interface_triangle_table.Append((edges[fg1].c0<<18)|(edges[fg1].c1<<15)|((i+12)<<10)|((j+12)<<5)|(k+12));
                                        if(edges[fg0].v1==j) edges[fg0].v1=k;
                                        else edges[fg0].v0=k;
                                        Insert_Face_Graph_Edge(face_graph,edges,fg0);
                                        fgij[m]=-1;
                                        if(m==0 && fgij[1]>=0) exchange(fgij[0],fgij[1]);
                                        fgjk[n]=-1;
                                        if(m==0 && fgjk[1]>=0) exchange(fgjk[0],fgjk[1]);
                                        progress=true;}}}}}}

    if(make_center){
        interface_triangle_table(mask_tri_index)=pt_mask|(1<<24);
        for(int i=0;i<6;i++)
            for(int j=0;j<6;j++){
                PHYSBAM_ASSERT(face_graph_copy[(1<<i)|(1<<j)][1]==-1);
                int fg=face_graph_copy[(1<<i)|(1<<j)][0];
                if(fg==-1) continue;
                interface_triangle_table.Append((edges[fg].c0<<18)|(edges[fg].c1<<15)|((i+12)<<10)|((j+12)<<5)|18);}}

    PHYSBAM_ASSERT(table_size!=interface_triangle_table.m);
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
// Function Enumerate_Interface_Cases_3D
//#####################################################################
void Enumerate_Interface_Cases_3D(int* colors, int i, int mc, int cs)
{
    if(i==8){
        interface_case_table(cs)=interface_triangle_table.m;
        if(cs) Emit_Interface_Triangles(colors);
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
#define EMIT(v0,v1,v2,c) boundary_triangle_table.Append((c<<18)|((v0)<<10)|((v1)<<5)|(v2))
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
#define EMIT(v0,v1,c0,c1) interface_segment_table.Append((c0<<18)|(c1<<15)|(v0<<10)|(v1<<5))
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
//#####################################################################
// Function Get_Interface_Elements_For_Cell
//#####################################################################
template<class T> void
Get_Interface_Elements_For_Cell(ARRAY<TRIPLE<TRIANGLE_3D<T>,int,int> >& surface,const VECTOR<int,8>& re_color,
    const VECTOR<int,8>& colors,const VECTOR<T,8>& phi,const HASHTABLE<int,int>& color_map,const int* color_list)
{
    typedef VECTOR<T,3> TV;
    int cs=0;
    for(int i=0;i<8;i++)
        cs=cs*(i+1)+re_color(i);
    if(!cs) return;

    int tri=interface_case_table(cs);
    if(interface_triangle_table(tri)&comparison_bit){
        int pat=interface_triangle_table(tri);
        if((pat>>18)&0x3f){
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

    const int num_pts=19;
    VECTOR<VECTOR<int,3>,8> bits=GRID<TV>::Binary_Counts(VECTOR<int,3>());
    TV pts[num_pts];
    for(int a=0,k=0;a<3;a++){
        int mask=1<<a;
        for(int v=0;v<8;v++)
            if(!(v&mask)){
                if(colors(v)!=colors(v|mask)){
                    T theta=phi(v)/(phi(v)+phi(v|mask));
                    pts[k]=(1-theta)*TV(bits(v))+theta*TV(bits(v|mask));}
                k++;}}

    int pt_mask=interface_triangle_table(tri);
    if(pt_mask&pts_bit){
        tri++;
        int num_center_points=0;
        if(0){
            for(int f=0;f<6;f++){
                int bits=(pt_mask>>(f*4))&0xf;
                if(!bits) continue;
                for(int i=0;i<4;i++)
                    if(bits&(1<<i))
                        pts[12+f]+=pts[face_edges[f][i]];
                pts[12+f]/=(3+(bits==0xf));
                if(pt_mask&(1<<24)){
                    num_center_points++;
                    pts[18]+=pts[12+f];}}
            if(pt_mask&(1<<24))
                pts[18]/=num_center_points;}
        else{
            for(int a=0;a<3;a++){
                T total[2]={0};
                pts[12+2*a]=TV();
                pts[12+2*a+1]=TV();
                for(int v=0;v<8;v++){
                    total[(v>>a)&1]+=1/phi(v);
                    pts[12+2*a+((v>>a)&1)]+=TV(bits(v))/phi(v);}
                pts[12+2*a]/=total[0];
                pts[12+2*a+1]/=total[1];}

            T total=0;
            for(int v=0;v<8;v++){
                total+=1/phi(v);
                pts[18]+=TV(bits(v))/phi(v);}
            pts[18]/=total;}}

    int pat;
    do{
        pat=interface_triangle_table(tri++);
        TRIANGLE_3D<T> triangle(pts[(pat>>10)&0x1f],pts[pat&0x1f],pts[(pat>>5)&0x1f]);
        TRIPLE<TRIANGLE_3D<T>,int,int> triple(triangle,color_list[(pat>>18)&0x7],color_list[(pat>>15)&0x7]);
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
    const int* colors,const T* phi,int s,const HASHTABLE<int,int>& color_map,const int* color_list)
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
        TRIANGLE_3D<T> triangle(pts[(pat>>10)&0x1f],pts[(pat>>5)&0x1f],pts[pat&0x1f]);
        PAIR<TRIANGLE_3D<T>,int> pair(triangle,color_list[(pat>>18)&0x7]);
        if(s) exchange(pair.x.X.y,pair.x.X.z);
        boundary.Append(pair);
    } while(!(pat&last_tri_bit));
}
//#####################################################################
// Function Get_Interface_Elements_For_Cell
//#####################################################################
template<class T> void
Get_Interface_Elements_For_Cell(ARRAY<TRIPLE<SEGMENT_2D<T>,int,int> >& surface,const VECTOR<int,4>& re_color,
    const VECTOR<int,4>& colors,const VECTOR<T,4>& phi,const HASHTABLE<int,int>& color_map,const int* color_list)
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
        SEGMENT_2D<T> segment(pts[(pat>>10)&0x1f],pts[(pat>>5)&0x1f]);
        TRIPLE<SEGMENT_2D<T>,int,int> triple(segment,color_list[(pat>>18)&0x7],color_list[(pat>>15)&0x7]);
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
    const int* colors,const T* phi,int s,const HASHTABLE<int,int>& color_map,const int* color_list)
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
    if(TV::m==3) Initialize_Case_Table_3D();
    if(TV::m==2) Initialize_Case_Table_2D();
}
//#####################################################################
// Function Get_Elements_For_Cell
//#####################################################################
template<class TV> void MARCHING_CUBES_COLOR<TV>::
Get_Elements_For_Cell(ARRAY<TRIPLE<T_FACE,int,int> >& surface,ARRAY<PAIR<T_FACE,int> >& boundary,
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

    Get_Interface_Elements_For_Cell(surface,re_color,colors,phi,color_map,color_list);
    Get_Boundary_Elements_For_Cell(boundary,re_color.array,colors.array,phi.array,0,color_map,color_list);

    next_color=0;
    color_map.Remove_All();
    int re_color2[num_corners/2];
    for(int i=0;i<num_corners/2;i++)
        if(!color_map.Get(colors(i+num_corners/2),re_color2[i])){
            re_color2[i]=next_color;
            color_list[next_color]=colors(i+num_corners/2);
            color_map.Set(colors(i+num_corners/2),next_color++);}
    Get_Boundary_Elements_For_Cell(boundary,re_color2,colors.array+num_corners/2,phi.array+num_corners/2,1,color_map,color_list);
}
template class MARCHING_CUBES_COLOR<VECTOR<float,2> >;
template class MARCHING_CUBES_COLOR<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MARCHING_CUBES_COLOR<VECTOR<double,2> >;
template class MARCHING_CUBES_COLOR<VECTOR<double,3> >;
#endif

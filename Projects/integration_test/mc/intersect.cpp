//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/SEGMENT_3D_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES_COLOR.h>
using namespace PhysBAM;

int Edges_Intersect(int a,int b,int c,int d);
void Initialize_Edges_Intersect();

typedef double T;
typedef float RW;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,3> IV;

HASHTABLE<VECTOR<IV,4> > hash;

IV cy(IV a){return IV(a.y,a.z,a.x);}
IV sw(IV a){return IV(a.y,a.x,a.z);}
IV fl(IV a){return IV(2-a.x,a.y,a.z);}

void Add(const VECTOR<IV,4>& z)
{
    if(hash.Contains(z)){return;}
    hash.Insert(z);
    Add(VECTOR<IV,4>(z(2),z(3),z(0),z(1)));
    Add(VECTOR<IV,4>(z(0),z(1),z(3),z(2)));
    Add(VECTOR<IV,4>(fl(z(0)),fl(z(1)),fl(z(2)),fl(z(3))));
    Add(VECTOR<IV,4>(cy(z(0)),cy(z(1)),cy(z(2)),cy(z(3))));
    Add(VECTOR<IV,4>(sw(z(0)),sw(z(1)),sw(z(2)),sw(z(3))));
}

int main(int argc, char* argv[])
{
    Initialize_Edges_Intersect();
    TV all_pts[256][19];
    for(int i=0;i<256;i++){
        TV* pts=all_pts[i];
        TV cor[8];
        for(int j=0;j<8;j++)
            cor[j]=TV(j&1,j/2&1,j/4&1);
        T phi[8];
        for(int j=0;j<8;j++)
            phi[j]=(i&(1<<j))?10:1;

        for(int a=0,k=0;a<3;a++){
            int mask=1<<a;
            for(int v=0;v<8;v++)
                if(!(v&mask)){
                    T theta=phi[v]/(phi[v]+phi[v|mask]);
                    pts[k++]=(1-theta)*TV(cor[v])+theta*TV(cor[v|mask]);}}

        for(int a=0;a<3;a++){
            T total[2]={0};
            pts[12+2*a]=TV();
            pts[12+2*a+1]=TV();
            for(int v=0;v<8;v++){
                total[(v>>a)&1]+=1/phi[v];
                pts[12+2*a+((v>>a)&1)]+=TV(cor[v])/phi[v];}
            pts[12+2*a]/=total[0];
            pts[12+2*a+1]/=total[1];}

        T total=0;
        for(int v=0;v<8;v++){
            total+=1/phi[v];
            pts[18]+=TV(cor[v])/phi[v];}
        pts[18]/=total;}

    PLANE<T> plane[9];
    plane[0].Specify_Three_Points(TV(0,0,0),TV(0,0,1),TV(1,1,1));
    plane[1].Specify_Three_Points(TV(0,0,0),TV(0,1,0),TV(1,1,1));
    plane[2].Specify_Three_Points(TV(0,0,0),TV(1,0,0),TV(1,1,1));
    plane[3].Specify_Three_Points(TV(1,0,0),TV(1,0,1),TV(0,1,1));
    plane[4].Specify_Three_Points(TV(0,0,1),TV(0,1,1),TV(1,1,0));
    plane[5].Specify_Three_Points(TV(0,1,0),TV(1,1,0),TV(1,0,1));
    plane[6].Specify_Three_Points(TV((T).5,0,0),TV((T).5,0,1),TV((T).5,1,1));
    plane[7].Specify_Three_Points(TV(0,(T).5,0),TV(1,(T).5,0),TV(1,(T).5,1));
    plane[8].Specify_Three_Points(TV(0,0,(T).5),TV(0,1,(T).5),TV(1,1,(T).5));
    int skip_plane[9]={(1<<16)|(1<<17)|(1<<18), (1<<14)|(1<<15)|(1<<18), (1<<12)|(1<<13)|(1<<18), (1<<16)|(1<<17)|(1<<18), (1<<14)|(1<<15)|(1<<18), (1<<12)|(1<<13)|(1<<18)};

    for(int i=0;i<19;i++){
        for(int j=i+1;j<19;j++){
            for(int k=i+1;k<19;k++){
                for(int m=k+1;m<19;m++){
                    if(j!=k && j!=m){
                        TV I(all_pts[0][i]),J(all_pts[0][j]),K(all_pts[0][k]),M(all_pts[0][m]);
                        if(abs((I+J)-(T)1).Max()>.99 || abs((K+M)-(T)1).Max()>.99){
                            printf("OK%i  (%i %i) (%i %i) - face edge\n", Edges_Intersect(i,j,k,m), i, j, k, m);
                            continue;}
                        int status=0; // 1=ok, 2=bad
                        for(int a=0;a<6;a++){
                            if(((1<<i)|(1<<j)|(1<<k)|(1<<m))&skip_plane[a]) continue;
                            T ii=plane[a].Signed_Distance(I),jj=plane[a].Signed_Distance(J),kk=plane[a].Signed_Distance(K),mm=plane[a].Signed_Distance(M);
                            if((ii<1e-10 && jj<1e-10 && kk>-1e-10 && mm>-1e-10) || (ii>-1e-10 && jj>-1e-10 && kk<1e-10 && mm<1e-10)){
                                printf("OK%i  (%i %i) (%i %i) - plane %i  %g %g %g %g\n", Edges_Intersect(i,j,k,m), i, j, k, m, a, ii, jj, kk, mm);
                                status=1;break;}}
                        if(status) continue;

                        for(int a=6;a<9;a++){
                            int ii[4]={i, j, k, m};
                            T s[4]={plane[a].Signed_Distance(I),plane[a].Signed_Distance(J),plane[a].Signed_Distance(K),plane[a].Signed_Distance(M)};
                            int t[4]={0},tt[2]={0};
                            for(int b=0;b<4;b++){
                                if(s[b]>1e-10) t[b]=1;
                                else if(s[b]<-1e-10) t[b]=2;
                                else if(ii[b]<12) t[b]=3;
                                tt[b/2]|=1<<t[b];}
                            if(tt[1]<tt[0]) exchange(tt[1],tt[0]);
                            if(tt[1]&8) continue;
                            if(tt[0]==6 || tt[1]==6) continue;
                            if(tt[0]>1 && tt[1]==tt[0]) continue;
                            status=a;
                            break;}

                        if(status>=6) printf("OK%i  (%i %i) (%i %i) - surface %i\n", Edges_Intersect(i,j,k,m), i, j, k, m, status-6);
                        
                        if(status) continue;

                        bool has[2]={0,0};
                        int stored_case=-1;
                        for(int a=0;a<256;a++){
                            TV I(all_pts[a][i]),J(all_pts[a][j]),K(all_pts[a][k]),M(all_pts[a][m]);
                            T v=TV::Dot_Product(TV::Cross_Product(J-I,K-I),M-I);
                            if(v>1e-6*0) has[1]=1;
                            if(v<-1e-6*0) has[0]=1;
                            if(has[0] && has[1]){
                                status=2;
                                printf("BAD%i  (%i %i) (%i %i) - crossed sign %i %i\n", Edges_Intersect(i,j,k,m), i, j, k, m, a, stored_case);
                                break;}
                            stored_case=a;}
                        if(status) continue;

                        VECTOR<IV,4> KY(IV(rint(I*2)),IV(rint(J*2)),IV(rint(K*2)),IV(rint(M*2)));
                        if(hash.Contains(KY)) continue;
                        Add(KY);

                        printf("UNKNOWN%i (%i %i) (%i %i)\n", Edges_Intersect(i,j,k,m), i, j, k, m);
                    }
                }
            }
        }
    }

    return 0;
}

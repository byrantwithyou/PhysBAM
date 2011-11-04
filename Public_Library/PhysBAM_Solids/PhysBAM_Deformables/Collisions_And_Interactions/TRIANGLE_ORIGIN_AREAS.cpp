#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_ORIGIN_AREAS.h>
#include <PhysBAM_Dynamics/Read_Write/EPS_FILE_GEOMETRY.h>

//extern PhysBAM::ARRAY<int> trap_cases;
namespace PhysBAM{
namespace TRIANGLE_ORIGIN_AREAS{

template<class T,int n> void Clear(VOL_DATA<T,n>& data)
{
    data.V=T();
    for(int i=0;i<n;i++) data.G[i]=VECTOR<T,3>();
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) data.H[i][j]=MATRIX<T,3>();
}

template<class T,class TV> void Data_From_Dof(PT_DATA<T>& data,const TV& A)
{
    data.V=A;
    data.G[0]=MATRIX<T,3>::Identity_Matrix();
    for(int i=0;i<2;i++) data.H[i][0][0]=MATRIX<T,3>();
}

template<class TV> POINT_CASE Classify_Point(const TV& A,const TV& B,const TV& C,const TV& P)
{
    if(TV::Triple_Product(A,B,P)<=0) return outside;
    if(TV::Triple_Product(A,P,C)<=0) return outside;
    if(TV::Triple_Product(P,B,C)<=0) return outside;
    if(TV::Triple_Product(A-P,B-P,C-P)<=0) return beyond;
    return inside;
}

template<class T,class TV> void Volume_From_Points(VOL_DATA<T,3>& data,const TV& A,const TV& B,const TV& C)
{
    data.V.x=(T)(1./6)*TV::Triple_Product(A,B,C);
    // TODO: Gradient and Hessian
}

template<class T,class TV> void Intersect_Triangle_Point(PT_DATA<T>& data,const TV& A,const TV& B,const TV& C,const TV& P)
{
    TV n=TV::Cross_Product(B-A,C-A);
    data.V=TV::Dot_Product(n,A)/TV::Dot_Product(n,P)*P;
    // TODO: Gradient and Hessian
}

template<class T,class TV> void Intersect_Triangle_Segment(PT_DATA<T>& data,const TV& A,const TV& B,const TV& C,const TV& P,const TV& Q)
{
    PT_DATA<T> tdata;
    Intersect_Segment_Point(tdata,A-Q,B-Q,C-Q,P-Q);
    data.V=tdata.V+Q;
    data.G[4]=MATRIX<T,3>::Identity_Matrix();
    for(int i=0;i<4;i++){data.G[i]=tdata.G[i];data.G[4]-=tdata.G[i];}

    for(int s=0;s<3;s++) for(int i=0;i<5;i++){data.H[s][i][4]=data.H[s][4][i]=MATRIX<T,3>();}
    for(int s=0;s<3;s++) for(int i=0;i<4;i++) for(int j=0;j<4;j++){
        data.H[s][i][j]=tdata.H[s][i][j];
        data.H[s][i][4]-=tdata.H[s][i][j];
        data.H[s][4][j]-=tdata.H[s][i][j];
        data.H[s][4][4]+=tdata.H[s][i][j];}
}

template<class T,class TV> void Intersect_Segment_Segment(PT_DATA<T>& data,const TV& A,const TV& B,const TV& P,const TV& Q)
{
    PT_DATA<T> tdata;
    Intersect_Segment_Point(tdata,-Q,A-Q,B-Q,P-Q);
    data.V=tdata.V+Q;
    data.G[4]=MATRIX<T,3>::Identity_Matrix();
    for(int i=0;i<4;i++){
        if(i) data.G[i-1]=tdata.G[i];
        data.G[3]-=tdata.G[i];}

    for(int s=0;s<3;s++) for(int i=0;i<4;i++){data.H[s][i][3]=data.H[s][3][i]=MATRIX<T,3>();}
    for(int s=0;s<3;s++) for(int i=0;i<4;i++) for(int j=0;j<4;j++){
        if(i && j) data.H[s][i-1][j-1]=tdata.H[s][i][j];
        if(i) data.H[s][i-1][3]-=tdata.H[s][i][j];
        if(j) data.H[s][3][j-1]-=tdata.H[s][i][j];
        data.H[s][3][3]+=tdata.H[s][i][j];}
}

// template<class T,class TV> void Volume_From_Points(PT_DATA<T>& data,const TV& A,const TV& B,const TV& C)
// {
//     int index[6]={0,1,2,3,4,5};
//     T sign=1;
//     T OABC=TV::Triple_Product(A,B,C);
//     T ODEF=TV::Triple_Product(D,E,F);

//     if(OABC<0){exchange(A,B);exchange(index[0],index[1]);sign=-sign;}
//     if(ODEF<0){exchange(D,E);exchange(index[3],index[4]);sign=-sign;}

//     POINT_CASE case_a=Classify_Point(D,E,F,A);
//     POINT_CASE case_b=Classify_Point(D,E,F,B);
//     POINT_CASE case_c=Classify_Point(D,E,F,C);
//     POINT_CASE case_d=Classify_Point(A,B,C,D);
//     POINT_CASE case_e=Classify_Point(A,B,C,E);
//     POINT_CASE case_f=Classify_Point(A,B,C,F);
//     if(case_b<case_a){exchange(case_b,case_a);exchange(A,B);exchange(index[0],index[1]);sign=-sign;}
//     if(case_c<case_a){exchange(case_c,case_a);exchange(A,C);exchange(index[0],index[2]);sign=-sign;}
//     if(case_c<case_b){exchange(case_c,case_b);exchange(B,C);exchange(index[1],index[2]);sign=-sign;}
//     if(case_e<case_d){exchange(case_e,case_d);exchange(D,E);exchange(index[3],index[4]);sign=-sign;}
//     if(case_f<case_d){exchange(case_f,case_d);exchange(D,F);exchange(index[3],index[5]);sign=-sign;}
//     if(case_f<case_e){exchange(case_f,case_e);exchange(E,F);exchange(index[4],index[5]);sign=-sign;}
//     if(case_a<case_d || (case_a==case_d && (case_b<case_e || (case_b==case_e && case_c<case_f)))){
//         exchange(case_a,case_d);exchange(case_b,case_e);exchange(case_c,case_f);
//         exchange(A,D);exchange(B,E);exchange(C,F);
//         exchange(index[0],index[3]);exchange(index[1],index[4]);exchange(index[2],index[5]);}

//     DATA<T,1,6> tdata;
//     Clear(tdata);

//     if(case_a==outside && case_b==outside && case_c==outside && case_d==outside && case_e==outside && case_f==outside) return;
//     if(case_a==outside && case_b==outside && case_c==outside && case_d==beyond && case_e==beyond && case_f==beyond) Case_CCCBBB(tdata,A,B,C,D,E,F);
//     if(case_a==outside && case_b==outside && case_c==outside && case_d==inside && case_e==inside && case_f==inside) Case_CCCAAA(tdata,A,B,C,D,E,F);
//     else printf("X: %c %c %c  %c %c %c\n", "ABC"[case_a], "ABC"[case_b], "ABC"[case_c], "ABC"[case_d], "ABC"[case_e], "ABC"[case_f]);

//     // TODO: Enumerate cases

//     data.V=sign*tdata.V;
//     for(int i=0;i<6;i++) data.G[index[i]]=sign*tdata.G[i];
//     for(int i=0;i<6;i++) for(int k=0;k<6;k++) data.H[0][index[i]][index[k]]=sign*tdata.H[0][i][k];
// }

template<class T> void Combine_Data(VOL_DATA<T,6>& data,const VOL_DATA<T,3>& V,const PT_DATA<T>& data_m,const PT_DATA<T>& data_n,const PT_DATA<T>& data_p)
{
    data.V+=V.V;
    PT_DATA<T> pd[3] = {&data_m, &data_n, &data_p};

    for(int z=0;z<3;z++) for(int j=0;j<pd[z].n;j++) data.G[pd[z].index[j]]+=V.G[0]*pd[z].G[j];

    for(int z=0;z<3;z++)
        for(int j=0;j<pd[z].n;j++)
            for(int s=0;s<pd[z].n;s++)
                for(int i=0;i<3;i++)
                    data.H[pd[z].index[j]][pd[z].index[s]]+=V.G[z](1,i+1)*pd[z].H[i][j][s];

    for(int y=0;y<3;y++)
        for(int z=0;z<3;z++)
            for(int j=0;j<pd[y].n;j++)
                for(int s=0;s<pd[z].n;s++)
                    data.H[pd[y].index[j]][pd[z].index[s]]+=pd[y].G[j].Transpose_Times(V.H[y][z]*pd[z].G[s]);
}

// const int vec_d[1]={0}, vec_e[1]={0}, vec_f[1]={0}, vec_abcd[4]={0,1,2,3}, vec_abce[4]={0,1,2,4}, vec_abcf[4]={0,1,2,5};
// template<class T,class TV> void Case_CCCAAA(DATA<T,1,6>& data,const TV& A,const TV& B,const TV& C,const TV& D,const TV& E,const TV& F)
// {
//     LOG::cout<<__FUNCTION__<<std::endl;
//     DATA<T,3,1> DD,DE,DF;
//     Data_From_Dof(DD,D);
//     Data_From_Dof(DE,E);
//     Data_From_Dof(DF,F);
//     DATA<T,1,3> V;
//     Volume_From_Points(V,DD.V,DE.V,DF.V);
//     Combine_Data(data,V,DD,DE,DF,vec_d,vec_e,vec_f);
// }

// template<class T,class TV> void Case_CCCBBB(DATA<T,1,6>& data,const TV& A,const TV& B,const TV& C,const TV& D,const TV& E,const TV& F)
// {
//     LOG::cout<<__FUNCTION__<<std::endl;
//     DATA<T,3,4> P1,P2,P3;
//     Intersect_Triangle_Point(P1,A,B,C,D);
//     Intersect_Triangle_Point(P2,A,B,C,E);
//     Intersect_Triangle_Point(P3,A,B,C,F);
//     DATA<T,1,3> V;
//     Volume_From_Points(V,P1.V,P2.V,P3.V);
//     Combine_Data(data,V,P1,P2,P3,vec_abcd,vec_abce,vec_abcf);
// }

int opp_pairs[5][2] = {{0,0},{1,2},{2,0},{0,0},{0,1}};
int opp_pt[7] = {0,0,0,2,0,1,0};

template<class TV> TV Point_From_Planes(int planes,TV pts[6])
{
    int a=planes&7,b=planes/8;
    if(a){
        TV n=TV::Cross_Product(pts[opp_pairs[b][0]+3],pts[opp_pairs[b][1]+3]);
        TV Q=pts[opp_pairs[a][1]],D=pts[opp_pairs[a][0]]-Q;
        return Q-TV::Dot_Product(n,Q)/TV::Dot_Product(n,D)*D;}
    TV n=TV::Cross_Product(pts[1]-pts[0],pts[2]-pts[0]);
    return TV::Dot_Product(n,pts[0])/TV::Dot_Product(n,pts[opp_pt[b]])*pts[opp_pt[b]];
}

template<class T,class TV> void Data_From_Dof_Helper(PT_DATA<T>& data, int* i, const TV* pts)
{
    data.index[0]=i[0];
    Data_From_Dof(data, pts[i[0]]);
}

template<class T,class TV> void Intersect_Triangle_Point_Helper(PT_DATA<T>& data, int* i, const TV* pts)
{
    for(int j=0;j<4;j++) data.index[j]=i[j];
    Intersect_Triangle_Point(data, pts[i[0]], pts[i[1]], pts[i[2]], pts[i[3]]);
}

template<class T,class TV> void Intersect_Triangle_Segment_Helper(PT_DATA<T>& data, int* i, const TV* pts)
{
    for(int j=0;j<5;j++) data.index[j]=i[j];
    Intersect_Triangle_Segment(data, pts[i[0]], pts[i[1]], pts[i[2]], pts[i[3]], pts[i[4]]);
}

template<class T,class TV> void Intersect_Segment_Segment_Helper(PT_DATA<T>& data, int* i, const TV* pts)
{
    for(int j=0;j<4;j++) data.index[j]=i[j];
    Intersect_Segment_Segment(data, pts[i[0]], pts[i[1]], pts[i[2]], pts[i[3]]);
}

template<class T,class TV> int Init_Data_From_Planes(void (*funcs[256])(PT_DATA<T>& data, int* ind, const TV* pts), int indices[256][5])
{
    for(int i=0;i<3;i++){
        int k=64+7-(1<<i);
        funcs[k]=&Data_From_Dof_Helper;
        indices[k][0]=i;
    }
    return 0;
}

template<class T,class TV> void Volume_From_Tetrahedron(VOL_DATA<T,6>& data,TV pts[6],int va,int vb,int vc)
{
    static void (*funcs[256])(PT_DATA<T>& data, int* ind, const TV* pts);
    static int indices[256][5];
    static int filled = Init_Data_From_Planes(funcs,indices);
    (void) filled;
    // TODO
}

template<class T,class TV> void Volume_From_Triangles_Cut(VOL_DATA<T,6>& data,TV pts[6])
{
    struct LIST
    {
        int planes;
        TV pt;
        bool inside;
    };
    int n=3;
    LIST alist[50], *list=alist;
    for(int i=0;i<n;i++){list[i].planes=7&~(1<<i);list[i].pt=pts[i];}

    // TODO: Robustness to in out in out.
    for(int p=3;p<6;p++){
        TV N=TV::Cross_Product(pts[(p+1)%3+3],pts[(p+2)%3+3]);
        bool has[2]={false,false};
        for(int i=0;i<n;i++){
            list[i].inside=TV::Dot_Product(N,list[i].pt)>0;
            has[list[i].inside]=true;}
        if(!has[1]) return;
        if(!has[0]) continue;
        while(list[n-1].inside || !list[0].inside){list[n]=list[0];list++;}
        int f=1;
        while(list[f].inside) f++;
        int verta=(list[n-1].planes&list[0].planes)|(1<<p),vertb=(list[f-1].planes&list[f].planes)|(1<<p);
        list[f].planes=verta;
        list[f+1].planes=vertb;
        list[f].pt=Point_From_Planes(list[f].planes,pts);
        list[f+1].pt=Point_From_Planes(list[f+1].planes,pts);
        n=f+2;}

    TV N=TV::Cross_Product(pts[4]-pts[3],pts[5]-pts[3]);

    bool has[2]={false,false};
    for(int i=0;i<n;i++){
        list[i].inside=TV::Dot_Product(N,list[i].pt-pts[3])>0;
        has[list[i].inside]=true;}

    if(!has[0]){
        for(int f=2;f<n;f++) Volume_From_Tetrahedron(data,pts,list[0].planes|64,list[f-1].planes|64,list[f].planes|64);
        return;}
    if(!has[1]){
        for(int f=2;f<n;f++) Volume_From_Tetrahedron(data,pts,list[0].planes|128,list[f-1].planes|128,list[f].planes|128);
        return;}

    while(list[n-1].inside || !list[0].inside){list[n]=list[0];list++;}
    int f=1;
    while(list[f].inside) f++;
    int verta=(list[n-1].planes&list[0].planes)|192,vertb=(list[f-1].planes&list[f].planes)|192;

    for(int i=1;i<f;i++) Volume_From_Tetrahedron(data,pts,verta|64,list[i-1].planes|64,list[i].planes|64);
    Volume_From_Tetrahedron(data,pts,verta|64,list[f-1].planes|64,vertb|64);

    for(int i=f+1;i<n;i++) Volume_From_Tetrahedron(data,pts,verta|128,list[i-1].planes|128,list[i].planes|128);
    Volume_From_Tetrahedron(data,pts,verta|128,list[n-1].planes|128,vertb|128);
}

//template void Volume_From_Triangles<float,VECTOR<float,3> >(DATA<float,1,6>&,VECTOR<float,3>,VECTOR<float,3>,VECTOR<float,3>,VECTOR<float,3>,VECTOR<float,3>,VECTOR<float,3>);
template void Volume_From_Triangles_Cut<float,VECTOR<float,3> >(VOL_DATA<float,6>& data,VECTOR<float,3> pts[6]);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
//template void Volume_From_Triangles<double,VECTOR<double,3> >(DATA<double,1,6>&,VECTOR<double,3>,VECTOR<double,3>,VECTOR<double,3>,VECTOR<double,3>,VECTOR<double,3>,VECTOR<double,3>);
#endif
}
}

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_ORIGIN_AREAS.h>
#include <PhysBAM_Dynamics/Read_Write/EPS_FILE_GEOMETRY.h>

//extern PhysBAM::ARRAY<int> trap_cases;
namespace PhysBAM{
namespace TRIANGLE_ORIGIN_AREAS{

template<class T,int m,int n> void Clear(DATA<T,m,n>& data)
{
    data.V=VECTOR<T,m>();
    for(int i=0;i<n;i++) data.G[i]=MATRIX<T,m,3>();
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) for(int k=0;k<m;k++) data.H[k][i][j]=MATRIX<T,3>();
}

template<class T,class TV> void Data_From_Dof(DATA<T,3,1>& data,const TV& A)
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

template<class T,class TV> void Volume_From_Points(DATA<T,1,3>& data,const TV& A,const TV& B,const TV& C)
{
    data.V.x=(T)(1./6)*TV::Triple_Product(A,B,C);
    // TODO: Gradient and Hessian
}


template<class T,class TV> void Intersect_Triangle_Point(DATA<T,3,4>& data,const TV& A,const TV& B,const TV& C,const TV& P)
{
    TV n=TV::Cross_Product(B-A,C-A);
    data.V=TV::Dot_Product(n,A)/TV::Dot_Product(n,P)*P;
    // TODO: Gradient and Hessian
}

template<class T,class TV> void Intersect_Triangle_Segment(DATA<T,3,5>& data,const TV& A,const TV& B,const TV& C,const TV& P,const TV& Q)
{
    DATA<T,3,4> tdata;
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

template<class T,class TV> void Intersect_Segment_Segment(DATA<T,3,4>& data,const TV& A,const TV& B,const TV& P,const TV& Q)
{
    DATA<T,3,4> tdata;
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

template<class T,class TV> void Volume_From_Triangles(DATA<T,1,6>& data,TV A,TV B,TV C,TV D,TV E,TV F)
{
    int index[6]={0,1,2,3,4,5};
    T sign=1;
    T OABC=TV::Triple_Product(A,B,C);
    T ODEF=TV::Triple_Product(D,E,F);

    if(OABC<0){exchange(A,B);exchange(index[0],index[1]);sign=-sign;}
    if(ODEF<0){exchange(D,E);exchange(index[3],index[4]);sign=-sign;}

    POINT_CASE case_a=Classify_Point(D,E,F,A);
    POINT_CASE case_b=Classify_Point(D,E,F,B);
    POINT_CASE case_c=Classify_Point(D,E,F,C);
    POINT_CASE case_d=Classify_Point(A,B,C,D);
    POINT_CASE case_e=Classify_Point(A,B,C,E);
    POINT_CASE case_f=Classify_Point(A,B,C,F);
    if(case_b<case_a){exchange(case_b,case_a);exchange(A,B);exchange(index[0],index[1]);sign=-sign;}
    if(case_c<case_a){exchange(case_c,case_a);exchange(A,C);exchange(index[0],index[2]);sign=-sign;}
    if(case_c<case_b){exchange(case_c,case_b);exchange(B,C);exchange(index[1],index[2]);sign=-sign;}
    if(case_e<case_d){exchange(case_e,case_d);exchange(D,E);exchange(index[3],index[4]);sign=-sign;}
    if(case_f<case_d){exchange(case_f,case_d);exchange(D,F);exchange(index[3],index[5]);sign=-sign;}
    if(case_f<case_e){exchange(case_f,case_e);exchange(E,F);exchange(index[4],index[5]);sign=-sign;}
    if(case_a<case_d || (case_a==case_d && (case_b<case_e || (case_b==case_e && case_c<case_f)))){
        exchange(case_a,case_d);exchange(case_b,case_e);exchange(case_c,case_f);
        exchange(A,D);exchange(B,E);exchange(C,F);
        exchange(index[0],index[3]);exchange(index[1],index[4]);exchange(index[2],index[5]);}

    DATA<T,1,6> tdata;
    Clear(tdata);

    if(case_a==outside && case_b==outside && case_c==outside && case_d==beyond && case_e==beyond && case_f==beyond) Case_CCCBBB(tdata,A,B,C,D,E,F);
    else printf("X: %c %c %c  %c %c %c\n", "ABC"[case_a], "ABC"[case_b], "ABC"[case_c], "ABC"[case_d], "ABC"[case_e], "ABC"[case_f]);

    // TODO: Enumerate cases

    data.V=sign*tdata.V;
    for(int i=0;i<6;i++) data.G[index[i]]=sign*tdata.G[i];
    for(int i=0;i<6;i++) for(int k=0;k<6;k++) data.H[0][index[i]][index[k]]=sign*tdata.H[0][i][k];
}

template<class T,int m,int n,int p> void Combine_Data(DATA<T,1,6>& data,const DATA<T,1,3>& V,const DATA<T,3,m>& data_m,const DATA<T,3,n>& data_n,const DATA<T,3,p>& data_p,
    const int index_m[m],const int index_n[n],const int index_p[p])
{
    data.V+=V.V;
    for(int j=0;j<m;j++) data.G[index_m[j]]+=V.G[0]*data_m.G[j];
    for(int j=0;j<n;j++) data.G[index_n[j]]+=V.G[1]*data_n.G[j];

    for(int j=0;j<m;j++)
        for(int s=0;s<m;s++)
            for(int i=0;i<3;i++)
                data.H[0][index_m[j]][index_m[s]]+=V.G[0](1,i+1)*data_m.H[i][j][s];

    for(int j=0;j<n;j++)
        for(int s=0;s<n;s++)
            for(int i=0;i<3;i++)
                data.H[0][index_n[j]][index_n[s]]+=V.G[1](1,i+1)*data_n.H[i][j][s];

    for(int j=0;j<p;j++)
        for(int s=0;s<p;s++)
            for(int i=0;i<3;i++)
                data.H[0][index_p[j]][index_p[s]]+=V.G[1](1,i+1)*data_p.H[i][j][s];

    for(int j=0;j<m;j++)
        for(int s=0;s<m;s++)
            data.H[0][index_m[j]][index_m[s]]+=data_m.G[j].Transpose_Times(V.H[0][0][0]*data_m.G[s]);

    for(int j=0;j<n;j++)
        for(int s=0;s<n;s++)
            data.H[0][index_n[j]][index_n[s]]+=data_n.G[j].Transpose_Times(V.H[0][1][1]*data_n.G[s]);

    for(int j=0;j<p;j++)
        for(int s=0;s<p;s++)
            data.H[0][index_p[j]][index_p[s]]+=data_p.G[j].Transpose_Times(V.H[0][1][1]*data_p.G[s]);

    for(int j=0;j<m;j++)
        for(int s=0;s<n;s++){
            MATRIX<T,3> x=data_m.G[j].Transpose_Times(V.H[0][0][1]*data_n.G[s]);
            data.H[0][index_m[j]][index_n[s]]+=x;
            data.H[0][index_n[s]][index_m[j]]+=x.Transposed();}

    for(int j=0;j<m;j++)
        for(int s=0;s<p;s++){
            MATRIX<T,3> x=data_m.G[j].Transpose_Times(V.H[0][0][1]*data_p.G[s]);
            data.H[0][index_m[j]][index_p[s]]+=x;
            data.H[0][index_p[s]][index_m[j]]+=x.Transposed();}

    for(int j=0;j<p;j++)
        for(int s=0;s<n;s++){
            MATRIX<T,3> x=data_p.G[j].Transpose_Times(V.H[0][0][1]*data_n.G[s]);
            data.H[0][index_p[j]][index_n[s]]+=x;
            data.H[0][index_n[s]][index_p[j]]+=x.Transposed();}
}

const int vec_a[1]={0}, vec_abcd[4]={0,1,2,3}, vec_abce[4]={0,1,2,4}, vec_abcf[4]={0,1,2,5};
template<class T,class TV> void Case_CCCAAA(DATA<T,1,6>& data,const TV& A,const TV& B,const TV& C,const TV& D,const TV& E,const TV& F)
{
    // TODO: Fill in
}

template<class T,class TV> void Case_CCCBBB(DATA<T,1,6>& data,const TV& A,const TV& B,const TV& C,const TV& D,const TV& E,const TV& F)
{
    LOG::cout<<__FUNCTION__<<std::endl;
    DATA<T,3,4> P1,P2,P3;
    Intersect_Triangle_Point(P1,A,B,C,D);
    Intersect_Triangle_Point(P2,A,B,C,E);
    Intersect_Triangle_Point(P3,A,B,C,F);
    DATA<T,1,3> V;
    Volume_From_Points(V,P1.V,P2.V,P3.V);
    Combine_Data(data,V,P1,P2,P3,vec_abcd,vec_abce,vec_abcf);

    static int cnt=0;cnt++;
    for(int i=1;i<=3;i++){
        char file[100];
        sprintf(file, "dump-%c-%i.eps", 'x'+i-1, cnt);
        EPS_FILE_GEOMETRY<T> epsx(file);
        epsx.Line_Color(VECTOR<T,3>(.5,.5,.5));
        epsx.Draw_Point(VECTOR<T,2>(0,0));
        epsx.Draw_Point(VECTOR<T,2>(1,1));
        epsx.Draw_Point(VECTOR<T,2>(-1,-1));
        epsx.Line_Color(VECTOR<T,3>(1,0,0));
        epsx.Draw_Line(A.Remove_Index(i),B.Remove_Index(i));
        epsx.Draw_Line(B.Remove_Index(i),C.Remove_Index(i));
        epsx.Draw_Line(C.Remove_Index(i),A.Remove_Index(i));
        epsx.Line_Color(VECTOR<T,3>(0,1,0));
        epsx.Draw_Line(D.Remove_Index(i),E.Remove_Index(i));
        epsx.Draw_Line(E.Remove_Index(i),F.Remove_Index(i));
        epsx.Draw_Line(F.Remove_Index(i),D.Remove_Index(i));
        
        epsx.Line_Color(VECTOR<T,3>(1,0,0));
        epsx.Draw_Point(A.Remove_Index(i));
        epsx.Draw_Point(D.Remove_Index(i));

        epsx.Line_Color(VECTOR<T,3>(0,1,0));
        epsx.Draw_Point(B.Remove_Index(i));
        epsx.Draw_Point(E.Remove_Index(i));

        epsx.Line_Color(VECTOR<T,3>(0,0,1));
        epsx.Draw_Point(C.Remove_Index(i));
        epsx.Draw_Point(F.Remove_Index(i));

        epsx.Line_Color(VECTOR<T,3>(0,1,1));
        epsx.Draw_Point(P1.V.Remove_Index(i));
        epsx.Draw_Point(P2.V.Remove_Index(i));
        epsx.Draw_Point(P3.V.Remove_Index(i));

        epsx.Draw_Line(D.Remove_Index(i),VECTOR<T,2>());
        epsx.Draw_Line(E.Remove_Index(i),VECTOR<T,2>());
        epsx.Draw_Line(F.Remove_Index(i),VECTOR<T,2>());

}
}

template void Volume_From_Triangles<float,VECTOR<float,3> >(DATA<float,1,6>&,VECTOR<float,3>,VECTOR<float,3>,VECTOR<float,3>,VECTOR<float,3>,VECTOR<float,3>,VECTOR<float,3>);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void Volume_From_Triangles<double,VECTOR<double,3> >(DATA<double,1,6>&,VECTOR<double,3>,VECTOR<double,3>,VECTOR<double,3>,VECTOR<double,3>,VECTOR<double,3>,VECTOR<double,3>);
#endif
}
}

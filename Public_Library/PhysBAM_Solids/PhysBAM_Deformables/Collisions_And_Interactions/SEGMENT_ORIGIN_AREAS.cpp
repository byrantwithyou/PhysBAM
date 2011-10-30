#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/SEGMENT_ORIGIN_AREAS.h>

extern PhysBAM::ARRAY<int> trap_cases;
namespace PhysBAM{
namespace SEGMENT_ORIGIN_AREAS{

template<class T,int m,int n> void Clear(DATA<T,m,n>& data)
{
    data.V=VECTOR<T,m>();
    for(int i=0;i<n;i++) data.G[i]=MATRIX<T,m,2>();
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) for(int k=0;k<m;k++) data.H[k][i][j]=MATRIX<T,2>();
}

template<class T,class TV> void Data_From_Dof(DATA<T,2,1>& data,const TV& A)
{
    data.V=A;
    data.G[0]=MATRIX<T,2>::Identity_Matrix();
    for(int i=0;i<2;i++) data.H[i][0][0]=MATRIX<T,2>();
}

template<class TV> POINT_CASE Classify_Point(const TV& A,const TV& B,const TV& P)
{
    if(TV::Cross_Product(A,P).x<0) return outside;
    if(TV::Cross_Product(P,B).x<0) return outside;
    if(TV::Cross_Product(A-P,B-P).x<0) return beyond;
    return inside;
}

template<class T,class TV> void Intersect_Segment_Point(DATA<T,2,3>& data,const TV& A,const TV& B,const TV& P)
{
    T AxB=TV::Cross_Product(A,B).x,cross2=TV::Cross_Product(A-B,P-B).x,PxB=TV::Cross_Product(P,B).x,PxA=TV::Cross_Product(P,A).x;
    T den=1/(AxB+cross2),sden=sqr(den),cden=den*sden;
    data.V=AxB*den*P;

    TV orthAB=(A-B).Orthogonal_Vector(),orthP=P.Orthogonal_Vector();
    MATRIX<T,2> M=MATRIX<T,2>::Outer_Product(P,sden*orthAB);
    data.G[0]=PxB*M;
    data.G[1]=-PxA*M;
    data.G[2]=-AxB*M+AxB*den;

    MATRIX<T,2> H00=(MATRIX<T,2>::Outer_Product(orthAB,orthP)+MATRIX<T,2>::Outer_Product(orthP,orthAB))*PxB*cden;
    data.H[0][0][0]=H00*P.x;
    data.H[1][0][0]=H00*P.y;
    MATRIX<T,2> H11=(MATRIX<T,2>::Outer_Product(orthAB,orthP)+MATRIX<T,2>::Outer_Product(orthP,orthAB))*PxA*cden;
    data.H[0][1][1]=H11*P.x;
    data.H[1][1][1]=H11*P.y;
    MATRIX<T,2> H22=(MATRIX<T,2>::Outer_Product(orthAB,orthP)+MATRIX<T,2>::Outer_Product(orthP,orthAB))*AxB*cden;
    data.H[0][2][2]=H22*(A.x-B.x);
    data.H[1][2][2]=H22*(A.y-B.y);
    MATRIX<T,2> H01=(-PxA*MATRIX<T,2>::Outer_Product(orthAB,orthP)-PxB*MATRIX<T,2>::Outer_Product(orthP,orthAB))*cden;
    data.H[0][0][1]=H01*P.x;
    data.H[1][0][1]=H01*P.y;
    data.H[0][0][2]=((B.x*(PxA-PxB)-2*(A.x-B.x)*PxB)*MATRIX<T,2>::Outer_Product(orthAB,orthP))*cden;
    data.H[1][0][2]=((B.y*(PxA-PxB)-2*(A.y-B.y)*PxB)*MATRIX<T,2>::Outer_Product(orthAB,orthP))*cden;
    data.H[0][1][2]=((A.x*(PxA-PxB)+2*P.x*AxB)*MATRIX<T,2>::Outer_Product(orthAB,orthP))*cden;
    data.H[1][1][2]=((A.y*(PxA-PxB)+2*P.y*AxB)*MATRIX<T,2>::Outer_Product(orthAB,orthP))*cden;
    data.H[0][1][0]=data.H[0][0][1].Transposed();
    data.H[1][1][0]=data.H[1][0][1].Transposed();
    data.H[0][2][0]=data.H[0][0][2].Transposed();
    data.H[1][2][0]=data.H[1][0][2].Transposed();
    data.H[0][2][1]=data.H[0][1][2].Transposed();
    data.H[1][2][1]=data.H[1][1][2].Transposed();
}

template<class T,class TV> void Intersect_Segments(DATA<T,2,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    DATA<T,2,3> tdata;
    Intersect_Segment_Point(tdata,A-D,B-D,C-D);
    data.V=tdata.V+D;
    data.G[3]=MATRIX<T,2>::Identity_Matrix();
    for(int i=0;i<3;i++){data.G[i]=tdata.G[i];data.G[3]-=tdata.G[i];}

    for(int s=0;s<2;s++) for(int i=0;i<4;i++){data.H[s][i][3]=data.H[s][3][i]=MATRIX<T,2>();}
    for(int s=0;s<2;s++) for(int i=0;i<3;i++) for(int j=0;j<3;j++){data.H[s][i][j]=tdata.H[s][i][j];data.H[s][i][3]-=tdata.H[s][i][j];data.H[s][3][j]-=tdata.H[s][i][j];data.H[s][3][3]+=tdata.H[s][i][j];}
}

template<class T,class TV> void Area_From_Points(DATA<T,1,2>& data,const TV& A,const TV& B)
{
    data.V=(T).5*TV::Cross_Product(A,B);
    data.G[0]=MATRIX<T,1,2>::Cross_Product_Matrix(-(T).5*B);
    data.G[1]=MATRIX<T,1,2>::Cross_Product_Matrix((T).5*A);
    MATRIX<T,2> M(0,-(T).5,(T).5,0);
    data.H[0][0][0]=data.H[0][1][1]=MATRIX<T,2>();
    data.H[0][0][1]=M;
    data.H[0][1][0]=M.Transposed();
}

template<class T,class TV> void Area_From_Segments(DATA<T,1,4>& data,TV A,TV B,TV C,TV D)
{
    int index[4]={0,1,2,3};
    T sign=1;
    T OAB=TV::Cross_Product(A,B).x;
    T OCD=TV::Cross_Product(C,D).x;

    if(OAB<0){exchange(A,B);exchange(index[0],index[1]);sign=-sign;}
    if(OCD<0){exchange(C,D);exchange(index[2],index[3]);sign=-sign;}

    POINT_CASE case_a=Classify_Point(C,D,A);
    POINT_CASE case_b=Classify_Point(C,D,B);
    POINT_CASE case_c=Classify_Point(A,B,C);
    POINT_CASE case_d=Classify_Point(A,B,D);
    if(case_b<case_a){exchange(case_b,case_a);exchange(A,B);exchange(index[0],index[1]);sign=-sign;}
    if(case_d<case_c){exchange(case_c,case_d);exchange(C,D);exchange(index[2],index[3]);sign=-sign;}
    if(case_a<case_c || (case_a==case_c && case_b<case_d)){
        exchange(case_a,case_c);exchange(case_b,case_d);
        exchange(A,C);exchange(B,D);exchange(index[0],index[2]);exchange(index[1],index[3]);}

    DATA<T,1,4> tdata;
    Clear(tdata);
    if(case_a==outside){
        if(case_b==outside && case_c==outside && case_d==outside) return;
        PHYSBAM_ASSERT(case_b==outside && case_c!=outside && case_d!=outside);
        if(case_c==beyond) Case_CCBB(tdata,A,B,C,D);
        else if(case_d==beyond) Case_CCAB(tdata,A,B,C,D);
        else Case_CCAA(tdata,A,B,C,D);}
    else if(case_a==beyond){
        PHYSBAM_ASSERT(case_b==outside && case_d==outside);
        if(case_c==inside){if(TV::Cross_Product(A,B).x<0) sign=-sign;Case_BCAC(tdata,A,B,C,D);}
        else{if(TV::Cross_Product(A,B).x<0) sign=-sign;Case_BCBC(tdata,A,B,C,D);}}
    else{
        PHYSBAM_ASSERT(case_b==outside && case_c==inside && case_d==outside);
        if(TV::Cross_Product(A,B).x<0) sign=-sign;Case_ACAC(tdata,A,B,C,D);}

    for(int i=0;i<4;i++) LOG::cout<<index[i]<<" ";LOG::cout<<std::endl;
    LOG::cout<<"sign "<<sign<<std::endl;

    data.V=sign*tdata.V;
    for(int i=0;i<4;i++) data.G[index[i]]=sign*tdata.G[i];
    for(int i=0;i<4;i++) for(int k=0;k<4;k++) data.H[0][index[i]][index[k]]=sign*tdata.H[0][i][k];
}

template<class T,int m,int n> void Combine_Data(DATA<T,1,4>& data,const DATA<T,1,2>& V,const DATA<T,2,m>& data_m,const DATA<T,2,n>& data_n,const int index_m[m],
    const int index_n[n])
{
    data.V+=V.V;
    for(int j=0;j<m;j++) data.G[index_m[j]]+=V.G[0]*data_m.G[j];
    for(int j=0;j<n;j++) data.G[index_n[j]]+=V.G[1]*data_n.G[j];

    for(int j=0;j<m;j++)
        for(int s=0;s<m;s++)
            for(int i=0;i<2;i++)
                data.H[0][index_m[j]][index_m[s]]+=V.G[0](1,i+1)*data_m.H[i][j][s];

    for(int j=0;j<n;j++)
        for(int s=0;s<n;s++)
            for(int i=0;i<2;i++)
                data.H[0][index_n[j]][index_n[s]]+=V.G[1](1,i+1)*data_n.H[i][j][s];

    for(int j=0;j<m;j++)
        for(int s=0;s<m;s++)
            data.H[0][index_m[j]][index_m[s]]+=data_m.G[j].Transpose_Times(V.H[0][0][0]*data_m.G[s]);

    for(int j=0;j<n;j++)
        for(int s=0;s<n;s++)
            data.H[0][index_n[j]][index_n[s]]+=data_n.G[j].Transpose_Times(V.H[0][1][1]*data_n.G[s]);

    for(int j=0;j<m;j++)
        for(int s=0;s<n;s++){
            MATRIX<T,2> x=data_m.G[j].Transpose_Times(V.H[0][0][1]*data_n.G[s]);
            data.H[0][index_m[j]][index_n[s]]+=x;
            data.H[0][index_n[s]][index_m[j]]+=x.Transposed();}
}

const int vec_a[1]={0}, vec_c[1]={2}, vec_d[1]={3}, vec_abc[3]={0,1,2}, vec_abd[3]={0,1,3}, vec_cda[3]={2,3,0}, vec_cdb[3]={2,3,1}, vec_abcd[4]={0,1,2,3};
template<class T,class TV> void Case_CCAA(DATA<T,1,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    trap_cases.Append(1);
    DATA<T,2,1> DC;
    Data_From_Dof(DC,C);

    DATA<T,2,1> DD;
    Data_From_Dof(DD,D);

    DATA<T,1,2> V;
    Area_From_Points(V,DC.V,DD.V);
    Combine_Data(data,V,DC,DD,vec_c,vec_d);
}

template<class T,class TV> void Case_CCAB(DATA<T,1,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    trap_cases.Append(2);
    DATA<T,2,4> Q;
    Intersect_Segments(Q,A,B,C,D);

    DATA<T,2,3> P;
    Intersect_Segment_Point(P,A,B,D);

    DATA<T,2,1> DC;
    Data_From_Dof(DC,C);

    DATA<T,1,2> V;
    Area_From_Points(V,DC.V,Q.V);
    Combine_Data(data,V,DC,Q,vec_c,vec_abcd);

    Area_From_Points(V,Q.V,P.V);
    Combine_Data(data,V,Q,P,vec_abcd,vec_abd);
}

template<class T,class TV> void Case_CCBB(DATA<T,1,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    trap_cases.Append(3);
    DATA<T,2,3> P1;
    Intersect_Segment_Point(P1,A,B,C);

    DATA<T,2,3> P2;
    Intersect_Segment_Point(P2,A,B,D);

    DATA<T,1,2> V;
    Area_From_Points(V,P1.V,P2.V);
    Combine_Data(data,V,P1,P2,vec_abc,vec_abd);
}

template<class T,class TV> void Case_BCAC(DATA<T,1,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    trap_cases.Append(4);
    DATA<T,2,3> P;
    Intersect_Segment_Point(P,C,D,A);

    DATA<T,2,1> DC;
    Data_From_Dof(DC,C);

    DATA<T,1,2> V;
    Area_From_Points(V,DC.V,P.V);
    Combine_Data(data,V,DC,P,vec_c,vec_cda);
}

template<class T,class TV> void Case_BCBC(DATA<T,1,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    trap_cases.Append(5);
    DATA<T,2,4> Q;
    Intersect_Segments(Q,A,B,C,D);

    DATA<T,2,3> P1;
    Intersect_Segment_Point(P1,C,D,A);

    DATA<T,2,3> P2;
    Intersect_Segment_Point(P2,A,B,C);

    DATA<T,1,2> V;
    Area_From_Points(V,Q.V,P1.V);
    Combine_Data(data,V,Q,P1,vec_abcd,vec_cda);

    Area_From_Points(V,P2.V,Q.V);
    Combine_Data(data,V,P2,Q,vec_abc,vec_abcd);
}

template<class T,class TV> void Case_ACAC(DATA<T,1,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    trap_cases.Append(6);
    DATA<T,2,4> Q;
    Intersect_Segments(Q,A,B,C,D);

    DATA<T,2,1> DA;
    Data_From_Dof(DA,A);

    DATA<T,2,1> DC;
    Data_From_Dof(DC,C);

    DATA<T,1,2> V;
    Area_From_Points(V,DC.V,Q.V);
    Combine_Data(data,V,DC,Q,vec_c,vec_abcd);

    Area_From_Points(V,Q.V,DA.V);
    Combine_Data(data,V,Q,DA,vec_abcd,vec_a);
}

template void Area_From_Segments<float,VECTOR<float,2> >(DATA<float,1,4>&,VECTOR<float,2>,VECTOR<float,2>,VECTOR<float,2>,
    VECTOR<float,2>);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void Area_From_Segments<double,VECTOR<double,2> >(DATA<double,1,4>&,VECTOR<double,2>,VECTOR<double,2>,VECTOR<double,2>,
    VECTOR<double,2>);
#endif
}
}

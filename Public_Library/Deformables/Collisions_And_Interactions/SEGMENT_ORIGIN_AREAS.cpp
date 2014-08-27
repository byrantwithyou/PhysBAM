#include <Tools/Arrays/ARRAY.h>
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Utilities/PROCESS_UTILITIES.h>
#include <Deformables/Collisions_And_Interactions/SEGMENT_ORIGIN_AREAS.h>

//extern PhysBAM::ARRAY<int> trap_cases;
using namespace PhysBAM;
using namespace ORIGIN_AREAS;

template<class T,int m,int n>
struct DATA
{
    VECTOR<T,m> V;
    MATRIX<T,m,2> G[n];
    MATRIX<T,2> H[m][n][n];
};

template<class T,int m,int n> void PhysBAM::ORIGIN_AREAS::Clear(VOL_DATA<T,m,n>& data)
{
    data.V=T();
    for(int i=0;i<n;i++) data.G[i]=VECTOR<T,m>();
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) data.H[i][j]=MATRIX<T,m>();
}

template<class T,class TV> void Triangle_Area(
    VOL_DATA<T,2,3>& data,
    const TV& A,const TV& B,const TV& C)
{
    // Computes the area [A,B,C], along with the 1st and 2nd derivatives.
    const TV AC_L=(A-C).Orthogonal_Vector();
    const TV AB_L=(A-B).Orthogonal_Vector();
    const TV BC_L=(B-C).Orthogonal_Vector();
    typedef MATRIX<T,2> TM;
    const TM J=TM(0,+1,-1,0)/2;
    data.V=TV::Cross_Product(A-C,B-C).x/2;
    data.G[0]=-BC_L/2;
    data.G[1]= AC_L/2;
    data.G[2]=-AB_L/2;
    data.H[0][0]=TM();
    data.H[1][1]=TM();
    data.H[2][2]=TM();
    data.H[0][1]=-J;
    data.H[0][2]= J;
    data.H[1][0]= J;
    data.H[1][2]=-J;
    data.H[2][0]=-J;
    data.H[2][1]= J;
}

template<class T,class TV> void Triangle_Area(
    VOL_DATA<T,2,4>& data,
    const TV& A,const TV& B,const TV& C,const TV& D)
{
    // Computes the area [B,C,P], along with 1st and 2nd derivatives, where P is
    // the intersection between segment A-B and segment C-D.
    const T ABC=TV::Cross_Product(A-C,B-C).x;
    const T ADB=TV::Cross_Product(A-B,D-B).x;
    const T ADC=TV::Cross_Product(A-C,D-C).x;
    const T BCD=TV::Cross_Product(B-D,C-D).x;
    const T ADBC=TV::Cross_Product(A-B,D-C).x;
    const TV AB_L=(A-B).Orthogonal_Vector();
    const TV AC_L=(A-C).Orthogonal_Vector();
    const TV AD_L=(A-D).Orthogonal_Vector();
    const TV BC_L=(B-C).Orthogonal_Vector();
    const TV BD_L=(B-D).Orthogonal_Vector();
    const TV CD_L=(C-D).Orthogonal_Vector();
    typedef MATRIX<T,2> TM;
    const TM J=TM(0,+1,-1,0)/2;
    const TM ABxAB=TM::Outer_Product(AB_L,AB_L);
    const TM ABxAC=TM::Outer_Product(AB_L,AC_L);
    const TM ABxAD=TM::Outer_Product(AB_L,AD_L);
    const TM ABxBC=TM::Outer_Product(AB_L,BC_L);
    const TM ABxBD=TM::Outer_Product(AB_L,BD_L);
    const TM ABxCD=TM::Outer_Product(AB_L,CD_L);
    const TM ACxAB=ABxAC.Transposed();
    const TM ACxBC=TM::Outer_Product(AC_L,BC_L);
    const TM ACxBD=TM::Outer_Product(AC_L,BD_L);
    const TM ACxCD=TM::Outer_Product(AC_L,CD_L);
    const TM ADxAB=ABxAD.Transposed();
    const TM ADxCD=TM::Outer_Product(AD_L,CD_L);
    const TM BCxAB=ABxBC.Transposed();
    const TM BCxBC=TM::Outer_Product(BC_L,BC_L);
    const TM BCxBD=TM::Outer_Product(BC_L,BD_L);
    const TM BDxAB=ABxBD.Transposed();
    const TM BCxCD=TM::Outer_Product(BC_L,CD_L);
    const TM BDxCD=TM::Outer_Product(BD_L,CD_L);
    const TM CDxAB=ABxCD.Transposed();
    const TM CDxAC=ACxCD.Transposed();
    const TM CDxAD=ADxCD.Transposed();
    const TM CDxBC=BCxCD.Transposed();
    const TM CDxBD=BDxCD.Transposed();
    const TM CDxCD=TM::Outer_Product(CD_L,CD_L);
    const T ADBC2=ADBC*ADBC;
    data.V=ABC*BCD/(2*ADBC);
    data.G[0]=-(BCD/(2*ADBC))*(BC_L+(ABC/ADBC)*CD_L);
    data.G[1]=(1/(2*ADBC))*(BCD*AC_L-(ABC*ADC/ADBC)*CD_L);
    data.G[2]=(1/(2*ADBC))*(ABC*BD_L-(ADB*BCD/ADBC)*AB_L);
    data.G[3]=-(ABC/(2*ADBC))*(BC_L+(BCD/ADBC)*AB_L);
    data.H[0][0]=+(BCD/ADBC2)*((ABC/ADBC)*CDxCD+(BCxCD+CDxBC)/2);
    data.H[1][1]=-(ADC/ADBC2)*((ABC/ADBC)*CDxCD+(ACxCD+CDxAC)/2);
    data.H[2][2]=-(ADB/ADBC2)*((BCD/ADBC)*ABxAB+(ABxBD+BDxAB)/2);
    data.H[3][3]=+(ABC/ADBC2)*((BCD/ADBC)*ABxAB+(ABxBC+BCxAB)/2);
    data.H[0][1]=-(BCD/ADBC)*J+(1/(4*ADBC2))*(ADC*(BCxCD+BDxCD)-BCD*(CDxAC+CDxAD)+((ABC-ADB)*(ADC-BCD)/ADBC)*CDxCD);
    data.H[2][3]=-(ABC/ADBC)*J+(1/(4*ADBC2))*(ADB*(ABxAC+ABxBC)-ABC*(ADxAB+BDxAB)+((ABC-ADB)*(ADC-BCD)/ADBC)*ABxAB);
    data.H[1][2]=-(1-ADB*ADC/ADBC2)*J+(1/(4*ADBC2))*(2*ADBC*ACxBD+ABC*(CDxAD+CDxBD)+BCD*(ACxAB+ADxAB)+(4*ADB*ADC/ADBC-(ADB+ADC))*CDxAB);
    data.H[0][2]=(ADB*BCD/ADBC2)*((1/ADBC)*CDxAB+J)-(1/(2*ADBC2))*(ADBC*BCxBD+ABC*CDxBD+BCD*BDxAB);
    data.H[0][3]=(ABC*BCD/ADBC2)*((1/ADBC)*CDxAB+J)+(1/(2*ADBC2))*(ADBC*BCxBC+ABC*CDxBC+BCD*BCxAB);
    data.H[1][3]=(ABC*ADC/ADBC2)*((1/ADBC)*CDxAB+J)-(1/(2*ADBC2))*(ADBC*ACxBC+ABC*CDxAC+BCD*ACxAB);
    for(int i=1;i!=4;++i)
        for(int j=0;j!=i;++j)
            data.H[i][j]=data.H[j][i].Transposed();
}

template<class T,class TV> void Data_From_Dof(DATA<T,2,1>& data,const TV& A)
{
    data.V=A;
    data.G[0]=MATRIX<T,2>::Identity_Matrix();
    for(int i=0;i<2;i++) data.H[i][0][0]=MATRIX<T,2>();
}

template<class T,class TV> void Intersect_Segment_Point(DATA<T,2,3>& data,const TV& A,const TV& B,const TV& P)
{
    T AxB=TV::Cross_Product(A,B).x,PxB=TV::Cross_Product(P,B).x,PxA=TV::Cross_Product(P,A).x,PxBmA=TV::Cross_Product(P,B-A).x;
    T den=(T)1/PxBmA,sden=sqr(den),cden=den*sden;
    data.V=AxB*den*P;

    TV orthAB=(A-B).Orthogonal_Vector(),orthP=P.Orthogonal_Vector();
    MATRIX<T,2> M=MATRIX<T,2>::Outer_Product(P,sden*orthAB);
    data.G[0]=PxB*M;
    data.G[1]=-PxA*M;
    data.G[2]=-AxB*M+AxB*den;

//    if(fabs(PxB-PxA)<1e-5){
//        LOG::cout<<"GRAD: "<<data.G[0]<<"  "<<data.G[1]<<"  "<<data.G[2]<<std::endl;}

    MATRIX<T,2> ABP=MATRIX<T,2>::Outer_Product(orthAB,orthP*cden),PAB=MATRIX<T,2>::Outer_Product(orthP*cden,orthAB),PAB_ABP=ABP+PAB;
    MATRIX<T,2> H00=PAB_ABP*PxB,H11=PAB_ABP*PxA,H22=PAB_ABP*AxB,H01=-PxA*ABP-PxB*PAB;
    data.H[0][0][0]=H00*P.x;
    data.H[1][0][0]=H00*P.y;
    data.H[0][1][1]=H11*P.x;
    data.H[1][1][1]=H11*P.y;
    data.H[0][2][2]=H22*(A.x-B.x);
    data.H[1][2][2]=H22*(A.y-B.y);
    data.H[0][0][1]=H01*P.x;
    data.H[1][0][1]=H01*P.y;
    data.H[0][0][2]=(-B.x*PxBmA-(T)2*(A.x-B.x)*PxB)*ABP;
    data.H[1][0][2]=(-B.y*PxBmA-(T)2*(A.y-B.y)*PxB)*ABP;
    data.H[0][1][2]=(-A.x*PxBmA+(T)2*P.x*AxB)*ABP;
    data.H[1][1][2]=(-A.y*PxBmA+(T)2*P.y*AxB)*ABP;
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

template<class T,class TV> void Area_From_Points(VOL_DATA<T,2,2>& data,const TV& A,const TV& B)
{
    data.V=(T).5*TV::Cross_Product(A,B).x;
    data.G[0]=-(T).5*B.Orthogonal_Vector();
    data.G[1]=(T).5*A.Orthogonal_Vector();
    MATRIX<T,2> M(0,-(T).5,(T).5,0);
    data.H[0][0]=data.H[1][1]=MATRIX<T,2>();
    data.H[0][1]=M;
    data.H[1][0]=M.Transposed();
}

template<class T,class TV> void Area_From_Points(VOL_DATA<T,2,3>& data,const TV& A,const TV& B,const TV& C)
{
    data.V=(T).5*TV::Cross_Product(B,C-A).x;
    data.G[0]=-(T).5*B.Orthogonal_Vector();
    data.G[1]=(T).5*(A-C).Orthogonal_Vector();
    data.G[2]=(T).5*B.Orthogonal_Vector();
    const MATRIX<T,2> M(0,-(T).5,(T).5,0);
    data.H[0][0]=data.H[0][2]=data.H[1][1]=data.H[2][0]=data.H[2][2]=MATRIX<T,2>();
    data.H[0][1]=data.H[1][2]=M;
    data.H[1][0]=data.H[2][1]=M.Transposed();
}

//template<class T,int m,int n> void Combine_Data(VOL_DATA<T,2,4>& data,const VOL_DATA<T,2,2>& V,const DATA<T,2,m>& data_m,const DATA<T,2,n>& data_n,const int index_m[m],
//    const int index_n[n])
template<class T,int m,int n> void Combine_Data(VOL_DATA<T,2,4>& data,const VOL_DATA<T,2,2>& V,const DATA<T,2,m>& data_m,const DATA<T,2,n>& data_n,const int* index_m,
    const int* index_n)
{
    data.V+=V.V;
    for(int j=0;j<m;j++) data.G[index_m[j]]+=data_m.G[j].Transpose_Times(V.G[0]);
    for(int j=0;j<n;j++) data.G[index_n[j]]+=data_n.G[j].Transpose_Times(V.G[1]);

    for(int j=0;j<m;j++)
        for(int s=0;s<m;s++)
            for(int i=0;i<2;i++)
                data.H[index_m[j]][index_m[s]]+=V.G[0](i+1)*data_m.H[i][j][s];

    for(int j=0;j<n;j++)
        for(int s=0;s<n;s++)
            for(int i=0;i<2;i++)
                data.H[index_n[j]][index_n[s]]+=V.G[1](i+1)*data_n.H[i][j][s];

    for(int j=0;j<m;j++)
        for(int s=0;s<m;s++)
            data.H[index_m[j]][index_m[s]]+=data_m.G[j].Transpose_Times(V.H[0][0]*data_m.G[s]);

    for(int j=0;j<n;j++)
        for(int s=0;s<n;s++)
            data.H[index_n[j]][index_n[s]]+=data_n.G[j].Transpose_Times(V.H[1][1]*data_n.G[s]);

    for(int j=0;j<m;j++)
        for(int s=0;s<n;s++){
            MATRIX<T,2> x=data_m.G[j].Transpose_Times(V.H[0][1]*data_n.G[s]);
            data.H[index_m[j]][index_n[s]]+=x;
            data.H[index_n[s]][index_m[j]]+=x.Transposed();}
}

template<class T,int n1,int n2,int n3>
void Combine_Data(
    VOL_DATA<T,2,4>& data,
    const VOL_DATA<T,2,3>& V,
    const DATA<T,2,n1>& data1,
    const DATA<T,2,n2>& data2,
    const DATA<T,2,n3>& data3,
    const int (&index1)[n1],
    const int (&index2)[n2],
    const int (&index3)[n3])
{
    data.V+=V.V;
    for(int i=0;i!=n1;++i) data.G[index1[i]]+=data1.G[i].Transpose_Times(V.G[0]);
    for(int i=0;i!=n2;++i) data.G[index2[i]]+=data2.G[i].Transpose_Times(V.G[1]);
    for(int i=0;i!=n3;++i) data.G[index3[i]]+=data3.G[i].Transpose_Times(V.G[2]);

    for(int i=0;i!=n1;++i)
        for(int j=0;j!=n1;++j)
            for(int d=0;d!=2;++d)
                data.H[index1[i]][index1[j]]+=V.G[0](d+1)*data1.H[d][i][j];
    for(int i=0;i!=n2;++i)
        for(int j=0;j!=n2;++j)
            for(int d=0;d!=2;++d)
                data.H[index2[i]][index2[j]]+=V.G[1](d+1)*data2.H[d][i][j];
    for(int i=0;i!=n3;++i)
        for(int j=0;j!=n3;++j)
            for(int d=0;d!=2;++d)
                data.H[index3[i]][index3[j]]+=V.G[2](d+1)*data3.H[d][i][j];

    for(int i=0;i!=n1;++i)
        for(int j=0;j!=n1;++j)
            data.H[index1[i]][index1[j]]+=
                data1.G[i].Transpose_Times(V.H[0][0]*data1.G[j]);
    for(int i=0;i!=n2;++i)
        for(int j=0;j!=n2;++j)
            data.H[index2[i]][index2[j]]+=
                data2.G[i].Transpose_Times(V.H[1][1]*data2.G[j]);
    for(int i=0;i!=n3;++i)
        for(int j=0;j!=n3;++j)
            data.H[index3[i]][index3[j]]+=
                data3.G[i].Transpose_Times(V.H[2][2]*data3.G[j]);

    for(int i1=0;i1!=n1;++i1)
        for(int i2=0;i2!=n2;++i2){
            const MATRIX<T,2> x=data1.G[i1].Transpose_Times(V.H[0][1]*data2.G[i2]);
            data.H[index1[i1]][index2[i2]]+=x;
            data.H[index2[i2]][index1[i1]]+=x.Transposed();}
    for(int i1=0;i1!=n1;++i1)
        for(int i3=0;i3!=n3;++i3){
            const MATRIX<T,2> x=data1.G[i1].Transpose_Times(V.H[0][2]*data3.G[i3]);
            data.H[index1[i1]][index3[i3]]+=x;
            data.H[index3[i3]][index1[i1]]+=x.Transposed();}
    for(int i2=0;i2!=n2;++i2)
        for(int i3=0;i3!=n3;++i3){
            const MATRIX<T,2> x=data2.G[i2].Transpose_Times(V.H[1][2]*data3.G[i3]);
            data.H[index2[i2]][index3[i3]]+=x;
            data.H[index3[i3]][index2[i2]]+=x.Transposed();}
}

template<class T,int m,int n>
void Add_Data(VOL_DATA<T,2,4>& data,const VOL_DATA<T,2,m>& V,int const (&datai)[n],int const (&Vi)[n])
{
    data.V+=V.V;
    for(int i=0;i!=n;++i){
        data.G[datai[i]]+=V.G[Vi[i]];
        for(int j=0;j!=n;++j)
            data.H[datai[i]][datai[j]]+=V.H[Vi[i]][Vi[j]];
    }
}

const int vec_a[1]={0}, vec_c[1]={2}, vec_d[1]={3}, vec_abc[3]={0,1,2}, vec_abd[3]={0,1,3}, vec_cda[3]={2,3,0}, vec_cdb[3]={2,3,1}, vec_abcd[4]={0,1,2,3};
template<class T,class TV> void Case_CCAA(VOL_DATA<T,2,4>& data,const TV& X0,const TV& /*A*/,const TV& /*B*/,const TV& C,const TV& D)
{
    // A X   X   X B
    //     D   C
    //       O
//    trap_cases.Append(1);
#if 0
    DATA<T,2,1> DC;
    Data_From_Dof(DC,C-X0);
    DATA<T,2,1> DD;
    Data_From_Dof(DD,D-X0);

    VOL_DATA<T,2,2> V;
    Area_From_Points(V,DC.V,DD.V);
    Combine_Data(data,V,DC,DD,vec_c,vec_d);
#else // #if 0|1
    static const int cd[]={2,3};
    static const int _01[]={0,1};
    {VOL_DATA<T,2,3> V;
    Triangle_Area(V,C,D,X0);
    Add_Data(data,V,cd,_01);}
#endif // #if 0|1
}

template<class T,class TV> void Case_CCAB(VOL_DATA<T,2,4>& data,const TV& X0,const TV& A,const TV& B,const TV& C,const TV& D)
{
    // D
    //   X
    // A P Q   X B
    //       C
    //     O
//    trap_cases.Append(2);
#if 0
    DATA<T,2,4> Q;
    Intersect_Segments(Q,A-X0,B-X0,C-X0,D-X0);
    DATA<T,2,3> P;
    Intersect_Segment_Point(P,A-X0,B-X0,D-X0);
    DATA<T,2,1> DC;
    Data_From_Dof(DC,C-X0);

    VOL_DATA<T,2,3> V;
    Area_From_Points(V,DC.V,Q.V,P.V);
    Combine_Data(data,V,DC,Q,P,vec_c,vec_abcd,vec_abd);
#else // #if 0|1
    static const int cd[]={2,3};
    static const int adb[]={0,3,1};
    static const int dab[]={3,0,1};
    static const int abdc[]={0,1,3,2};
    static const int _01[]={0,1};
    static const int _012[]={0,1,2};
    static const int _123[]={1,2,3};
    static const int _0123[]={0,1,2,3};
    {VOL_DATA<T,2,3> V;
    Triangle_Area(V,A,D,B);
    Add_Data(data,V,adb,_012);}
    {VOL_DATA<T,2,4> V;
    Triangle_Area(V,X0,D,A,B);
    Add_Data(data,V,dab,_123);}
    {VOL_DATA<T,2,4> V;
    Triangle_Area(V,A,B,D,C);
    Add_Data(data,V,abdc,_0123);}
    {VOL_DATA<T,2,3> V;
    Triangle_Area(V,C,D,X0);
    Add_Data(data,V,cd,_01);}
#endif // #if 0|1
}

template<class T,class TV> void Case_CCBB(VOL_DATA<T,2,4>& data,const TV& X0,const TV& A,const TV& B,const TV& C,const TV& D)
{
    //   D       C
    // A   P   P   B
    //       O
//    trap_cases.Append(3);
#if 0
    DATA<T,2,3> P1;
    Intersect_Segment_Point(P1,A-X0,B-X0,C-X0);
    DATA<T,2,3> P2;
    Intersect_Segment_Point(P2,A-X0,B-X0,D-X0);

    VOL_DATA<T,2,2> V;
    Area_From_Points(V,P1.V,P2.V);
    Combine_Data(data,V,P1,P2,vec_abc,vec_abd);
#else // #if 0|1
    static const int ba[]={1,0};
    static const int abc[]={0,1,2};
    static const int dab[]={3,0,1};
    static const int _01[]={0,1};
    static const int _023[]={0,2,3};
    static const int _013[]={0,1,3};
    {VOL_DATA<T,2,4> V;
    Triangle_Area(V,D,X0,A,B);
    Add_Data(data,V,dab,_023);}
    {VOL_DATA<T,2,4> V;
    Triangle_Area(V,A,B,X0,C);
    Add_Data(data,V,abc,_013);}
    {VOL_DATA<T,2,3> V;
    Triangle_Area(V,B,A,X0);
    Add_Data(data,V,ba,_01);}
#endif // #if 0|1
}

template<class T,class TV> void Case_BCAC(VOL_DATA<T,2,4>& data,const TV& X0,const TV& A,const TV& /*B*/,const TV& C,const TV& D)
{
    //   A         B
    // D   P   C
    //       O
//    trap_cases.Append(4);
#if 0
    DATA<T,2,3> P;
    Intersect_Segment_Point(P,C-X0,D-X0,A-X0);
    DATA<T,2,1> DC;
    Data_From_Dof(DC,C-X0);

    VOL_DATA<T,2,2> V;
    Area_From_Points(V,DC.V,P.V);
    Combine_Data(data,V,DC,P,vec_c,vec_cda);
#else // #if 0|1
    static const int acd[]={0,2,3};
    static const int _023[]={0,2,3};
    {VOL_DATA<T,2,4> V;
    Triangle_Area(V,A,X0,C,D);
    Add_Data(data,V,acd,_023);}
#endif // #if 0|1
}

template<class T,class TV> void Case_BCBC(VOL_DATA<T,2,4>& data,const TV& X0,const TV& A,const TV& B,const TV& C,const TV& D)
{
    // A
    //   X
    // D P Q   C
    //       P
    //     O   B
//    trap_cases.Append(5);
#if 0
    DATA<T,2,4> Q;
    Intersect_Segments(Q,A-X0,B-X0,C-X0,D-X0);
    DATA<T,2,3> P1;
    Intersect_Segment_Point(P1,C-X0,D-X0,A-X0);
    DATA<T,2,3> P2;
    Intersect_Segment_Point(P2,A-X0,B-X0,C-X0);

    VOL_DATA<T,2,3> V;
    Area_From_Points(V,P2.V,Q.V,P1.V);
    Combine_Data(data,V,P2,Q,P1,vec_abc,vec_abcd,vec_cda);
#else // #if 0|1
    static const int ca[]={2,0};
    static const int bac[]={1,0,2};
    static const int acd[]={0,2,3};
    static const int dcab[]={3,2,0,1};
    static const int _01[]={0,1};
    static const int _012[]={0,1,2};
    static const int _123[]={1,2,3};
    static const int _0123[]={0,1,2,3};
    {VOL_DATA<T,2,4> V;
    Triangle_Area(V,B,A,C,X0);
    Add_Data(data,V,bac,_012);}
    {VOL_DATA<T,2,4> V;
    Triangle_Area(V,X0,A,C,D);
    Add_Data(data,V,acd,_123);}
    {VOL_DATA<T,2,4> V;
    Triangle_Area(V,D,C,A,B);
    Add_Data(data,V,dcab,_0123);}
    {VOL_DATA<T,2,3> V;
    Triangle_Area(V,C,A,X0);
    Add_Data(data,V,ca,_01);}
#endif // #if 0|1
}

template<class T,class TV> void Case_ACAC(VOL_DATA<T,2,4>& data,const TV& X0,const TV& A,const TV& B,const TV& C,const TV& D)
{
    // D
    //   A  Q  X  X  X  X  B
    //            C
    //      O
//    trap_cases.Append(6);
#if 0
    DATA<T,2,4> Q;
    Intersect_Segments(Q,A-X0,B-X0,C-X0,D-X0);
    DATA<T,2,1> DA;
    Data_From_Dof(DA,A-X0);
    DATA<T,2,1> DC;
    Data_From_Dof(DC,C-X0);

    VOL_DATA<T,2,3> V;
    Area_From_Points(V,DC.V,Q.V,DA.V);
    Combine_Data(data,V,DC,Q,DA,vec_c,vec_abcd,vec_a);
#else // #if 0|1
    static const int ca[]={2,0};
    static const int bacd[]={1,0,2,3};
    static const int _01[]={0,1};
    static const int _0123[]={0,1,2,3};
    {VOL_DATA<T,2,3> V;
    Triangle_Area(V,C,A,X0);
    Add_Data(data,V,ca,_01);}
    {VOL_DATA<T,2,4> V;
    Triangle_Area(V,B,A,C,D);
    Add_Data(data,V,bacd,_0123);}
#endif // #if 0|1
}


template<class T,class TV> void PhysBAM::ORIGIN_AREAS::Volume_From_Simplices(VOL_DATA<T,2,4>& data,TV const & X0,TV const (&LA)[4])
{
    TV A=LA[0],B=LA[1],C=LA[2],D=LA[3];
    if(A==B && C==D) return;
    int index[4]={0,1,2,3};
    T sign=1,tol=(T)1e-10;

    // Enforce consistency with: AD BC + AB CD - AC BD = 0
    T AB=TV::Cross_Product(A-X0,B-X0).x;
    T CD=TV::Cross_Product(C-X0,D-X0).x;
    T AC=TV::Cross_Product(A-X0,C-X0).x;
    T BD=TV::Cross_Product(B-X0,D-X0).x;
    T AD=TV::Cross_Product(A-X0,D-X0).x;
    T BC=TV::Cross_Product(B-X0,C-X0).x;
    int sab=AB<-tol?-1:AB>tol;
    int scd=CD<-tol?-1:CD>tol;
    int sac=AC<-tol?-1:AC>tol;
    int sbd=BD<-tol?-1:BD>tol;
    int sad=AD<-tol?-1:AD>tol;
    int sbc=BC<-tol?-1:BC>tol;
    bool bins[3]={false,false,false};
    bins[sad*sbc+1]=true;
    bins[sab*scd+1]=true;
    bins[-sac*sbd+1]=true;

    // Inconsistent; force each term to zero
    if(bins[0] != bins[2]){
        if(sad && sbc){
            if(fabs(AD)<fabs(BC)) sad=0;
            else sbc=0;}
        if(sab && scd){
            if(fabs(AB)<fabs(CD)) sab=0;
            else scd=0;}
        if(sac && sbd){
            if(fabs(AC)<fabs(BD)) sac=0;
            else sbd=0;}}

    if(sab<0){
        exchange(A,B);
        exchange(index[0],index[1]);
        sign=-sign;
        exchange(sac,sbc);
        exchange(sbd,sad);
        sab=-sab;}
    if(scd<0){
        exchange(C,D);
        exchange(index[2],index[3]);
        sign=-sign;
        exchange(sac,sad);
        exchange(sbd,sbc);
        scd=-scd;}

    if(A==B || sac<0 || (sac==0 && sbd<0)){
        exchange(A,C);
        exchange(B,D);
        exchange(index[0],index[2]);
        exchange(index[1],index[3]);
        exchange(sab,scd);
        exchange(sad,sbc);
        sad=-sad;
        sbc=-sbc;
        sac=-sac;
        sbd=-sbd;}

    VOL_DATA<T,2,4> tdata;
    Clear(tdata);

    if(C==D){
        if(sbc>0 || sbd>0 || sac<0 || sad<0) return;
        if(TV::Cross_Product(A-C,B-C).x>-tol) Case_CCAA(tdata,X0,A,B,C,D); // CCAA
        else Case_CCBB(tdata,X0,A,B,C,D); // CCBB

        data.V=sign*tdata.V;
        for(int i=0;i<4;i++) data.G[index[i]]=sign*tdata.G[i];
        for(int i=0;i<4;i++) for(int k=0;k<4;k++) data.H[index[i]][index[k]]=sign*tdata.H[i][k];}

    // NOTE: (DxA) (CxB) + (CxD) (AxB) + (AxC) (DxB) = 0
    // (after subtracing off X0 from each of A,B,C,D)

    if(sbc>=0) return; // CCCC
    if(TV::Cross_Product(A-C,B-C).x>-tol){
        if(sbd>=0){
            if(TV::Cross_Product(C-B,D-B).x>tol){ // CAAC
                exchange(A,B);exchange(index[0],index[1]);sign=-sign;
                if(TV::Cross_Product(A-X0,B-X0).x<0) sign=-sign;
                Case_ACAC(tdata,X0,A,B,C,D);}
            else{ // CBAC
                exchange(A,B);
                exchange(index[0],index[1]);sign=-sign;
                if(TV::Cross_Product(A-X0,B-X0).x<0) sign=-sign;
                Case_BCAC(tdata,X0,A,B,C,D);}}
        else{
            if(TV::Cross_Product(A-D,B-D).x>-tol) Case_CCAA(tdata,X0,A,B,C,D); // CCAA
            else Case_CCAB(tdata,X0,A,B,C,D);}} // CCAB
    else{
        if(sbd>0){
            if(TV::Cross_Product(C-B,D-B).x>-tol){ // CABC
                exchange(A,B);
                exchange(index[0],index[1]);sign=-sign;
                exchange(A,C);exchange(B,D);
                exchange(index[0],index[2]);exchange(index[1],index[3]);
                if(TV::Cross_Product(A-X0,B-X0).x<0) sign=-sign;
                Case_BCAC(tdata,X0,A,B,C,D);}
            else{ // CBBC
                exchange(A,B);
                exchange(index[0],index[1]);sign=-sign;
                if(TV::Cross_Product(A-X0,B-X0).x<0) sign=-sign;
                Case_BCBC(tdata,X0,A,B,C,D);}}
        else{
            if(TV::Cross_Product(A-D,B-D).x>-tol){ // CCBA
                exchange(C,D);
                exchange(index[2],index[3]);sign=-sign;
                Case_CCAB(tdata,X0,A,B,C,D);}
            else Case_CCBB(tdata,X0,A,B,C,D);}} // CCBB

    data.V=sign*tdata.V;
    for(int i=0;i<4;i++) data.G[index[i]]=sign*tdata.G[i];
    for(int i=0;i<4;i++) for(int k=0;k<4;k++) data.H[index[i]][index[k]]=sign*tdata.H[i][k];
}
namespace PhysBAM{
template void ORIGIN_AREAS::Clear<float,3,6>(VOL_DATA<float,3,6>&);
template void ORIGIN_AREAS::Clear<float,2,4>(VOL_DATA<float,2,4>&);
template void ORIGIN_AREAS::Volume_From_Simplices<float,VECTOR<float,2> >(VOL_DATA<float,2,4>&,VECTOR<float,2> const &,VECTOR<float,2> const (&)[4]);
template void ORIGIN_AREAS::Clear<double,3,6>(VOL_DATA<double,3,6>&);
template void ORIGIN_AREAS::Clear<double,2,4>(VOL_DATA<double,2,4>&);
template void ORIGIN_AREAS::Volume_From_Simplices<double,VECTOR<double,2> >(VOL_DATA<double,2,4>&,VECTOR<double,2> const &,VECTOR<double,2> const (&)[4]);
}

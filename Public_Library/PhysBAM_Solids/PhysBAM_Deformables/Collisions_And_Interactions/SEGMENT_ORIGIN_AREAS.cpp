#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/SEGMENT_ORIGIN_AREAS.h>
using namespace PhysBAM;
using namespace SEGMENT_ORIGIN_AREAS;

template<class T,int n> void Clear(DATA<T,n>& data)
{
    data.V=T();
    for(int i=0;i<2*n;i++) data.G[i]=T();
    for(int i=0;i<2*n;i++) for(j=0;j<2*n;j++) data.H[i][j]=T();
}

template<class TV> void Data_From_Dof(DATA<TV,1>& data,const TV& a)
{
    data.V=a;
    data.G[0]=TV(1,0);
    data.G[1]=TV(0,1);
    for(int i=0;i<2;i++) for(j=0;j<2;j++) data.H[i][j]=TV();
}

template<class T,class TV> void Case_CCAA(DATA<T,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    Clear(data);

    DATA<TV,3> DC;
    Data_From_Dof(DC,C);

    DATA<TV,3> DD;
    Data_From_Dof(DD,D);

    DATA<T,2> V;
    Area_From_Points(V,DC.V,DD.V);
    Combine_Data(data,V,DC,DD,VECTOR<int,1>(3),VECTOR<int,1>(4));
}

template<class T,class TV> void Case_CCAB(DATA<T,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    Clear(data);

    DATA<TV,4> Q;
    Intersect_Segments(Q,A,B,C,D);

    DATA<TV,3> P;
    Intersect_Segment_Point(P,A,B,D);

    DATA<TV,3> DC;
    Data_From_Dof(DC,C);

    DATA<T,2> V;
    Area_From_Points(V,P.V,Q.V);
    Combine_Data(data,V,P,Q,VECTOR<int,3>(1,2,4),VECTOR<int,4>(1,2,3,4));

    Area_From_Points(V,Q.V,DC.V);
    Combine_Data(data,V,Q,DC,VECTOR<int,4>(1,2,3,4),VECTOR<int,1>(3));
}

template<class T,class TV> void Case_CCBB(DATA<T,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    Clear(data);

    DATA<TV,3> P1;
    Intersect_Segment_Point(P1,A,B,C);

    DATA<TV,3> P2;
    Intersect_Segment_Point(P2,A,B,D);

    DATA<T,2> V;
    Area_From_Points(V,P1.V,P2.V);
    Combine_Data(data,V,P1,P2,VECTOR<int,3>(1,2,3),VECTOR<int,3>(1,2,4));
}

template<class T,class TV> void Case_BCAC(DATA<T,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    Clear(data);

    DATA<TV,3> P;
    Intersect_Segment_Point(P,A,B,C);

    DATA<TV,3> DA;
    Data_From_Dof(DA,A);

    DATA<T,2> V;
    Area_From_Points(V,DA.V,P.V);
    Combine_Data(data,V,DA,P,VECTOR<int,1>(1),VECTOR<int,3>(1,2,3));
}

template<class T,class TV> void Case_BCBC(DATA<T,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    Clear(data);

    DATA<TV,4> Q;
    Intersect_Segments(Q,A,B,C,D);

    DATA<TV,3> P;
    Intersect_Segment_Point(P,C,D,B);

    DATA<T,2> V;
    Area_From_Points(V,Q.V,P.V);

    Combine_Data(data,V,Q,P,VECTOR<int,4>(1,2,3,4),VECTOR<int,3>(3,4,2));
}

template<class T,class TV> void Case_ACAC(DATA<T,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    Clear(data);

    DATA<TV,4> Q;
    Intersect_Segments(Q,A,B,C,D);

    DATA<TV,3> DA;
    Data_From_Dof(DA,A);

    DATA<TV,3> DC;
    Data_From_Dof(DC,C);

    DATA<T,2> V;
    Area_From_Points(V,DA.V,Q.V);
    Combine_Data(data,V,DA,Q,VECTOR<int,1>(1),VECTOR<int,4>(1,2,3,4));

    Area_From_Points(V,Q.V,DC.V);
    Combine_Data(data,V,Q,DC,VECTOR<int,4>(1,2,3,4),VECTOR<int,1>(3));
}





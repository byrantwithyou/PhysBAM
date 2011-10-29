#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/SEGMENT_ORIGIN_AREAS.h>
using namespace PhysBAM;
using namespace SEGMENT_ORIGIN_AREAS;

template<class T,int n> void Clear(DATA<T,n>& data)
{
    data.V=T();
    for(int i=0;i<2*n;i++) data.G[i]=T();
    for(int i=0;i<2*n;i++) for(int j=0;j<2*n;j++) data.H[i][j]=T();
}

template<class TV> void Data_From_Dof(DATA<TV,1>& data,const TV& a)
{
    data.V=a;
    data.G[0]=TV(1,0);
    data.G[1]=TV(0,1);
    for(int i=0;i<2;i++) for(int j=0;j<2;j++) data.H[i][j]=TV();
}

template<class TV> POINT_CASE Classify_Point(const TV& A,const TV& B,const TV& P)
{
    if(TV::Cross_Product(A,P).x<0) return outside;
    if(TV::Cross_Product(P,B).x<0) return outside;
    if(TV::Cross_Product(A-P,B-P).x<0) return beyond;
    return inside;
}

template<class TV> void Intersect_Segment_Point(DATA<TV,3>& data,const TV& A,const TV& B,const TV& P)
{
    bool Do_Intersect = false;
    typename TV::value_type EPS = 1e-10;
    typename TV::value_type mua,mub;
    typename TV::value_type denom,numera,numerb;
    typename TV::value_type denom_square;

    denom=(0.0-P.y)*(B.x-A.x)-(0.0-P.x)*(B.y-A.y);
    numera=(0.0-P.x)*(A.y-P.y)-(0.0-P.y)*(A.x-P.x);
    numerb=(B.x-A.x)*(A.y-P.y)-(B.y-A.y)*(A.x-P.x);

    // Are the segments coincident?
    if(std::abs(numera)<EPS && std::abs(numerb)<EPS && std::abs(denom)<EPS){
        data.V.x=(A.x+B.x)/2.0;
        data.V.y=(A.y+B.y)/2.0;
    }
    // Are the segments parallel?
    else if(std::abs(denom)<EPS) {
        data.V=TV();
    }
    else{
        mua=numera/denom;
        mub=numerb/denom;
        // Do they really intersect?
        if(mua<0.0||mua>1.0||mub<0.0||mub>1.0){
            data.V=TV();
        }
        else{
            Do_Intersect=true;
            data.V.x=A.x+mua*(B.x-A.x);
            data.V.y=A.y+mua*(B.y-A.y);
        }
    }
    
    if(Do_Intersect==false){
        //
        // What to do with G and H??
        //
    }
    else{
        denom_square=denom*denom;
        data.G[0].x = P.x*(A.y-B.y)*(B.x*P.y-P.x*B.y)/denom_square;
        data.G[1].x =  -P.x*(A.x-B.x)*(B.x*P.y-P.x*B.y)/denom_square;
        data.G[2].x = -P.x*(A.y-B.y)*(A.x*P.y-P.x*A.y)/denom_square;
        data.G[3].x = P.x*(A.x-B.x)*(A.x*P.y-P.x*A.y)/denom_square;
        data.G[4].x =  P.y*(A.x-B.x)*(A.x*B.y-B.x*A.y)/denom_square;
        data.G[5].x = -P.x*(A.x-B.x)*(A.x*B.y-B.x*A.y)/denom_square;

        data.G[0].y = P.y*(A.y-B.y)*(B.x*P.y-P.x*B.y)/denom_square;
        data.G[1].y = -P.y*(A.x-B.x)*(B.x*P.y-P.x*B.y)/denom_square;
        data.G[2].y = -P.y*(A.y-B.y)*(A.x*P.y-P.x*A.y)/denom_square;
        data.G[3].y = P.y*(A.x-B.x)*(A.x*P.y-P.x*A.y)/denom_square;
        data.G[4].y =  P.y*(A.y-B.y)*(A.x*B.y-B.x*A.y)/denom_square;
        data.G[5].y = -P.x*(A.y-B.y)*(A.x*B.y-B.x*A.y)/denom_square;
    }

    
    
}

template<class TV> void Intersect_Segments(DATA<TV,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    bool Do_Intersect = false;
    typename TV::value_type EPS = 1e-10;
    typename TV::value_type mua,mub;
    typename TV::value_type denom,numera,numerb;
    typename TV::value_type denom_square;

    denom=(D.y-C.y)*(B.x-A.x)-(D.x-C.x)*(B.y-A.y);
    numera=(D.x-C.x)*(A.y-C.y)-(D.y-C.y)*(A.x-C.x);
    numerb=(B.x-A.x)*(A.y-C.y)-(B.y-A.y)*(A.x-C.x);

    // Are the segments coincident?
    if(std::abs(numera)<EPS && std::abs(numerb)<EPS && std::abs(denom)<EPS){
        data.V.x=(A.x+B.x)/2.0;
        data.V.y=(A.y+B.y)/2.0;
    }
    // Are the segments parallel?
    else if(std::abs(denom)<EPS) {
        data.V=TV();
    }
    else{
        mua=numera/denom;
        mub=numerb/denom;
        // Do they really intersect?
        if(mua<0.0||mua>1.0||mub<0.0||mub>1.0){
            data.V=TV();
        }
        else{
            Do_Intersect=true;
            data.V.x=A.x+mua*(B.x-A.x);
            data.V.y=A.y+mua*(B.y-A.y);
        }
    }
    
    if(Do_Intersect==false){
        //
        // What to do with G and H??
        //
    }
    else{
        denom_square=denom*denom;
        data.G[0].x = (A.y-B.y)*(C.x-D.x)*(B.y*D.x-C.y*D.x+B.x*C.y-C.x*B.y-B.x*D.y+C.x*D.y)/denom_square;
        data.G[1].x = -(A.x-B.x)*(C.x-D.x)*(B.y*D.x-C.y*D.x+B.x*C.y-C.x*B.y-B.x*D.y+C.x*D.y)/denom_square;
        data.G[2].x = -(A.y-B.y)*(C.x-D.x)*(A.y*D.x-C.y*D.x+A.x*C.y-C.x*A.y-A.x*D.y+C.x*D.y)/denom_square;
        data.G[3].x = (A.x-B.x)*(C.x-D.x)*(A.y*D.x-C.y*D.x+A.x*C.y-C.x*A.y-A.x*D.y+C.x*D.y)/denom_square;
        data.G[4].x = (A.x-B.x)*(C.y-D.y)*(A.y*D.x-B.y*D.x+A.x*B.y-B.x*A.y-A.x*D.y+B.x*D.y)/denom_square;
        data.G[5].x = -(A.x-B.x)*(C.x-D.x)*(A.y*D.x-B.y*D.x+A.x*B.y-B.x*A.y-A.x*D.y+B.x*D.y)/denom_square;
        data.G[6].x = -(A.x-B.x)*(C.y-D.y)*(A.x*B.y-B.x*A.y-A.x*C.y+C.x*A.y+B.x*C.y-C.x*B.y)/denom_square;
        data.G[7].x = (A.x-B.x)*(C.x-D.x)*(A.x*B.y-B.x*A.y-A.x*C.y+C.x*A.y+B.x*C.y-C.x*B.y)/denom_square;
        
        data.G[0].y = (A.y-B.y)*(C.y-D.y)*(B.y*D.x-C.y*D.x+B.x*C.y-C.x*B.y-B.x*D.y+C.x*D.y)/denom_square;
        data.G[1].y = -(A.x-B.x)*(C.y-D.y)*(B.y*D.x-C.y*D.x+B.x*C.y-C.x*B.y-B.x*D.y+C.x*D.y)/denom_square;
        data.G[2].y = -(A.y-B.y)*(C.y-D.y)*(A.y*D.x-C.y*D.x+A.x*C.y-C.x*A.y-A.x*D.y+C.x*D.y)/denom_square;
        data.G[3].y = (A.x-B.x)*(C.y-D.y)*(A.y*D.x-C.y*D.x+A.x*C.y-C.x*A.y-A.x*D.y+C.x*D.y)/denom_square;
        data.G[4].y = (A.y-B.y)*(C.y-D.y)*(A.y*D.x-B.y*D.x+A.x*B.y-B.x*A.y-A.x*D.y+B.x*D.y)/denom_square;
        data.G[5].y = -(A.y-B.y)*(C.x-D.x)*(A.y*D.x-B.y*D.x+A.x*B.y-B.x*A.y-A.x*D.y+B.x*D.y)/denom_square;
        data.G[6].y = -(A.y-B.y)*(C.y-D.y)*(A.x*B.y-B.x*A.y-A.x*C.y+C.x*A.y+B.x*C.y-C.x*B.y)/denom_square;
        data.G[7].y = (A.y-B.y)*(C.x-D.x)*(A.x*B.y-B.x*A.y-A.x*C.y+C.x*A.y+B.x*C.y-C.x*B.y)/denom_square;
    }
}

template<class T,class TV> void Area_From_Points(DATA<T,2>& data,const TV& A,const TV& B)
{
    data.v=Cross_Product(A,B);
    
    data.G[0]=B.y; data.G[1]=-B.x; data.G[2]=-A.y; data.G[3]=A.x;
    
    data.H[0][0]=0.0;   data.H[0][1]=0.0;   data.H[0][2]=0.0;   data.H[0][3]=1.0;
    data.H[1][0]=0.0;   data.H[1][1]=0.0;   data.H[1][2]=-1.0;  data.H[1][3]=0.0;
    data.H[2][0]=0.0;   data.H[2][1]=-1.0;  data.H[2][2]=0.0;   data.H[2][3]=0.0;
    data.H[3][0]=1.0;   data.H[3][1]=0.0;   data.H[3][2]=0.0;   data.H[3][3]=0.0;
}

template<class T,class TV> void Area_From_Segments(DATA<T,4>& data,TV A,TV B,TV C,TV D)
{
    int index[4]={1,2,3,4};
    T sign=1;
    T OAB=TV::Cross_Product(A,B).x;
    T OCD=TV::Cross_Product(C,D).x;

    if(OAB<0){exchange(A,B);exchange(index[0],index[1]);sign=-sign;}
    if(OCD<0){exchange(C,D);exchange(index[2],index[3]);sign=-sign;}

    POINT_CASE case_a=Classify_Point(C,D,A);
    POINT_CASE case_b=Classify_Point(C,D,B);
    POINT_CASE case_c=Classify_Point(A,B,C);
    POINT_CASE case_d=Classify_Point(A,B,D);
    if(case_b<case_a){exchange(A,B);exchange(index[0],index[1]);sign=-sign;}
    if(case_d<case_c){exchange(C,D);exchange(index[2],index[3]);sign=-sign;}
    if(case_a<case_c || (case_a==case_c && case_b<case_d)){exchange(A,C);exchange(B,D);exchange(index[0],index[2]);exchange(index[1],index[3]);}

    DATA<T,4> tdata;
    if(case_a==outside){
        if(case_b==outside && case_c==outside && case_d==outside) return;
        PHYSBAM_ASSERT(case_b==outside && case_c!=outside && case_d!=outside);
        if(case_c==beyond) Case_CCBB(tdata,A,B,C,D);
        else if(case_d==beyond) Case_CCAB(tdata,A,B,C,D);
        else Case_CCAA(tdata,A,B,C,D);}
    else if(case_a==beyond){
        PHYSBAM_ASSERT(case_b==outside && case_d==outside);
        if(case_c==inside) Case_BCAC(tdata,A,B,C,D);
        else Case_BCBC(tdata,A,B,C,D);}
    else{
        PHYSBAM_ASSERT(case_b==outside && case_c==inside && case_d==outside);
        Case_ACAC(tdata,A,B,C,D);}

    data.V=sign*tdata.V;
    for(int i=0;i<4;i++) for(int j=0;j<2;j++) data.G[2*index[i]+j]=sign*data.G[2*i+j];
    for(int i=0;i<4;i++) for(int k=0;k<4;k++) for(int j=0;j<2;j++) for(int m=0;m<2;m++) data.G[2*index[i]+j][2*index[k]+m]=sign*data.G[2*i+j][2*k+m];
}

template<class T,class TV,int m,int n> void
Combine_Data(DATA<T,4>& data,const DATA<T,2>& V,const DATA<TV,m>& data_m,const DATA<TV,n>& data_n,const int index_m[m],const int index_n[n])
{
    data.V=V.V;
    for(int j=0;j<m;j++)
        for(int k=0;k<2;k++)
            for(int i=0;i<2;i++)
                data.G[2*index_m[j]+k]+=V.G[i]*data_m.G[2*j+k](i+1);

    for(int j=0;j<n;j++)
        for(int k=0;k<2;k++)
            for(int i=0;i<2;i++)
                data.G[2*index_n[j]+k]+=V.G[i+2]*data_n.G[2*j+k](i+1);

    for(int j=0;j<m;j++)
        for(int k=0;k<2;k++)
            for(int s=0;j<m;j++)
                for(int t=0;k<2;k++)
                    for(int i=0;i<2;i++)
                        for(int r=0;i<2;i++)
                            data.H[2*index_m[j]+k][2*index_m[s]+t]+=V.H[i][r]*data_m.G[2*j+k](i+1)*data_m.G[2*s+t](r+1);

    for(int j=0;j<m;j++)
        for(int k=0;k<2;k++)
            for(int s=0;j<m;j++)
                for(int t=0;k<2;k++)
                    for(int i=0;i<2;i++)
                        data.H[2*index_m[j]+k][2*index_m[s]+t]+=V.G[i]*data_m.H[2*j+k][2*s+t](i+1);

    for(int j=0;j<n;j++)
        for(int k=0;k<2;k++)
            for(int s=0;j<n;j++)
                for(int t=0;k<2;k++)
                    for(int i=0;i<2;i++)
                        for(int r=0;i<2;i++)
                            data.H[2*index_n[j]+k][2*index_n[s]+t]+=V.H[i+2][r+2]*data_n.G[2*j+k](i+1)*data_n.G[2*s+t](r+1);

    for(int j=0;j<n;j++)
        for(int k=0;k<2;k++)
            for(int s=0;j<n;j++)
                for(int t=0;k<2;k++)
                    for(int i=0;i<2;i++)
                        data.H[2*index_n[j]+k][2*index_n[s]+t]+=V.G[i+2]*data_n.H[2*j+k][2*s+t](i+1);
}

const int vec_a[1]={0}, vec_c[1]={2}, vec_d[1]={3}, vec_abc[3]={0,1,2}, vec_abd[3]={0,1,3}, vec_cdb[3]={2,3,1}, vec_abcd[4]={0,1,2,3};
template<class T,class TV> void Case_CCAA(DATA<T,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    DATA<TV,3> DC;
    Data_From_Dof(DC,C);

    DATA<TV,3> DD;
    Data_From_Dof(DD,D);

    DATA<T,2> V;
    Area_From_Points(V,DC.V,DD.V);
    Combine_Data(data,V,DC,DD,vec_c,vec_d);
}

template<class T,class TV> void Case_CCAB(DATA<T,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    DATA<TV,4> Q;
    Intersect_Segments(Q,A,B,C,D);

    DATA<TV,3> P;
    Intersect_Segment_Point(P,A,B,D);

    DATA<TV,3> DC;
    Data_From_Dof(DC,C);

    DATA<T,2> V;
    Area_From_Points(V,P.V,Q.V);
    Combine_Data(data,V,P,Q,vec_abd,vec_abcd);

    Area_From_Points(V,Q.V,DC.V);
    Combine_Data(data,V,Q,DC,vec_abcd,vec_c);
}

template<class T,class TV> void Case_CCBB(DATA<T,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    DATA<TV,3> P1;
    Intersect_Segment_Point(P1,A,B,C);

    DATA<TV,3> P2;
    Intersect_Segment_Point(P2,A,B,D);

    DATA<T,2> V;
    Area_From_Points(V,P1.V,P2.V);
    Combine_Data(data,V,P1,P2,vec_abc,vec_abd);
}

template<class T,class TV> void Case_BCAC(DATA<T,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    DATA<TV,3> P;
    Intersect_Segment_Point(P,A,B,C);

    DATA<TV,3> DA;
    Data_From_Dof(DA,A);

    DATA<T,2> V;
    Area_From_Points(V,DA.V,P.V);
    Combine_Data(data,V,DA,P,vec_a,vec_abc);
}

template<class T,class TV> void Case_BCBC(DATA<T,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    DATA<TV,4> Q;
    Intersect_Segments(Q,A,B,C,D);

    DATA<TV,3> P;
    Intersect_Segment_Point(P,C,D,B);

    DATA<T,2> V;
    Area_From_Points(V,Q.V,P.V);

    Combine_Data(data,V,Q,P,vec_abcd,vec_cdb);
}

template<class T,class TV> void Case_ACAC(DATA<T,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    DATA<TV,4> Q;
    Intersect_Segments(Q,A,B,C,D);

    DATA<TV,3> DA;
    Data_From_Dof(DA,A);

    DATA<TV,3> DC;
    Data_From_Dof(DC,C);

    DATA<T,2> V;
    Area_From_Points(V,DA.V,Q.V);
    Combine_Data(data,V,DA,Q,vec_a,vec_abcd);

    Area_From_Points(V,Q.V,DC.V);
    Combine_Data(data,V,Q,DC,vec_abcd,vec_c);
}





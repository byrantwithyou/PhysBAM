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

template<class TV> void Intersect_Segment_Point(DATA<TV,3>& data,const TV& A,const TV& B,const TV& P)
{

    
    
}

template<class TV> void Intersect_Segments(DATA<TV,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    // Reference for segment-segment intersection: http://paulbourke.net/geometry/lineline2d/
    bool Do_Intersect = false;
    TV::value_type EPS = 1e-10;
    TV::value_type mua,mub;
    TV::value_type denom,numera,numerb;
    TV::value_type denom_square;

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
        
        // H[0][0].x = -2(a-b)(c-d)(C-D)((bD-cD+Bc-Cb-Bd+Cd)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3))  
        // H[1][0].x = H[0][1].x = -(C-D)((bD-cD+Bc-Cb-Bd+Cd)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3)) (aD-bD-Ac-Ca+Ad+Bc+Cb-Bd)
        // H[1][1].x = -2(A-B)(C-D)2((bD-cD+Bc-Cb-Bd+Cd)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3))
        // H[2][0].x = H[0][2].x =  (a-b)(c-d)((C-D)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3)) (aD+bD-2cD+Ac-Ca-Ad+Bc-Cb-Bd+2Cd) 
        // H[2][1].x = H[1][2].x = -((C-D)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3))(A2c2-2A2cd+A2d2+ACcd-aACc-ACd2+aACd-Ac2D-BAc2+AcdD+2BAcd+aAcD-BAd2-aAdD+C2b2-C2bd-aC2b+aC2d-2Cb2D+CbcD-BCbc+CbdD+BCbd+2aCbD-BCcd-aCcD+2BaCc+BCd2-aCdD-2BaCd+b2D2-bcD2+BbcD-BbdD-abD2+Bc2D-BcdD+acD2-2BacD+2BadD)
        // H[2][2].x = -2(a-b)(c-d)(C-D)((aD-cD+Ac-Ca-Ad+Cd)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3))
        // H[3][0].x = H[0][3].x = ((C-D)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3))(B2c2-2B2cd+B2d2+BCcd-bBCc-BCd2+bBCd-Bc2D-ABc2+BcdD+2ABcd+bBcD-ABd2-bBdD+C2a2-C2ad-bC2a+bC2d-2Ca2D+CacD-ACac+CadD+ACad+2bCaD-ACcd-bCcD+2AbCc+ACd2-bCdD-2AbCd+a2D2-acD2+AacD-AadD-baD2+Ac2D-AcdD+bcD2-2AbcD+2AbdD)
        // H[3][1].x = H[1][3].x =  (A-B)(((C-D)2)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3)) (aD+bD-2cD+Ac-Ca-Ad+Bc-Cb-Bd+2Cd)
        // H[3][2].x = H[2][3].x = -(C-D)((aD-cD+Ac-Ca-Ad+Cd)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3)) (aD-bD-Ac-Ca+Ad+Bc+Cb-Bd)
        // H[3][3].x =  -2(A-B)(C-D)2((aD-cD+Ac-Ca-Ad+Cd)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3))
        // H[4][0].x = H[0][4].x = -(a-b)((c-d)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3)) (B2c-B2d-aD2+bD2-2AbD+BaD+AcD+BbD+CaD+AdD-BcD-CbD-BdD-ABc+2ACb-BCa+ABd-BCb-2ACd+2BCd)
        // H[4][1].x = H[1][4].x =  (A-B)((c-d)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3)) (B2c-B2d-aD2+bD2-2AbD+BaD+AcD+BbD+CaD+AdD-BcD-CbD-BdD-ABc+2ACb-BCa+ABd-BCb-2ACd+2BCd)
        // H[4][2].x = H[2][4].x = -(a-b)((c-d)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3)) (A2c-A2d+aD2-bD2+AaD+AbD-2BaD-AcD-CaD-AdD+BcD+CbD+BdD-ACa-ABc-ACb+2BCa+ABd+2ACd-2BCd)
        // H[4][3].x = H[3][4].x =  (A-B)((c-d)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3)) (A2c-A2d+aD2-bD2+AaD+AbD-2BaD-AcD-CaD-AdD+BcD+CbD+BdD-ACa-ABc-ACb+2BCa+ABd+2ACd-2BCd)
        // H[4][4].x =  2(A-B)(a-b)(c-d)((aD-bD+Ab-Ba-Ad+Bd)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3))
        // H[5][0].x = H[0][5].x =  (a-b)((C-D)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3)) (B2c-B2d-aD2+bD2-2AbD+BaD+AcD+BbD+CaD+AdD-BcD-CbD-BdD-ABc+2ACb-BCa+ABd-BCb-2ACd+2BCd)
        // H[5][1].x = H[1][5].x = -(A-B)((C-D)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3)) (B2c-B2d-aD2+bD2-2AbD+BaD+AcD+BbD+CaD+AdD-BcD-CbD-BdD-ABc+2ACb-BCa+ABd-BCb-2ACd+2BCd)
        // H[5][2].x = H[2][5].x =  (a-b)((C-D)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3)) (A2c-A2d+aD2-bD2+AaD+AbD-2BaD-AcD-CaD-AdD+BcD+CbD+BdD-ACa-ABc-ACb+2BCa+ABd+2ACd-2BCd)
        // H[5][3].x = H[3][5].x =  -(A-B)((C-D)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3)) (A2c-A2d+aD2-bD2+AaD+AbD-2BaD-AcD-CaD-AdD+BcD+CbD+BdD-ACa-ABc-ACb+2BCa+ABd+2ACd-2BCd)
        // H[5][4].x = H[4][5].x =  (A-B)((aD-bD+Ab-Ba-Ad+Bd)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3)) (aD-bD-Ac-Ca+Ad+Bc+Cb-Bd)
        // H[5][5].x =  2(A-B)2(C-D)((aD-bD+Ab-Ba-Ad+Bd)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3)) 
        // H[6][0].x = H[0][6].x =  (a-b)((c-d)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3)) (C2a+B2c-C2b-B2d-2AbD+BaD+2AcD+BbD-CaD-2BcD+CbD-ABc+2ACb-BCa+ABd-ACc-BCb-ACd+BCc+BCd)
        // H[6][1].x = H[1][6].x = -(A-B)((c-d)/((aD-bD+Ac-Ca-Ad-Bc+Cb+Bd)3)) (C2a+B2c-C2b-B2d-2AbD+BaD+2AcD+BbD-CaD-2BcD+CbD-ABc+2ACb-BCa+ABd-ACc-BCb-ACd+BCc+BCd)
        // H[6][2].x = H[2][6].x = 
        // H[6][3].x = H[3][6].x = 
        // H[6][4].x = H[4][6].x = 
        // H[6][5].x = H[5][6].x = 
        // H[6][6].x = 
        // H[7][0].x = H[0][7].x = 
        // H[7][1].x = H[1][7].x = 
        // H[7][2].x = H[2][7].x = 
        // H[7][3].x = H[3][7].x = 
        // H[7][4].x = H[4][7].x = 
        // H[7][5].x = H[5][7].x = 
        // H[7][6].x = H[6][7].x = 
        // H[7][7].x = 
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





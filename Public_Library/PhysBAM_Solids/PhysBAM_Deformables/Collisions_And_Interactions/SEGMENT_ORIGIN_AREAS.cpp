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
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) for(int k=0;k<m;k++) data.H[m][i][j]=MATRIX<T,2>();
}

template<class T,class TV> void Data_From_Dof(DATA<T,2,1>& data,const TV& A)
{
    data.V=A;
    data.G[0]=MATRIX<T,2>::Identity_Matrix();
    for(int i=0;i<2;i++) for(int j=0;j<2;j++) data.H[0][i][j]=MATRIX<T,2>();
}

template<class TV> POINT_CASE Classify_Point(const TV& A,const TV& B,const TV& P)
{
    if(TV::Cross_Product(A,P).x<0) return outside;
    if(TV::Cross_Product(P,B).x<0) return outside;
    if(TV::Cross_Product(A-P,B-P).x<0) return beyond;
    puts("INSIDE");
    return inside;
}

template<class T,class TV> void Intersect_Segment_Point(DATA<T,2,3>& data,const TV& A,const TV& B,const TV& P)
{
    T cross1=TV::Cross_Product(A,B).x,cross2=TV::Cross_Product(A-B,P-B).x,den=1/(cross1+cross2),PxB=TV::Cross_Product(P,B).x,PxA=TV::Cross_Product(P,A).x;
    data.V=cross1*den*P;

    TV orthAB=(A-B).Orthogonal_Vector();
    MATRIX<T,2> M=MATRIX<T,2>::Outer_Product(P,sqr(den)*orthAB);
    MATRIX<T,2> GA=PxB*M;
    MATRIX<T,2> GB=-PxA*M;
    MATRIX<T,2> GP=-cross1*M;
    data.G[0]=GA;
    data.G[1]=GB;
    data.G[2]=GP;
}

template<class T,class TV> void Intersect_Segments(DATA<T,2,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    // V = A + z (B - A); compute z
    T ABxCD=TV::Cross_Product(A-B,C-D).x,ADxCD=TV::Cross_Product(A-D,C-D).x,ADxBD=TV::Cross_Product(A-D,B-D).x,ACxBC=TV::Cross_Product(A-C,B-C).x,z=ABxCD/ADxCD;
    data.V=A+z*(B-A);

    // Gradients of z:
    TV oAB=(A-B).Orthogonal_Vector(),oCD=(C-D).Orthogonal_Vector(),g=oCD/sqr(ADxCD);
    TV dzdA=(ABxCD-ADxCD)*g,dzdB=-ABxCD*g,dzdC=ADxBD*g,dzdD=-ACxBC*g;
    MATRIX<T,2> dVdA=MATRIX<T,2>::Outer_Product(B-A,dzdA)+(1-z),dVdB=MATRIX<T,2>::Outer_Product(B-A,dzdB)+z;
    MATRIX<T,2> dVdC=MATRIX<T,2>::Outer_Product(B-A,dzdC),dVdD=MATRIX<T,2>::Outer_Product(B-A,dzdD);
    data.G[0]=dVdA;
    data.G[1]=dVdB;
    data.G[2]=dVdC;
    data.G[3]=dVdD;

    // Hessian components of z:
    MATRIX<T,2> d2zdAdA=(2*(ABxCD-ADxCD)/(ADxCD*sqr(ADxCD)))*MATRIX<T,2>::Outer_Product(oCD,oCD);
    MATRIX<T,2> d2zdBdB=(2*ADxCD/(ADxCD*sqr(ADxCD)))*MATRIX<T,2>::Outer_Product(oCD,oCD);
    MATRIX<T,2> d2zdCdC=(-ADxBD/(ADxCD*sqr(ADxCD)))*(MATRIX<T,2>::Outer_Product(oCD,oAB)+MATRIX<T,2>::Outer_Product(oAB,oCD));
    MATRIX<T,2> d2zdDdD=(-ACxBC/(ADxCD*sqr(ADxCD)))*(MATRIX<T,2>::Outer_Product(oCD,oAB)+MATRIX<T,2>::Outer_Product(oAB,oCD));


    // Hessian 

    bool Do_Intersect = false;
    typename TV::SCALAR EPS = 1e-10;
    typename TV::SCALAR mua,mub;
    typename TV::SCALAR denom,numera,numerb;
    typename TV::SCALAR denom_square;

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
/*        data.G[0].x = (A.y-B.y)*(C.x-D.x)*(B.y*D.x-C.y*D.x+B.x*C.y-C.x*B.y-B.x*D.y+C.x*D.y)/denom_square;
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
*/    }
}

template<class T,class TV> void Area_From_Points(DATA<T,1,2>& data,const TV& A,const TV& B)
{
    data.V=(T).5*TV::Cross_Product(A,B);
    data.G[0]=MATRIX<T,1,2>::Cross_Product_Matrix(B);
    data.G[1]=MATRIX<T,1,2>::Cross_Product_Matrix(-A);
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
        if(case_c==inside) Case_BCAC(tdata,A,B,C,D);
        else Case_BCBC(tdata,A,B,C,D);}
    else{
        PHYSBAM_ASSERT(case_b==outside && case_c==inside && case_d==outside);
        Case_ACAC(tdata,A,B,C,D);}

    data.V=sign*tdata.V;
    for(int i=0;i<4;i++) data.G[index[i]]=sign*tdata.G[i];
    for(int i=0;i<4;i++) for(int k=0;k<4;k++) data.H[0][index[i]][index[k]]=sign*tdata.H[0][i][k];
}

template<class T,int m,int n> void Combine_Data(DATA<T,1,4>& data,const DATA<T,1,2>& V,const DATA<T,2,m>& data_m,const DATA<T,2,n>& data_n,const int index_m[m],
    const int index_n[n])
{
    data.V=V.V;
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

const int vec_a[1]={0}, vec_c[1]={2}, vec_d[1]={3}, vec_abc[3]={0,1,2}, vec_abd[3]={0,1,3}, vec_cdb[3]={2,3,1}, vec_abcd[4]={0,1,2,3};
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
    Area_From_Points(V,P.V,Q.V);
    Combine_Data(data,V,P,Q,vec_abd,vec_abcd);

    Area_From_Points(V,Q.V,DC.V);
    Combine_Data(data,V,Q,DC,vec_abcd,vec_c);
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
    Intersect_Segment_Point(P,A,B,C);

    DATA<T,2,1> DA;
    Data_From_Dof(DA,A);

    DATA<T,1,2> V;
    Area_From_Points(V,DA.V,P.V);
    Combine_Data(data,V,DA,P,vec_a,vec_abc);
}

template<class T,class TV> void Case_BCBC(DATA<T,1,4>& data,const TV& A,const TV& B,const TV& C,const TV& D)
{
    trap_cases.Append(5);
    DATA<T,2,4> Q;
    Intersect_Segments(Q,A,B,C,D);

    DATA<T,2,3> P;
    Intersect_Segment_Point(P,C,D,B);

    DATA<T,1,2> V;
    Area_From_Points(V,Q.V,P.V);

    Combine_Data(data,V,Q,P,vec_abcd,vec_cdb);
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
    Area_From_Points(V,DA.V,Q.V);
    Combine_Data(data,V,DA,Q,vec_a,vec_abcd);

    Area_From_Points(V,Q.V,DC.V);
    Combine_Data(data,V,Q,DC,vec_abcd,vec_c);
}

template void Area_From_Segments<float,VECTOR<float,2> >(DATA<float,1,4>&,VECTOR<float,2>,VECTOR<float,2>,VECTOR<float,2>,
    VECTOR<float,2>);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void Area_From_Segments<double,VECTOR<double,2> >(DATA<double,1,4>&,VECTOR<double,2>,VECTOR<double,2>,VECTOR<double,2>,
    VECTOR<double,2>);
#endif
}
}

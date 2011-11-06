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
    data.n=1;
    data.V=A;
    data.G[0]=MATRIX<T,3>::Identity_Matrix();
    for(int i=0;i<3;i++) data.H[i][0][0]=MATRIX<T,3>();
}

template<class T,class TV> void Volume_From_Points(VOL_DATA<T,3>& data,const TV& A,const TV& B,const TV& C)
{
    data.V=(T)(1./6)*TV::Triple_Product(A,B,C);
    data.G[0]=(T)(1./6)*TV::Cross_Product(B,C);
    data.G[1]=(T)(1./6)*TV::Cross_Product(C,A);
    data.G[2]=(T)(1./6)*TV::Cross_Product(A,B);
    for(int i=0;i<3;i++) data.H[i][i]=MATRIX<T,3>();
    data.H[0][1]=-(T)(1./6)*MATRIX<T,3>::Cross_Product_Matrix(C);
    data.H[0][2]=(T)(1./6)*MATRIX<T,3>::Cross_Product_Matrix(B);
    data.H[1][2]=-(T)(1./6)*MATRIX<T,3>::Cross_Product_Matrix(A);
    for(int i=0;i<3;i++) for(int j=i+1;j<3;j++) data.H[j][i]=data.H[i][j].Transposed();
}

template<class T,class TV> void Intersect_Triangle_Point(PT_DATA<T>& data,const TV& A,const TV& B,const TV& C,const TV& P)
{
    data.n=4;
    TV n=TV::Cross_Product(B-A,C-A),U=TV::Cross_Product(B-C,P),V=TV::Cross_Product(C-A,P),W=TV::Cross_Product(A-B,P);
    T nP=TV::Dot_Product(n,P),nA=TV::Dot_Product(n,A),ABC=TV::Dot_Product(TV::Cross_Product(A,B),C);
    T BCP=TV::Dot_Product(TV::Cross_Product(B,C),P),CAP=TV::Dot_Product(TV::Cross_Product(C,A),P),ABP=TV::Dot_Product(TV::Cross_Product(A,B),P);
    data.V=(nA/nP)*P;
    MATRIX<T,3> onP=MATRIX<T,3>::Outer_Product(n,P),onU=MATRIX<T,3>::Outer_Product(n,U),onV=MATRIX<T,3>::Outer_Product(n,V),onW=MATRIX<T,3>::Outer_Product(n,W);
    MATRIX<T,3> onn=MATRIX<T,3>::Outer_Product(n,n);
    data.G[0]=BCP/(nP*nP)*onP;
    data.G[1]=CAP/(nP*nP)*onP;
    data.G[2]=ABP/(nP*nP)*onP;
    data.G[3]=-ABC/(nP*nP)*(onP-nP);

    // T a1=A.x;
    // T a2=A.y;
    // T a3=A.z;
    // T b1=B.x;
    // T b2=B.y;
    // T b3=B.z;
    // T c1=C.x;
    // T c2=C.y;
    // T c3=C.z;
    // T d1=P.x;
    // T d2=P.y;
    // T d3=P.z;
    // double cg[12][3];
    //for(int i=0;i<4;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) printf("g %f\n", data.G[i](j+1,k+1)-cg[i*3+j][k]);

    MATRIX<T,3> H00=-BCP/(nP*nP*nP)*(onU+onU.Transposed());
    data.H[0][0][0]=P.x*H00;
    data.H[1][0][0]=P.y*H00;
    data.H[2][0][0]=P.z*H00;
    MATRIX<T,3> H11=-CAP/(nP*nP*nP)*(onV+onV.Transposed());
    data.H[0][1][1]=P.x*H11;
    data.H[1][1][1]=P.y*H11;
    data.H[2][1][1]=P.z*H11;
    MATRIX<T,3> H22=-ABP/(nP*nP*nP)*(onW+onW.Transposed());
    data.H[0][2][2]=P.x*H22;
    data.H[1][2][2]=P.y*H22;
    data.H[2][2][2]=P.z*H22;
    MATRIX<T,3> H01=BCP*onW.Transposed()+BCP*onU.Transposed()-CAP*onU;
    data.H[0][0][1]=P.x*H01;
    data.H[1][0][1]=P.y*H01;
    data.H[2][0][1]=P.z*H01;
    MATRIX<T,3> H12=CAP*onU.Transposed()+CAP*onV.Transposed()-ABP*onV;
    data.H[0][1][2]=P.x*H12;
    data.H[1][1][2]=P.y*H12;
    data.H[2][1][2]=P.z*H12;
    MATRIX<T,3> H20=ABP*onV.Transposed()+ABP*onW.Transposed()-BCP*onW;
    data.H[0][2][0]=P.x*H20.Transposed();
    data.H[1][2][0]=P.y*H20.Transposed();
    data.H[2][2][0]=P.z*H20.Transposed();
    MATRIX<T,3> H33=(T)2*onn;
    for(int i=0;i<3;i++) data.H[i][3][3]=P(i+1)*H33-nP*(MATRIX<T,3>::Outer_Product(TV::Axis_Vector(i+1),n)+MATRIX<T,3>::Outer_Product(n,TV::Axis_Vector(i+1)));
    MATRIX<T,3> H03=ABC*onU-BCP*onn;
    for(int i=0;i<3;i++) data.H[i][0][3]=P(i+1)*H03+MATRIX<T,3>::Outer_Product(TV::Axis_Vector(i+1),n)*nP*BCP;
    MATRIX<T,3> H13=ABC*onV-CAP*onn;
    for(int i=0;i<3;i++) data.H[i][1][3]=P(i+1)*H13+MATRIX<T,3>::Outer_Product(TV::Axis_Vector(i+1),n)*nP*CAP;
    MATRIX<T,3> H23=ABC*onW-ABP*onn;
    for(int i=0;i<3;i++) data.H[i][2][3]=P(i+1)*H23+MATRIX<T,3>::Outer_Product(TV::Axis_Vector(i+1),n)*nP*ABP;
    for(int s=0;s<3;s++) for(int i=0;i<4;i++) for(int j=i+1;j<4;j++) data.H[s][j][i]=data.H[s][i][j].Transposed();
}

template<class T,class TV> void Intersect_Triangle_Segment(PT_DATA<T>& data,const TV& A,const TV& B,const TV& C,const TV& P,const TV& Q)
{
    data.n=5;
    PT_DATA<T> tdata;
    Intersect_Triangle_Point(tdata,A-Q,B-Q,C-Q,P-Q);
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
    data.n=4;
    PT_DATA<T> tdata;
    Intersect_Triangle_Point(tdata,-Q,A-Q,B-Q,P-Q);
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

template<class T> void Combine_Data(VOL_DATA<T,6>& data,const VOL_DATA<T,3>& V,const PT_DATA<T>& data_m,const PT_DATA<T>& data_n,const PT_DATA<T>& data_p)
{
    data.V+=V.V;
    const PT_DATA<T>* pd[3] = {&data_m, &data_n, &data_p};

    for(int z=0;z<3;z++) for(int j=0;j<pd[z]->n;j++) data.G[pd[z]->index[j]]+=V.G[z]*pd[z]->G[j];

    for(int z=0;z<3;z++)
        for(int j=0;j<pd[z]->n;j++)
            for(int s=0;s<pd[z]->n;s++)
                for(int i=0;i<3;i++)
                    data.H[pd[z]->index[j]][pd[z]->index[s]]+=V.G[z](i+1)*pd[z]->H[i][j][s];

    for(int y=0;y<3;y++)
        for(int z=0;z<3;z++)
            for(int j=0;j<pd[y]->n;j++)
                for(int s=0;s<pd[z]->n;s++)
                    data.H[pd[y]->index[j]][pd[z]->index[s]]+=pd[y]->G[j].Transpose_Times(V.H[y][z]*pd[z]->G[s]);
}

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
    return TV::Dot_Product(n,pts[0])/TV::Dot_Product(n,pts[opp_pt[b]+3])*pts[opp_pt[b]+3];
}

template<class T,class TV> void Data_From_Dof_Helper(PT_DATA<T>& data, int* i, const TV* pts)
{
    data.index[0]=i[0];
    printf("case (%i)\n", i[0]);
    Data_From_Dof(data, pts[i[0]]);
}

template<class T,class TV> void Intersect_Triangle_Point_Helper(PT_DATA<T>& data, int* i, const TV* pts)
{
    for(int j=0;j<4;j++) data.index[j]=i[j];
    printf("case (%i %i %i : %i)\n", i[0], i[1], i[2], i[3]);
    Intersect_Triangle_Point(data, pts[i[0]], pts[i[1]], pts[i[2]], pts[i[3]]);
}

template<class T,class TV> void Intersect_Triangle_Segment_Helper(PT_DATA<T>& data, int* i, const TV* pts)
{
    for(int j=0;j<5;j++) data.index[j]=i[j];
    printf("case (%i %i %i : %i %i)\n", i[0], i[1], i[2], i[3], i[4]);
    Intersect_Triangle_Segment(data, pts[i[0]], pts[i[1]], pts[i[2]], pts[i[3]], pts[i[4]]);
}

template<class T,class TV> void Intersect_Segment_Segment_Helper(PT_DATA<T>& data, int* i, const TV* pts)
{
    for(int j=0;j<4;j++) data.index[j]=i[j];
    printf("case (%i %i : %i %i)\n", i[0], i[1], i[2], i[3]);
    Intersect_Segment_Segment(data, pts[i[0]], pts[i[1]], pts[i[2]], pts[i[3]]);
}

template<class T,class TV> void Init_Data_From_Planes(void (*funcs[256])(PT_DATA<T>& data, int* ind, const TV* pts), int indices[256][5])
{
    for(int i=0;i<3;i++){
        int k=7-(1<<i);
        funcs[k|64]=&Data_From_Dof_Helper;
        indices[k|64][0]=i;
        funcs[k*8|128]=&Data_From_Dof_Helper;
        indices[k*8|128][0]=i+3;}

    for(int i=0;i<3;i++){
        int k=7-(1<<i);
        funcs[k|128]=&Intersect_Triangle_Point_Helper;
        indices[k|128][0]=3;
        indices[k|128][1]=4;
        indices[k|128][2]=5;
        indices[k|128][3]=i;
        funcs[k*8|64]=&Intersect_Triangle_Point_Helper;
        indices[k*8|64][0]=0;
        indices[k*8|64][1]=1;
        indices[k*8|64][2]=2;
        indices[k*8|64][3]=i+3;}

    for(int i=0;i<3;i++){
        int k=1<<i;
        funcs[k|192]=&Intersect_Triangle_Segment_Helper;
        indices[k|192][0]=3;
        indices[k|192][1]=4;
        indices[k|192][2]=5;
        indices[k|192][3]=(i+1)%3;
        indices[k|192][4]=(i+2)%3;
        funcs[k*8|192]=&Intersect_Triangle_Segment_Helper;
        indices[k*8|192][0]=0;
        indices[k*8|192][1]=1;
        indices[k*8|192][2]=2;
        indices[k*8|192][3]=(i+1)%3+3;
        indices[k*8|192][4]=(i+2)%3+3;}

    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++){
            int k=(1<<i)|(1<<(j+3));
            funcs[k|128]=&Intersect_Segment_Segment_Helper;
            indices[k|128][0]=(i+1)%3;
            indices[k|128][1]=(i+2)%3;
            indices[k|128][2]=(j+1)%3+3;
            indices[k|128][3]=(j+2)%3+3;
            funcs[k|64]=&Intersect_Segment_Segment_Helper;
            indices[k|64][0]=(j+1)%3+3;
            indices[k|64][1]=(j+2)%3+3;
            indices[k|64][2]=(i+1)%3;
            indices[k|64][3]=(i+2)%3;}
}

template<class T,class TV> void Volume_From_Tetrahedron(VOL_DATA<T,6>& data,TV pts[6],int va,int vb,int vc)
{
    static void (*funcs[256])(PT_DATA<T>& data, int* ind, const TV* pts);
    static int indices[256][5];
    static bool initialized=false;
    if(!initialized){initialized=true;Init_Data_From_Planes(funcs,indices);}

    printf("tet %i %i %i\n", va, vb, vc);
    PT_DATA<T> pd1,pd2,pd3;
    funcs[va](pd1,indices[va],pts);
    funcs[vb](pd2,indices[vb],pts);
    funcs[vc](pd3,indices[vc],pts);

    VOL_DATA<T,3> vd;
    Volume_From_Points(vd,pd1.V,pd2.V,pd3.V);

    Combine_Data(data,vd,pd1,pd2,pd3);
}

template<class T,class TV> void Volume_From_Triangles(VOL_DATA<T,6>& data,TV A,TV B,TV C,TV D,TV E,TV F)
{
    TV pts[6] = {A,B,C,D,E,F};
    Volume_From_Triangles_Cut(data,pts);
}

namespace{
    template<class TV>
    struct LIST
    {
        int planes;
        TV pt;
        bool inside;
    };
}

template<class T,class TV> void Volume_From_Triangles_Cut(VOL_DATA<T,6>& data,TV pts[6])
{
    int index[6]={0,1,2,3,4,5};
    T sign=1;
    if(TV::Triple_Product(pts[0],pts[1],pts[2])<0){exchange(pts[1],pts[2]);exchange(index[1],index[2]);sign=-sign;}
    if(TV::Triple_Product(pts[3],pts[4],pts[5])<0){exchange(pts[4],pts[5]);exchange(index[4],index[5]);sign=-sign;}
    int n=3;
    LIST<TV> alist[50], *list=alist;
    for(int i=0;i<n;i++){list[i].planes=7&~(1<<i);list[i].pt=pts[i];}
    VOL_DATA<T,6> tdata;

    // TODO: Robustness to in out in out.
    for(int p=3;p<6;p++){
        TV N=TV::Cross_Product(pts[(p+1)%3+3],pts[(p+2)%3+3]);
        bool has[2]={false,false};
        for(int i=0;i<n;i++){
            list[i].inside=TV::Dot_Product(N,list[i].pt)>0;
            has[list[i].inside]=true;}
        if(!has[1]){
            Clear(data);
            return;}
        if(!has[0]) continue;
        while(list[n-1].inside || !list[0].inside){list[n]=list[0];list++;}
        int f=1;
        while(list[f].inside) f++;
        int verta=(list[n-1].planes&list[0].planes)|(1<<p),vertb=(list[f-1].planes&list[f].planes)|(1<<p);
        list[f].planes=vertb;
        list[f+1].planes=verta;
        list[f].pt=Point_From_Planes(list[f].planes,pts);
        list[f+1].pt=Point_From_Planes(list[f+1].planes,pts);
        n=f+2;}

    TV N=TV::Cross_Product(pts[4]-pts[3],pts[5]-pts[3]);

    bool has[2]={false,false};
    for(int i=0;i<n;i++){
        list[i].inside=TV::Dot_Product(N,list[i].pt-pts[3])<0;
        has[list[i].inside]=true;}

    Clear(tdata);
    if(!has[0]) for(int f=2;f<n;f++) Volume_From_Tetrahedron(tdata,pts,list[0].planes|64,list[f-1].planes|64,list[f].planes|64);
    else if(!has[1]) for(int f=2;f<n;f++) Volume_From_Tetrahedron(tdata,pts,list[0].planes|128,list[f-1].planes|128,list[f].planes|128);
    else{
        while(list[n-1].inside || !list[0].inside){list[n]=list[0];list++;}
        int f=1;
        while(list[f].inside) f++;
        int verta=(list[n-1].planes&list[0].planes)|192,vertb=(list[f-1].planes&list[f].planes)|192;

        for(int i=1;i<f;i++) Volume_From_Tetrahedron(tdata,pts,verta|64,list[i-1].planes|64,list[i].planes|64);
        Volume_From_Tetrahedron(tdata,pts,verta|64,list[f-1].planes|64,vertb|64);

        for(int i=f+1;i<n;i++) Volume_From_Tetrahedron(tdata,pts,vertb|128,list[i-1].planes|128,list[i].planes|128);
        Volume_From_Tetrahedron(tdata,pts,vertb|128,list[n-1].planes|128,verta|128);}

    data.V=sign*tdata.V;
    for(int i=0;i<6;i++) data.G[index[i]]=sign*tdata.G[i];
    for(int i=0;i<6;i++) for(int k=0;k<6;k++) data.H[index[i]][index[k]]=sign*tdata.H[i][k];
}

template void Volume_From_Triangles<float,VECTOR<float,3> >(VOL_DATA<float,6>&,VECTOR<float,3>,VECTOR<float,3>,VECTOR<float,3>,VECTOR<float,3>,VECTOR<float,3>,VECTOR<float,3>);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void Volume_From_Triangles<double,VECTOR<double,3> >(VOL_DATA<double,6>&,VECTOR<double,3>,VECTOR<double,3>,VECTOR<double,3>,VECTOR<double,3>,VECTOR<double,3>,VECTOR<double,3>);
#endif
}
}

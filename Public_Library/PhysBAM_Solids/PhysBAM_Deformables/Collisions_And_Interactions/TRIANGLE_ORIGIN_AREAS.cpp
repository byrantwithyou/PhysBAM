#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_ORIGIN_AREAS.h>
#include <PhysBAM_Dynamics/Read_Write/EPS_FILE_GEOMETRY.h>

//extern PhysBAM::ARRAY<int> trap_cases;
using namespace PhysBAM;
using namespace ORIGIN_AREAS;

template<class T>
struct PT_DATA
{
    int n;
    int index[5];
    VECTOR<T,3> V;
    MATRIX<T,3> G[5];
    MATRIX<T,3> H[3][5][5];
};

template<class T,class TV> void Data_From_Dof(PT_DATA<T>& data,const TV& A)
{
    data.n=1;
    data.V=A;
    data.G[0]=MATRIX<T,3>::Identity_Matrix();
    for(int i=0;i<3;i++) data.H[i][0][0]=MATRIX<T,3>();
}

template<class T,class TV> void Volume_From_Points(VOL_DATA<T,3,3>& data,const TV& A,const TV& B,const TV& C)
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

template<class T,class TV> void Intersect_Triangle_Point(PT_DATA<T>& data,const TV A[3],const TV& P)
{
    data.n=4;
    TV n=TV::Cross_Product(A[1]-A[0],A[2]-A[0]),U[3];
    T nP=TV::Dot_Product(n,P),nA=TV::Dot_Product(n,A[0]),ABC=TV::Dot_Product(TV::Cross_Product(A[0],A[1]),A[2]),XYP[3];
    MATRIX<T,3> oPn=MATRIX<T,3>::Outer_Product(P,n),onn=MATRIX<T,3>::Outer_Product(n,n),onU[3],H33=(T)2/(nP*nP*nP)*ABC*onn;
    data.V=(nA/nP)*P;
    data.G[3]=-ABC/(nP*nP)*(oPn-nP);
    for(int k=0;k<3;k++){
        int r=(k+1)%3,s=(k+2)%3;
        U[k]=TV::Cross_Product(A[r]-A[s],P);
        XYP[k]=TV::Dot_Product(TV::Cross_Product(A[r],A[s]),P)/(nP*nP);
        data.G[k]=XYP[k]*oPn;
        XYP[k]/=nP;
        onU[k]=MATRIX<T,3>::Outer_Product(n,U[k]);
        data.H[k][3][3]=P(k)*H33-nP*ABC/(nP*nP*nP)*(MATRIX<T,3>::Outer_Product(TV::Axis_Vector(k),n)+MATRIX<T,3>::Outer_Product(n,TV::Axis_Vector(k)));}
    for(int k=0;k<3;k++){
        int r=(k+1)%3,s=(k+2)%3;
        MATRIX<T,3> Hkk=-XYP[k]*(onU[k]+onU[k].Transposed());
        MATRIX<T,3> Hrs=XYP[r]*onU[k].Transposed()+XYP[r]*onU[r].Transposed()-XYP[s]*onU[r];
        MATRIX<T,3> Hk3=ABC/(nP*nP*nP)*onU[k]-XYP[k]*onn;
        for(int i=0;i<3;i++){
            data.H[i][k][k]=P(i)*Hkk;
            data.H[i][r][s]=P(i)*Hrs;
            data.H[i][s][r]=data.H[i][r][s].Transposed();
            data.H[i][k][3]=P(i)*Hk3+MATRIX<T,3>::Outer_Product(n,TV::Axis_Vector(i))*nP*XYP[k];
            data.H[i][3][k]=data.H[i][k][3].Transposed();}}
}

template<class T,class TV> void Intersect_Triangle_Segment(PT_DATA<T>& data,const TV& A,const TV& B,const TV& C,const TV& P,const TV& Q)
{
    data.n=5;
    PT_DATA<T> tdata;
    TV pt3[3]={A-Q,B-Q,C-Q};
    Intersect_Triangle_Point(tdata,pt3,P-Q);
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
    TV pt3[3]={-Q,A-Q,B-Q};
    Intersect_Triangle_Point(tdata,pt3,P-Q);
    data.V=tdata.V+Q;
    data.G[3]=MATRIX<T,3>::Identity_Matrix();
    for(int i=0;i<4;i++){
        if(i) data.G[i]=tdata.G[i];
        data.G[3]-=tdata.G[i];}

    for(int s=0;s<3;s++) for(int i=0;i<4;i++){data.H[s][i][3]=data.H[s][3][i]=MATRIX<T,3>();}
    for(int s=0;s<3;s++) for(int i=0;i<4;i++) for(int j=0;j<4;j++){
        if(i && j) data.H[s][i][j]=tdata.H[s][i][j];
        if(i) data.H[s][i][3]-=tdata.H[s][i][j];
        if(j) data.H[s][3][j]-=tdata.H[s][i][j];
        data.H[s][3][3]+=tdata.H[s][i][j];}
}

template<class T> void Combine_Data(VOL_DATA<T,3,6>& data,const VOL_DATA<T,3,3>& V,const PT_DATA<T> pd[3])
{
    data.V+=V.V;

    for(int z=0;z<3;z++) for(int j=0;j<pd[z].n;j++) data.G[pd[z].index[j]]+=V.G[z]*pd[z].G[j];

    for(int z=0;z<3;z++)
        for(int j=0;j<pd[z].n;j++)
            for(int s=0;s<pd[z].n;s++)
                for(int i=0;i<3;i++)
                    data.H[pd[z].index[j]][pd[z].index[s]]+=V.G[z](i)*pd[z].H[i][j][s];

    for(int y=0;y<3;y++)
        for(int z=0;z<3;z++)
            for(int j=0;j<pd[y].n;j++)
                for(int s=0;s<pd[z].n;s++)
                    data.H[pd[y].index[j]][pd[z].index[s]]+=pd[y].G[j].Transpose_Times(V.H[y][z]*pd[z].G[s]);
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
//    printf("case (%i)\n", i[0]);
    Data_From_Dof(data, pts[i[0]]);
}

template<class T,class TV> void Intersect_Triangle_Point_Helper(PT_DATA<T>& data, int* i, const TV* pts)
{
    for(int j=0;j<4;j++) data.index[j]=i[j];
//    printf("case (%i %i %i : %i)\n", i[0], i[1], i[2], i[3]);
    TV pt3[3]={pts[i[0]], pts[i[1]], pts[i[2]]};
    Intersect_Triangle_Point(data, pt3, pts[i[3]]);
}

template<class T,class TV> void Intersect_Triangle_Segment_Helper(PT_DATA<T>& data, int* i, const TV* pts)
{
    for(int j=0;j<5;j++) data.index[j]=i[j];
//    printf("case (%i %i %i : %i %i)\n", i[0], i[1], i[2], i[3], i[4]);
    Intersect_Triangle_Segment(data, pts[i[0]], pts[i[1]], pts[i[2]], pts[i[3]], pts[i[4]]);
}

template<class T,class TV> void Intersect_Segment_Segment_Helper(PT_DATA<T>& data, int* i, const TV* pts)
{
    for(int j=0;j<4;j++) data.index[j]=i[j];
//    printf("case (%i %i : %i %i)\n", i[0], i[1], i[2], i[3]);
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

template<class T,class TV> void Volume_From_Tetrahedron(VOL_DATA<T,3,6>& data,TV pts[6],int va,int vb,int vc)
{
    static void (*funcs[256])(PT_DATA<T>& data, int* ind, const TV* pts);
    static int indices[256][5];
    static bool initialized=false;
    if(!initialized){initialized=true;Init_Data_From_Planes(funcs,indices);}

//    printf("tet %i %i %i\n", va, vb, vc);
    PT_DATA<T> pd[3];
    funcs[va](pd[0],indices[va],pts);
    funcs[vb](pd[1],indices[vb],pts);
    funcs[vc](pd[2],indices[vc],pts);

    VOL_DATA<T,3,3> vd;
    Volume_From_Points(vd,pd[0].V,pd[1].V,pd[2].V);

    Combine_Data(data,vd,pd);
}

template<class T,class TV> void Volume_From_Simplices_Helper(VOL_DATA<T,3,6>& data,TV const (&LA)[6])
{
    TV pts[]={LA[0],LA[1],LA[2],LA[3],LA[4],LA[5]};
    struct LIST
    {
        int planes;
        TV pt;
        bool inside;
    };
    int index[6]={0,1,2,3,4,5};
    T sign=1;
    if(TV::Triple_Product(pts[0],pts[1],pts[2])<0){exchange(pts[1],pts[2]);exchange(index[1],index[2]);sign=-sign;}
    if(TV::Triple_Product(pts[3],pts[4],pts[5])<0){exchange(pts[4],pts[5]);exchange(index[4],index[5]);sign=-sign;}
    int n=3;
    LIST alist[50], *list=alist;
    for(int i=0;i<n;i++){list[i].planes=7&~(1<<i);list[i].pt=pts[i];}
    VOL_DATA<T,3,6> tdata;

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

template<class T,class TV> void PhysBAM::ORIGIN_AREAS::Volume_From_Simplices(VOL_DATA<T,3,6>& data,const TV& X0,TV const A[6])
{
    TV const B[]={A[0]-X0,A[1]-X0,A[2]-X0,A[3]-X0,A[4]-X0,A[5]-X0};
    Volume_From_Simplices_Helper(data,B);
}

template void PhysBAM::ORIGIN_AREAS::Volume_From_Simplices<float,VECTOR<float,3> >(VOL_DATA<float,3,6>&,VECTOR<float,3> const &,VECTOR<float,3> const [6]);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void PhysBAM::ORIGIN_AREAS::Volume_From_Simplices<double,VECTOR<double,3> >(VOL_DATA<double,3,6>&,VECTOR<double,3> const &,VECTOR<double,3> const [6]);
#endif

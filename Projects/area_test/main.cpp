//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/SEGMENT_ORIGIN_AREAS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRAPEZOID_INTERSECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_INTERSECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/VOLUME_COLLISIONS.h>
#include <PhysBAM_Dynamics/Read_Write/EPS_FILE_GEOMETRY.h>
#include <cassert>
#include <cmath>
#include <iomanip>
using namespace PhysBAM;

ARRAY<int> trap_cases;

template<class TV>
bool Inside_Trapezoid(TV a,TV b,const TV& p)
{
    if(a.x>b.x) exchange(a,b);
    PHYSBAM_ASSERT(a.y>=0);
    PHYSBAM_ASSERT(b.y>=0);
    PHYSBAM_ASSERT(p.y>=0);

    if(p.x<a.x) return false;
    if(p.x>b.x) return false;
    if((-a.y*b.x-b.y*p.x+a.x*b.y+a.y*p.x+p.y*b.x-p.y*a.x)>0) return false;
    return true;
}

RANDOM_NUMBERS<double> rn;

bool Test()
{
    typedef double T;
    typedef VECTOR<double,2> TV;
    VECTOR<TV,4> a;
    rn.Fill_Uniform(a(0),0,1);
    rn.Fill_Uniform(a(1),0,1);
    rn.Fill_Uniform(a(2),0,1);
    rn.Fill_Uniform(a(3),0,1);

    T e=1e-5;
    VECTOR<TV,4> da;
    rn.Fill_Uniform(da(0),-e,e);
    rn.Fill_Uniform(da(1),-e,e);
    rn.Fill_Uniform(da(2),-e,e);
    rn.Fill_Uniform(da(3),-e,e);

    trap_cases.Remove_All();
//    fprintf(stderr, "1 1 1 setrgbcolor newpath 0 0 moveto 1000 0 lineto 1000 1000 lineto 0 1000 lineto closepath fill\n");
    VECTOR<TV,4> G1;
    VECTOR<VECTOR<MATRIX<T,2>,4>,4> H1;
    T A1 = Trapezoid_Intersection_Area(a(0),a(1),a(2),a(3),G1,H1);
    ARRAY<int> tmp_cases=trap_cases;

    trap_cases.Remove_All();
    VECTOR<TV,4> G2;
    VECTOR<VECTOR<MATRIX<T,2>,4>,4> H2;
    T A2 = Trapezoid_Intersection_Area(a(0)+da(0),a(1)+da(1),a(2)+da(2),a(3)+da(3),G2,H2);

    VECTOR<T,8>& Va=(VECTOR<T,8>&)da;
    VECTOR<T,8>& V1=(VECTOR<T,8>&)G1;
    VECTOR<T,8>& V2=(VECTOR<T,8>&)G2;

    T G=Va.Dot_Product(Va,(V1+V2)/(T)2);
    T aa=G/e;
    T bb=(A2-A1)/e;
    if((!aa != !bb) || fabs((aa-bb)/bb)>1e-5) if(tmp_cases==trap_cases){printf("ZG %g %g %g   ", aa, bb, fabs((aa-bb)/bb));LOG::cout<<tmp_cases<<"   "<<trap_cases<<std::endl;}

    MATRIX<T,8> M1,M2;
    for(int i=0;i<8;i++) for(int j=0;j<8;j++) M1(i,j)=H1(i/2)(j/2)(i%2,j%2);
    for(int i=0;i<8;i++) for(int j=0;j<8;j++) M2(i,j)=H2(i/2)(j/2)(i%2,j%2);

    VECTOR<T,8> dG1=(M1+M2)/(T)2*Va;
    VECTOR<T,8> dG2=V2-V1;
    T cc=dG2.Magnitude()/e;
    T dd=dG1.Magnitude()/e;
    T ee=(dG2-dG1).Magnitude()/e;
    if((!cc != !dd) || ee/cc>1e-5) if(tmp_cases==trap_cases){printf("ZH %g %g %g      ", cc, dd, ee/cc);LOG::cout<<tmp_cases<<"   "<<trap_cases<<std::endl;}

    return true;
}

bool Tri_Test()
{
    typedef double T;
    typedef VECTOR<double,2> TV;
    VECTOR<TV,6> a;
    rn.Fill_Uniform(a(0),0,1);
    rn.Fill_Uniform(a(1),0,1);
    rn.Fill_Uniform(a(2),0,1);
    rn.Fill_Uniform(a(3),0,1);
    rn.Fill_Uniform(a(4),0,1);
    rn.Fill_Uniform(a(5),0,1);

    T e=1e-5;
    VECTOR<TV,6> da;
    rn.Fill_Uniform(da(0),-e,e);
    rn.Fill_Uniform(da(1),-e,e);
    rn.Fill_Uniform(da(2),-e,e);
    rn.Fill_Uniform(da(3),-e,e);
    rn.Fill_Uniform(da(4),-e,e);
    rn.Fill_Uniform(da(5),-e,e);

    trap_cases.Remove_All();
//    fprintf(stderr, "1 1 1 setrgbcolor newpath 0 0 moveto 1000 0 lineto 1000 1000 lineto 0 1000 lineto closepath fill\n");
    VECTOR<TV,6> G1;
    VECTOR<VECTOR<MATRIX<T,2>,6>,6> H1;
    T A1 = Triangle_Intersection_Area(TRIANGLE_2D<T>(a(0),a(1),a(2)),TRIANGLE_2D<T>(a(3),a(4),a(5)),G1,H1);
    ARRAY<int> tmp_cases=trap_cases;

    trap_cases.Remove_All();
    VECTOR<TV,6> G2;
    VECTOR<VECTOR<MATRIX<T,2>,6>,6> H2;
    T A2 = Triangle_Intersection_Area(TRIANGLE_2D<T>(a(0)+da(0),a(1)+da(1),a(2)+da(2)),TRIANGLE_2D<T>(a(3)+da(3),a(4)+da(4),a(5)+da(5)),G2,H2);

    VECTOR<T,12>& Va=(VECTOR<T,12>&)da;
    VECTOR<T,12>& V1=(VECTOR<T,12>&)G1;
    VECTOR<T,12>& V2=(VECTOR<T,12>&)G2;

    T G=Va.Dot_Product(Va,(V1+V2)/(T)2);
    T aa=G/e;
    T bb=(A2-A1)/e;
    if((fabs(aa)>1e-5 || fabs(bb)>1e-5) && fabs((aa-bb)/bb)>1e-5){printf("TG %g %g %g   ", aa, bb, fabs((aa-bb)/bb));LOG::cout<<tmp_cases<<"   "<<trap_cases<<std::endl;}

    MATRIX<T,12> M1,M2;
    for(int i=0;i<12;i++) for(int j=0;j<12;j++) M1(i,j)=H1(i/2)(j/2)(i%2,j%2);
    for(int i=0;i<12;i++) for(int j=0;j<12;j++) M2(i,j)=H2(i/2)(j/2)(i%2,j%2);

    VECTOR<T,12> dG1=(M1+M2)/(T)2*Va/e;
    VECTOR<T,12> dG2=(V2-V1)/e;
    T cc=dG2.Magnitude();
    T dd=dG2.Magnitude();
    T ee=(dG2-dG1).Magnitude();
    if((fabs(cc)>1e-5 || fabs(dd)>1e-5) && ee/dd>1e-5){printf("TH %g %g %g   ", cc, dd, ee/dd);LOG::cout<<tmp_cases<<"   "<<trap_cases<<std::endl;}

    return true;
}

void Test_Triangulated_Areas()
{
    typedef double T;
    typedef VECTOR<double,2> TV;
    TRIANGULATED_AREA<T>* ta1=TESSELLATION::Generate_Triangles(SPHERE<TV>(TV(-(T).41,0),1.12),3);
    TRIANGULATED_AREA<T>* ta2=TESSELLATION::Generate_Triangles(SPHERE<TV>(TV((T).5,(T).02),1),3);

    ta2->mesh.elements.Flattened()+=ta1->particles.X.m;
    ta1->particles.array_collection->Append(*ta2->particles.array_collection);
    ta1->mesh.elements.Append_Elements(ta2->mesh.elements);
    ta1->Update_Number_Nodes();
    delete ta2;
    ta2=0;

    VOLUME_COLLISIONS<TV> vc;
    vc.objects.Append(ta1);

    trap_cases.Remove_All();
    vc.Compute_Collision_Triangles();

    T A1=vc.area;
    VECTOR_ND<T> G1(ta1->particles.X.m*2);
    MATRIX_MXN<T> H1(ta1->particles.X.m*2,ta1->particles.X.m*2);
    for(HASHTABLE<int,TV>::ITERATOR it(vc.gradient);it.Valid();it.Next())
        for(int s=0;s<2;s++)
            G1(it.Key()*2-2+s)=it.Data()(s);
    for(HASHTABLE<VECTOR<int,2>,MATRIX<T,2> >::ITERATOR it(vc.hessian);it.Valid();it.Next())
        for(int s=0;s<2;s++)
            for(int t=0;t<2;t++)
                H1(it.Key().x*2-2+s,it.Key().y*2-2+t)=it.Data()(s,t);

    T e=1e-5;
    VECTOR_ND<T> vdX(ta1->particles.X.m*2);
    rn.Fill_Uniform(vdX,-e,e);
    ARRAY_VIEW<TV> dX(ta1->particles.X.m,(TV*)vdX.Get_Array_Pointer());
    ta1->particles.X+=dX;

    ARRAY<int> tmp_cases=trap_cases;
    trap_cases.Remove_All();
    vc.Compute_Collision_Triangles();

    T A2=vc.area;
    VECTOR_ND<T> G2(ta1->particles.X.m*2);
    MATRIX_MXN<T> H2(ta1->particles.X.m*2,ta1->particles.X.m*2);
    for(HASHTABLE<int,TV>::ITERATOR it(vc.gradient);it.Valid();it.Next())
        for(int s=0;s<2;s++)
            G2(it.Key()*2-2+s)=it.Data()(s);
    for(HASHTABLE<VECTOR<int,2>,MATRIX<T,2> >::ITERATOR it(vc.hessian);it.Valid();it.Next())
        for(int s=0;s<2;s++)
            for(int t=0;t<2;t++)
                H2(it.Key().x*2-2+s,it.Key().y*2-2+t)=it.Data()(s,t);

    OCTAVE_OUTPUT<T> oo("m.txt");
    oo.Write("A1",A1);
    oo.Write("A2",A2);
    oo.Write("G1",G1);
    oo.Write("G2",G2);
    oo.Write("H1",H1);
    oo.Write("H2",H2);

    T dE1=vdX.Dot_Product(vdX,G1+G2)/(2*e);
    T dE2=(A2-A1)/e;
    LOG::cout<<tmp_cases<<std::endl;
    LOG::cout<<trap_cases<<std::endl;
    printf("AG %g %g %g\n", dE1, dE2, fabs((dE1-dE2)/dE2));

    VECTOR_ND<T> dG1((H1*vdX+H2*vdX)/(T)2);
    VECTOR_ND<T> dG2(G2-G1);
    T dG1m=dG1.Magnitude()/e;
    T dG2m=dG2.Magnitude()/e;
    T DGm=(dG2-dG1).Magnitude()/e;
    printf("AH %g %g %g\n", dG1m, dG2m, fabs(DGm/dG2m));
}

int fail_number=0;

template<class TV>
void Test_Triangle_Intersection(const VECTOR<int,3>& t1,const VECTOR<int,3>& t2,const ARRAY<TV>& X)
{
    typedef typename TV::SCALAR T;
    T size=500;
    VECTOR<TV,6> G;
    VECTOR<VECTOR<MATRIX<T,2>,6>,6> H;
    T A=Triangle_Intersection_Area(TRIANGLE_2D<T>(X.Subset(t1)),TRIANGLE_2D<T>(X.Subset(t2)),G,H);
    bool B=Topology_Aware_Triangle_Intersection_Test(t1,t2,ARRAY_VIEW<const TV>(X));

    if((fabs(A)>1e-10)==B) return;

    printf("%g %i (fail %d)\n", A, B, fail_number);

    char buff[100];
    sprintf(buff, "fail-%d.eps", fail_number++);
    FILE * F = fopen(buff, "w");

    fprintf(F, "%%!PS-Adobe-3.0 EPSF-3.0\n");
    fprintf(F, "%%%%BoundingBox: 0 0 500 500\n");
    fprintf(F, "1 1 1 setrgbcolor newpath 0 0 moveto 500 0 lineto 500 500 lineto 0 500 lineto closepath fill\n");
    fprintf(F, "1 0 0 setrgbcolor %g %g moveto %g %g lineto %g %g lineto closepath stroke\n", X(t1.x).x*size, X(t1.x).y*size, X(t1.y).x*size, X(t1.y).y*size,X(t1.z).x*size, X(t1.z).y*size);
    fprintf(F, "0 1 0 setrgbcolor %g %g moveto %g %g lineto %g %g lineto closepath stroke\n", X(t2.x).x*size, X(t2.x).y*size, X(t2.y).x*size, X(t2.y).y*size,X(t2.z).x*size, X(t2.z).y*size);

    fclose(F);
}

template<class TV>
void Test_Triangle_Intersection()
{
    ARRAY<TV> X;
    X.Append(rn.Get_Uniform_Vector(RANGE<TV>::Unit_Box()));
    X.Append(rn.Get_Uniform_Vector(RANGE<TV>::Unit_Box()));
    X.Append(rn.Get_Uniform_Vector(RANGE<TV>::Unit_Box()));
    X.Append(rn.Get_Uniform_Vector(RANGE<TV>::Unit_Box()));
    X.Append(rn.Get_Uniform_Vector(RANGE<TV>::Unit_Box()));
    X.Append(rn.Get_Uniform_Vector(RANGE<TV>::Unit_Box()));
    VECTOR<int,3> t1(0,1,2);
    VECTOR<int,3> t2(2,1,3);
    VECTOR<int,3> t3(2,3,4);
    VECTOR<int,3> t4(3,4,5);

    Test_Triangle_Intersection(t1,t2,X);
    Test_Triangle_Intersection(t1,t3,X);
    Test_Triangle_Intersection(t1,t4,X);
}

void Case_Test()
{
    typedef double T;
    typedef VECTOR<double,2> TV;

    TV a,b,c,d;
    rn.Fill_Uniform(a,-1,1);
    rn.Fill_Uniform(b,-1,1);
    rn.Fill_Uniform(c,-1,1);
    rn.Fill_Uniform(d,-1,1);

    {
        EPS_FILE_GEOMETRY<T> eps("case.eps");
        eps.Line_Color(VECTOR<T,3>(.5,.5,.5));
        eps.Draw_Point(TV());
        eps.Draw_Point(TV(-1,0));
        eps.Draw_Point(TV(1,0));
        eps.Draw_Point(TV(0,1));
        eps.Draw_Point(TV(0,-1));
        eps.Line_Color(VECTOR<T,3>(1,0,0));
        eps.Draw_Line(a,b);
        eps.Draw_Point(a);
        eps.Draw_Line(TV(),a/a.Max_Abs());
        eps.Draw_Line(TV(),b/b.Max_Abs());
        eps.Line_Color(VECTOR<T,3>(0,1,0));
        eps.Draw_Line(c,d);
        eps.Draw_Point(c);
        eps.Draw_Line(TV(),c/c.Max_Abs());
        eps.Draw_Line(TV(),d/d.Max_Abs());
    }

    trap_cases.Remove_All();
    ORIGIN_AREAS::VOL_DATA<T,2,4> data;
    TV ar[4]={a,b,c,d};
    ORIGIN_AREAS::Volume_From_Simplices(data,TV(),ar);

    VECTOR<TV,6> G1;
    VECTOR<VECTOR<MATRIX<T,2>,6>,6> H1;
    T A1 = Triangle_Intersection_Area(TRIANGLE_2D<T>(TV(),a,b),TRIANGLE_2D<T>(TV(),c,d),G1,H1);

    int ii[13]={0,1,1,1,2,3,4,5,5,5,6,7,8};
    VECTOR<T,12>& V12=(VECTOR<T,12>&)G1;
    VECTOR<T,8> V;
    for(int i=0;i<12;i++) V(ii[i])=V12(i);
    MATRIX<T,8> M;
    for(int i=0;i<12;i++) for(int j=0;j<12;j++) M(ii[i],ii[j])=H1(i/2)(j/2)(i%2,j%2);
    VECTOR<T,8>& W=(VECTOR<T,8>&)data.G;
    MATRIX<T,8> N;for(int i=0;i<4;i++) for(int j=0;j<4;j++) N.Set_Submatrix(2*i,2*j,data.H[i][j]);

    printf("ERRORS: (case %i)  %.4f (%.4f)  %.4f (%.4f)  %.4f (%.4f)\n", (trap_cases.m?trap_cases(0):0), fabs(data.V-A1), fabs(A1), (W-V).Magnitude(), V.Magnitude(),
        (M-N).Frobenius_Norm(),M.Frobenius_Norm());
    LOG::cout<<"Areas:  "<<data.V<<"   "<<A1<<std::endl;
    LOG::cout<<"Gradients: "<<W<<"    "<<V<<std::endl;
    LOG::cout<<"Hessians: "<<N<<"    "<<M<<std::endl;
}

template<class T>
void Volume_From_Simplices_Test_Gradient(
    const PhysBAM::ORIGIN_AREAS::VOL_DATA<T,2,4>& vol_data,
    PhysBAM::VECTOR<T,2> (&X)[4])
{
    typedef PhysBAM::VECTOR<T,2> TV;
    const T tol=(T)1/(1024*1024);
    const T dx=(T)1/1024;
    PhysBAM::ORIGIN_AREAS::VOL_DATA<T,2,4> vol_data_dx0;
    PhysBAM::ORIGIN_AREAS::VOL_DATA<T,2,4> vol_data_dx1;
    for(int i=0;i!=4;++i){
        for(int d=0;d<2;d++){
            X[i](d)+=dx;
            PhysBAM::ORIGIN_AREAS::Clear(vol_data_dx1);
            PhysBAM::ORIGIN_AREAS::Volume_From_Simplices(vol_data_dx1,TV(),X);
            X[i](d)-=dx;
            X[i](d)-=dx;
            PhysBAM::ORIGIN_AREAS::Clear(vol_data_dx0);
            PhysBAM::ORIGIN_AREAS::Volume_From_Simplices(vol_data_dx0,TV(),X);
            X[i](d)+=dx;
            const T Gida=vol_data.G[i](d);
            const T Gidb=(vol_data_dx1.V-vol_data_dx0.V)/(2*dx);
            PHYSBAM_ASSERT(std::abs(Gida-Gidb)<=tol);
        }
    }
}

template<class T>
void Volume_From_Simplices_Test_Hessian(
    const PhysBAM::ORIGIN_AREAS::VOL_DATA<T,2,4>& vol_data,
    PhysBAM::VECTOR<T,2> (&X)[4])
{
    typedef PhysBAM::VECTOR<T,2> TV;
    const T tol=(T)1/(1024*1024);
    const T dx=(T)1/1024;
    PhysBAM::ORIGIN_AREAS::VOL_DATA<T,2,4> vol_data_dx00;
    PhysBAM::ORIGIN_AREAS::VOL_DATA<T,2,4> vol_data_dx01;
    PhysBAM::ORIGIN_AREAS::VOL_DATA<T,2,4> vol_data_dx10;
    PhysBAM::ORIGIN_AREAS::VOL_DATA<T,2,4> vol_data_dx11;
    for(int i=0;i!=4;++i){
        for(int j=0;j!=4;++j){
            for(int d=0;d<2;d++){
                for(int e=0;e<2;e++){
                    const T Hijdea=vol_data.H[i][j](d,e);
                    T Hijdeb;
                    if(i==j&&d==e){
                        X[i](d)+=dx;
                        PhysBAM::ORIGIN_AREAS::Clear(vol_data_dx11);
                        PhysBAM::ORIGIN_AREAS::Volume_From_Simplices(vol_data_dx11,TV(),X);
                        X[i](d)-=dx;
                        X[i](d)-=dx;
                        PhysBAM::ORIGIN_AREAS::Clear(vol_data_dx00);
                        PhysBAM::ORIGIN_AREAS::Volume_From_Simplices(vol_data_dx00,TV(),X);
                        X[i](d)+=dx;
                        Hijdeb=((vol_data_dx11.V-vol_data.V)-(vol_data.V-vol_data_dx00.V))/(dx*dx);
                    }
                    else{
                        X[i](d)+=dx;
                        X[j](e)+=dx;
                        PhysBAM::ORIGIN_AREAS::Clear(vol_data_dx11);
                        PhysBAM::ORIGIN_AREAS::Volume_From_Simplices(vol_data_dx11,TV(),X);
                        X[j](e)-=dx;
                        X[j](e)-=dx;
                        PhysBAM::ORIGIN_AREAS::Clear(vol_data_dx10);
                        PhysBAM::ORIGIN_AREAS::Volume_From_Simplices(vol_data_dx10,TV(),X);
                        X[j](e)+=dx;
                        X[i](d)-=dx;
                        X[i](d)-=dx;
                        X[j](e)+=dx;
                        PhysBAM::ORIGIN_AREAS::Clear(vol_data_dx01);
                        PhysBAM::ORIGIN_AREAS::Volume_From_Simplices(vol_data_dx01,TV(),X);
                        X[j](e)-=dx;
                        X[j](e)-=dx;
                        PhysBAM::ORIGIN_AREAS::Clear(vol_data_dx00);
                        PhysBAM::ORIGIN_AREAS::Volume_From_Simplices(vol_data_dx00,TV(),X);
                        X[j](e)+=dx;
                        X[i](d)+=dx;
                        Hijdeb=((vol_data_dx11.V-vol_data_dx10.V)-(vol_data_dx01.V-vol_data_dx00.V))/(4*dx*dx);
                    }
                    PHYSBAM_ASSERT(std::abs(Hijdea-Hijdeb)<=tol);
                }
            }
        }
    }
}

template<class T>
void Volume_From_Simplices_Test()
{
    typedef PhysBAM::VECTOR<T,2> TV;
    PhysBAM::ORIGIN_AREAS::VOL_DATA<T,2,4> vol_data;
    TV X[4];

    // Case "CCAA"
    // A X   X   X B
    //     D   C
    //       O
    PhysBAM::ORIGIN_AREAS::Clear(vol_data);
    X[0]=TV(-3,+2);X[1]=TV(+3,+2);X[2]=TV(+1,+1);X[3]=TV(-1,+1);
    PhysBAM::ORIGIN_AREAS::Volume_From_Simplices(vol_data,TV(),X);
    PHYSBAM_ASSERT(vol_data.V==static_cast<T>(-1));
    Volume_From_Simplices_Test_Gradient(vol_data,X);
    Volume_From_Simplices_Test_Hessian(vol_data,X);

    // Case "CCAB"
    // D
    //   X
    // A P Q   X B
    //       C
    //     O
    PhysBAM::ORIGIN_AREAS::Clear(vol_data);
    X[0]=TV(-2,+2);X[1]=TV(+3,+2);X[2]=TV(+1,+1);X[3]=TV(-2,+4);
    PhysBAM::ORIGIN_AREAS::Volume_From_Simplices(vol_data,TV(),X);
    PHYSBAM_ASSERT(vol_data.V==static_cast<T>(-2));
    Volume_From_Simplices_Test_Gradient(vol_data,X);
    Volume_From_Simplices_Test_Hessian(vol_data,X);

    // Case "CCBB"
    //   D       C
    // A   P   P   B
    //       O
    PhysBAM::ORIGIN_AREAS::Clear(vol_data);
    X[0]=TV(-3,+1);X[1]=TV(+3,+1);X[2]=TV(+2,+2);X[3]=TV(-2,+2);
    PhysBAM::ORIGIN_AREAS::Volume_From_Simplices(vol_data,TV(),X);
    PHYSBAM_ASSERT(vol_data.V==static_cast<T>(-1));
    Volume_From_Simplices_Test_Gradient(vol_data,X);
    Volume_From_Simplices_Test_Hessian(vol_data,X);

    // Case "BCAC"
    //   A         B
    // D   P   C
    //       O
    PhysBAM::ORIGIN_AREAS::Clear(vol_data);
    X[0]=TV(-2,+2);X[1]=TV(+3,+2);X[2]=TV(+1,+1);X[3]=TV(-3,+1);
    PhysBAM::ORIGIN_AREAS::Volume_From_Simplices(vol_data,TV(),X);
    PHYSBAM_ASSERT(vol_data.V==static_cast<T>(-1));
    Volume_From_Simplices_Test_Gradient(vol_data,X);
    Volume_From_Simplices_Test_Hessian(vol_data,X);

    // Case "BCBC"
    // A
    //   X
    // D P Q   C
    //       P
    //     O   B
    PhysBAM::ORIGIN_AREAS::Clear(vol_data);
    X[0]=TV(-2,+4);X[1]=TV(+2, 0);X[2]=TV(+2,+2);X[3]=TV(-2,+2);
    PhysBAM::ORIGIN_AREAS::Volume_From_Simplices(vol_data,TV(),X);
    PHYSBAM_ASSERT(vol_data.V==static_cast<T>(-2));
    Volume_From_Simplices_Test_Gradient(vol_data,X);
    Volume_From_Simplices_Test_Hessian(vol_data,X);

    // Case "ACAC"
    // D
    //   A  Q  X  X  X  X  B
    //            C
    //      O
    PhysBAM::ORIGIN_AREAS::Clear(vol_data);
    X[0]=TV(-1,+2);X[1]=TV(+5,+2);X[2]=TV(+2,+1);X[3]=TV(-2,+3);
    PhysBAM::ORIGIN_AREAS::Volume_From_Simplices(vol_data,TV(),X);
    PHYSBAM_ASSERT(vol_data.V==static_cast<T>(-3));
    Volume_From_Simplices_Test_Gradient(vol_data,X);
    Volume_From_Simplices_Test_Hessian(vol_data,X);
}

int main(int argc,char *argv[])
{
    typedef double T;
    typedef double RW;
    typedef VECTOR<T,2> TV;
    LOG::cout<<std::setprecision(16);

#if 0

//    Test_Triangulated_Areas();

//    fprintf(stderr, "%%!PS-Adobe-3.0 EPSF-3.0\n");
//    fprintf(stderr, "%%%%BoundingBox: 0 0 1000 1000\n");

    for(int k=0;k<100000;k++)
        Tri_Test();

//    for(int k=0;k<1000000;k++)
//        Test();

//    for(int i=0;i<100;i++) Test_Triangle_Intersection<TV>();

    for(int i=0;i<10000;i++)
    {
        try{Case_Test();}catch(...){}
    }

#else // #if 0|1

    Volume_From_Simplices_Test<T>();

#endif // #if 0|1

    return 0;
}

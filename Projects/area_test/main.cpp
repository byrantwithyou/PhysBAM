//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRAPEZOID_INTERSECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_INTERSECTION.h>
using namespace PhysBAM;

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
    rn.Fill_Uniform(a(1),0,1);
    rn.Fill_Uniform(a(2),0,1);
    rn.Fill_Uniform(a(3),0,1);
    rn.Fill_Uniform(a(4),0,1);

    T e=1e-4;
    VECTOR<TV,4> da;
    rn.Fill_Uniform(da(1),-e,e);
    rn.Fill_Uniform(da(2),-e,e);
    rn.Fill_Uniform(da(3),-e,e);
    rn.Fill_Uniform(da(4),-e,e);

    fprintf(stderr, "1 1 1 setrgbcolor newpath 0 0 moveto 1000 0 lineto 1000 1000 lineto 0 1000 lineto closepath fill\n");
    VECTOR<TV,4> dA1;
    VECTOR<VECTOR<MATRIX<T,2>,4>,4> H1;
    T A1 = Trapezoid_Intersection_Area(a(1),a(2),a(3),a(4),dA1,H1);

    VECTOR<TV,4> dA2;
    VECTOR<VECTOR<MATRIX<T,2>,4>,4> H2;
    T A2 = Trapezoid_Intersection_Area(a(1)+da(1),a(2)+da(2),a(3)+da(3),a(4)+da(4),dA2,H2);

    VECTOR<T,8>& Va=(VECTOR<T,8>&)da;
    VECTOR<T,8>& V1=(VECTOR<T,8>&)dA1;
    VECTOR<T,8>& V2=(VECTOR<T,8>&)dA2;

    T dA=Va.Dot_Product(Va,(V1+V2)/(T)2);
    T aa=dA/e;
    T bb=(A2-A1)/e;
    if((!aa != !bb) || fabs((aa-bb)/bb)>1e-5) printf("AG %g %g %g\n", aa, bb, fabs((aa-bb)/bb));

    MATRIX<T,8> M1,M2;
    for(int i=1;i<=8;i++) for(int j=1;j<=8;j++) M1(i,j)=H1((i+1)/2)((j+1)/2)((i+1)%2+1,(j+1)%2+1);
    for(int i=1;i<=8;i++) for(int j=1;j<=8;j++) M2(i,j)=H2((i+1)/2)((j+1)/2)((i+1)%2+1,(j+1)%2+1);

    VECTOR<T,8> dG1=(M1+M2)/(T)2*Va;
    VECTOR<T,8> dG2=V2-V1;
    T cc=dG2.Magnitude()/e;
    T dd=(dG2-dG1).Magnitude()/e;
    if((!cc != !dd) || dd/cc>1e-5) printf("AH %g %g %g       ", cc, dd, dd/cc),LOG::cout<<dG1<<" "<<dG2<<std::endl;

    return true;
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
    VECTOR<int,3> t1(1,2,3);
    VECTOR<int,3> t2(3,2,4);
    VECTOR<int,3> t3(3,4,5);
    VECTOR<int,3> t4(4,5,6);

    Test_Triangle_Intersection(t1,t2,X);
    Test_Triangle_Intersection(t1,t3,X);
    Test_Triangle_Intersection(t1,t4,X);
}

int main(int argc,char *argv[])
{
    typedef double T;
    typedef double RW;
    typedef VECTOR<T,2> TV;


//    fprintf(stderr, "%%!PS-Adobe-3.0 EPSF-3.0\n");
//    fprintf(stderr, "%%%%BoundingBox: 0 0 1000 1000\n");


//    for(int k=0;k<10000;k++){
//        Test();
//    }

    Test_Triangle_Intersection<TV>();

    return 0;
}
